/**
 * @file OrthoMAP_toolbox.cpp
 * @brief Core utility functions for the OrthoMap pipeline
 * 
 * High-performance C++ interface for R via Rcpp. This script defines a 
 * collection of utility functions used across the OrthoMAP pipeline. 
 * While the pipeline’s behavior can vary depending on the input data, 
 * these functions represent the invariant logic — stable, reusable 
 * components that support consistent downstream analysis across all 
 * input configurations. Use this toolbox to avoid duplication, promote 
 * clarity, and maintain consistent processing once the input-specific 
 * handling is complete.
 * 
 * Intended as a sourceable script: Rcpp::sourceCpp("OrthoMAP_toolbox.cpp")
 *
 * @author David Blickenstorfer, Technau Group (2025)
 * @date   25/06/2025
 * 
 * Version: v1.0.0
 * Changelog:
 *  - v1.0.0: Initial release
 * 
 * (c) Technau Group, University of Vienna, 2025
 * (c) Boeva Lab, ETH Zurich, 2025
 */

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#include <omp.h>    // OMP shared memory multi-threading
#include "OrthoMAP_toolbox.h"

// =========================================================
// STEP 1:  Create OrthoGroup against Celltype matrix      #
// Functions                                               #
//      convert_gene_to_og_expression                      #
// =========================================================

// Convert gene expression matrix to orthogroup matrix
//
// Converts a sparse gene-level expression matrix (`dgRMatrix`) into an 
// orthogroup-level expression matrix. For each orthogroup, the maximum 
// expression value of its member genes is used per cell.
//
// @param gene_expression_matrix A \code{dgRMatrix} object (sparse), containing
//        gene expression values (genes × cells).
// @param genes_names A character vector of gene names corresponding to the rows 
//        of \code{gene_expression_matrix}.
// @param orthogroup_names A character vector of all orthogroup names.
// @param lookup_gene_to_orthogroup A named character vector mapping gene names 
//        to orthogroup names.
// @param lookup_orthogroup_to_idx A named integer vector mapping orthogroup names 
//        to row indices in the final matrix.
// @param verbose boolean, debugging and print information if true
//
// @return A new \code{dgTMatrix} object with orthogroups as rows and cells as columns.
//
// @examples
// \dontrun{
//   convert_gene_to_oma_expression(expr, genes, oma_names, gene2oma, oma2idx)
// }
//
// [[Rcpp::export]]
Rcpp::S4 convert_gene_to_og_expression(Rcpp::S4 gene_expression_matrix,
                                       Rcpp::CharacterVector genes_names,
                                       Rcpp::CharacterVector orthogroup_names,
                                       Rcpp::CharacterVector lookup_gene_to_orthogroup,
                                       Rcpp::IntegerVector lookup_orthogroup_to_idx,
                                       bool verbose)
{
    // Create hashmap and lookup tables
    std::unordered_map<std::string, std::string> gene_to_orthogroup 
                            = gene_to_orthogroup_map(lookup_gene_to_orthogroup);
    std::unordered_map<std::string, int> orthogroup_to_idx
                            = name_to_idx_map(lookup_orthogroup_to_idx);

    // Extract tool variables from gene expression matrix
    const Rcpp::NumericVector x = gene_expression_matrix.slot("x");       // non-zero values
    const Rcpp::IntegerVector j = gene_expression_matrix.slot("j");       // column indices
    const Rcpp::IntegerVector p = gene_expression_matrix.slot("p");       // row pointers
    const Rcpp::IntegerVector dim = gene_expression_matrix.slot("Dim");   // dimensions
    const Rcpp::List dimnames = gene_expression_matrix.slot("Dimnames");  

    const int n_genes = dim[0];
    const int n_cells = dim[1];

    if(verbose){
        Rcpp::Rcout << "Number of cells : " << n_cells << "\n";
        Rcpp::Rcout << "Number of genes : " << n_genes << "\n";
    }
    const int n_orthogroups = orthogroup_names.size();

    // Initialize orthogroup expression matrix variables
    // Strategy: vector (n = cells) of hashmap (OMA -> expression)
    // Important pre-allocate the memory to reduce the runtime
    std::vector<std::unordered_map<int, double>> cellular_OrthoMAP_vector(n_cells);

    // Iterate over the genes and thus orthogroup groups! (dgRMatrix format)
    Progress pb(n_genes, verbose);
    for(int gene = 0; gene < n_genes; gene++){
        // Check if gene is in Orthogroup
        if(Progress::check_abort()) return NULL;
        std::string gene_name = Rcpp::as<std::string>(genes_names[gene]);
        auto gene_it = gene_to_orthogroup.find(gene_name);
        if(gene_it != gene_to_orthogroup.end()){
            // get index of the gene expressions in the cells
            int start_it = p[gene];
            int end_it = p[gene + 1];
            // get orthogroup idx
            std::string orthogroup_name = gene_to_orthogroup[gene_name];
            int orthogroup_idx = orthogroup_to_idx[orthogroup_name];
            // iterate over the expressions of the gene
            // can be parallelized using omp
            #pragma omp parallel for
            for(int it = start_it; it < end_it; it++){
                // get cell idx and expression values
                int cell = j[it];
                double expression = x[it];
                // find the iterator from cellular hashmap to orthogorup
                if(cellular_OrthoMAP_vector[cell].count(orthogroup_idx)){
                    // update when hashmap contains value
                    if((cellular_OrthoMAP_vector[cell])[orthogroup_idx] < expression){
                        (cellular_OrthoMAP_vector[cell])[orthogroup_idx] = expression;
                    }
                } 
                else {
                    // insert in hashmap
                    cellular_OrthoMAP_vector[cell].insert({orthogroup_idx, expression});
                }
            }
        }
        if (verbose) pb.increment();
    }
    if(verbose) Rcpp::Rcout << "Generate Orthogroup Expression Matrix \n\n";

    // Create orthogroup expression matrix in form of COO
    // pre-allocate memory of triplet vectors for each cells (n_cells is known)
    // triplet<int, int, double> in format (i, j, x)
    Rcpp::S4 orthogroup_expression_matrix("dgTMatrix");
    std::vector<int> i_out;     // cols (orthogroup)
    std::vector<int> j_out;     // rows (cell)
    std::vector<double> x_out;  // expression

    // calculate nnz and pre-allocate for runtime enhancement
    int nnz = 0;
    #pragma omp parallel reduction(+:nnz)
    for(int cell_id = 0; cell_id < n_cells; cell_id++){
        auto cell_hashmap = cellular_OrthoMAP_vector[cell_id];
        nnz += cell_hashmap.size();
    }
    i_out.reserve(nnz);
    j_out.reserve(nnz);
    x_out.reserve(nnz);

    // iterate over all oma expressions
    for(int cell_id = 0; cell_id < n_cells; cell_id++){
        // pre-allocate memory for orthogroups foreach cell in triplet vector
        auto cell_hashmap = cellular_OrthoMAP_vector[cell_id];
        // traverse over all orthogroups in cellular hashmap
        for(auto it = cell_hashmap.begin(); it != cell_hashmap.end(); it++){
            i_out.push_back(it->first);     // orthogroup
            j_out.push_back(cell_id);       // cell
            x_out.push_back(it->second);    // expression
        }
    }

    // build sparse dgTMatrix
    orthogroup_expression_matrix.slot("i") = Rcpp::IntegerVector(i_out.begin(), i_out.end());
    orthogroup_expression_matrix.slot("j") = Rcpp::IntegerVector(j_out.begin(), j_out.end());
    orthogroup_expression_matrix.slot("x") = Rcpp::NumericVector(x_out.begin(), x_out.end());
    orthogroup_expression_matrix.slot("Dim") = Rcpp::IntegerVector::create(n_orthogroups, n_cells);
    orthogroup_expression_matrix.slot("Dimnames") = Rcpp::List::create(orthogroup_names, dimnames[1]);
    return orthogroup_expression_matrix;
}
