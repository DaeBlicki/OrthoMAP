# Resource management using SLURM

## Motivation
Computer evaluation uses resources such as **Memory** and **CPU**. SLURM manages resources and distibutes them on the
users. When asking for to less memory, the program may crash. When you ask for way too much, the job is greedy and
prevent other users job scheduling. The same can be said for number requested cores. In high-performance computing,
the number of allocated cores should ideally reduce the runtime be the same factor. In very optimized applications,
the runtime is even more reduced as the theoretical ideal runtime due to super-linear scalability (optimal cashe usage).
However, asking for too much, is again very greedy and other users that using the SLURM resourceful must wait.

Observe the weak and strong scalability would be asked to much. However, I can provide "mindful" allocation by constantly
observing the resource efficiency of my jobs. For that, use following command after the job is done:

```
seff <job-id>
```
