# UNSTEADY-INCOMPRESSIBLE-NAVIER-STOKES in 2D and 3D
This project aims to solve the unsteady, incompressible Navier-Stokes equations using the finite element method. The focus is on simulating the benchmark problem "**Flow past a Cylinder**" in two and three dimensions.

## Strong formulation
$$ \begin{cases} 
\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla) \mathbf{u} - \nu \Delta \mathbf{u} + \nabla p = \mathbf{f} & \text{in } \Omega \, \\ 
\nabla \cdot \mathbf{u} = 0 & \text{in } \Omega \,\\ 
\mathbf{u} = \mathbf{g} & \text{on } \Gamma_{\mathrm{D}} \subset \partial \Omega \, \\ 
\nu \nabla \mathbf{u} \, \mathbf{n} - p \mathbf{n} = \mathbf{h} & \text{on } \Gamma_{\mathrm{N}} = \partial \Omega \backslash \Gamma_{\mathrm{D}} \,\\ 
\mathbf{u}(t=0) = \mathbf{u}_0 & \text{in } \Omega \, . 
\end{cases} $$

### COMPILE AND RUN
To compile and run the project these are the steps that need to be followed:

+ create and move inside build:<br> `mkdir build` `cd build`
+ load the dealii modules:<br> `module load gcc-glibc dealii`
+ build: <br>`cmake ..` `make`
+ run:
  - 2D Flow past a cylinder  -> `./navier_stokes2D`
  - 3D Flow past a cylinder  -> `./navier_stokes3D`
  - 3D Ethier-Steinmann cube -> `./convergence`

Output are saved in the _/build/output_ directory

### IMPORTANT COMANDS TO WORK ON THE CLUSTER
After building the code we can submit a job using
`sbatch run3D.sh` --- _specific for 3D test_
+ List of running jobs: <br>
`squeue -u <username>`
+ Stop a specific job (you can retrieve the id by listing the jobs first):<br>
`scancel <job_id>`
+ Follow in real-time the output:<br>
`tail -f myJob.out`
+ Visualize all the output file (non real-time):<br>
`cat myJob.out`

### MODIFY NUM OF THREADS
If you want to modify the number of nodes and tasks per node open the script <br>
`nano run3D.sh` <br>
and modify *ONLY* these two lines: <br>
`#SBATCH --nodes=2           # number of nodes` <br>
`#SBATCH --ntasks-per-node=32 # number of tasks per node` 
