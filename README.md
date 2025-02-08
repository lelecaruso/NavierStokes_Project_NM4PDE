# UNSTEADY-INCOMPRESSIBLE-NAVIER-STOKES in 2D and 3D
This project aims to solve the unsteady, incompressible Navier-Stokes equations using the finite element method. The focus is on simulating the benchmark problem "flow past a cylinder" in two or three dimensions.

## Strong formulation
$$ \begin{cases}
  -\mu \nabla \vec{u} + (\rho \vec{u} \cdot \nabla) \vec{u} + \nabla p = \vec{f} & \text{in }\Omega \\
  \nabla \cdot \vec{u} = 0 & \text{in }\Omega \\
  \mu \nabla \vec{u} \cdot \widehat{\vec{n}} - p \cdot \widehat{\vec{n}} = -p_\text{out} \cdot \widehat{\vec{n}} & \text{on } \Gamma_0 \\
  \vec{u} = \vec{u}_\text{in} & \text{on } \Gamma_1 \\
  \vec{u} \cdot \widehat{\vec{n}} = \vec{0} & \text{on } \Gamma_2 \cup \Gamma_3 \cup \Gamma_4 \cup \Gamma_5 \\
  (\mu \nabla \vec{u} \cdot \widehat{\vec{n}} - p  \cdot \widehat{\vec{n}}) \cdot \widehat{\vec{t}}_i = \vec{0} & \text{on } \Gamma_i,\ i=2,3,4,5
\end{cases} $$

### COMPILE AND RUN
To compile and run the project these are the steps that need to be followed:

+ create build: `mkdir build`
+ move inside build
+ load the dealii modules:`module load gcc-glibc dealii`
+ execute: `cmake ..` `make`
+ run 2D problem: `./navier_stokes2D`
