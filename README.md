# UNSTEADY-INCOMPRESSIBLE-NAVIER-STOKES in 2D and 3D
This project aims to solve the unsteady, incompressible Navier-Stokes equations using the finite element method. The focus is on simulating the benchmark problem "flow past a cylinder" in two or three dimensions.

## Strong Formulation

$$
egin{cases}
  -\mu 
abla ec{u} + (
ho ec{u} \cdot 
abla) ec{u} + 
abla p = ec{f}, & 	ext{in } \Omega, \
  
abla \cdot ec{u} = 0, & 	ext{in } \Omega, \
  \mu 
abla ec{u} \cdot \widehat{ec{n}} - p \widehat{ec{n}} = -p_	ext{out} \widehat{ec{n}}, & 	ext{on } \Gamma_0, \
  ec{u} = ec{u}_	ext{in}, & 	ext{on } \Gamma_1, \
  ec{u} \cdot \widehat{ec{n}} = 0, & 	ext{on } \Gamma_2 \cup \Gamma_3 \cup \Gamma_4 \cup \Gamma_5, \
  (\mu 
abla ec{u} \cdot \widehat{ec{n}} - p \widehat{ec{n}}) \cdot \widehat{ec{t}}_i = 0, & 	ext{on } \Gamma_i,\ i = 2, 3, 4, 5.
\end{cases}
$$

### COMPILE AND RUN
To compile and run the project these are the steps that need to be followed:

+ create build: `mkdir build`
+ move inside build
+ load the dealii modules:`module load gcc-glibc dealii`
+ execute: `cmake ..` `make`
+ run 2D problem: `./navier_stokes2D`
