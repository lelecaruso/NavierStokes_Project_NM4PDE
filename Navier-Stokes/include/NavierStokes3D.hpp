#ifndef NAVIER_STOKES_HPP
#define NAVIER_STOKES_HPP

#include "Preconditioners.hpp"
#include "IncludesFile.hpp"

using namespace dealii;

// Class implementing a solver for the Stokes problem.
class NavierStokes
{
public:
  // Physical dimension (3D)
  static constexpr unsigned int dim = 3;

  // Function for inlet velocity. This actually returns an object with four dimensions
  class InletVelocity : public Function<dim>
  {
  public:
  InletVelocity(unsigned int test_case_ = 2)
    : Function<dim>(dim + 1) , test_case(test_case_)
  {
  }
    
    virtual void
    vector_value(const Point<dim> & p, Vector<double> &values) const override
    {
      switch(this->test_case) {
        case 1:
          values[0] = 0.0;  
          break;
        case 3:
          values[0] = 16.0 * u_m * p[1] * p[2] * ( H - p[2] ) * (H - p[1]) * std::sin(M_PI * get_time() / 8.0) / (H*H*H*H);
          break;
        case 2:
        default:
          values[0] = 16.0 * u_m * p[1] * p[2] * ( H - p[2] ) * (H - p[1]) / (H*H*H*H);
          break;
      }
      
      for (unsigned int i = 1; i < dim + 1; ++i)
        values[i] = 0.0;
    }
    
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override
    {
      if (component == 0) {
        switch(this->test_case) {
          case 1:
            return 0.0; 
          case 3:
            return 16.0 * u_m * p[1] * p[2] * ( H - p[2] ) * (H - p[1]) * std::sin(M_PI * get_time() / 8.0) / (H*H*H*H);
          case 2:
          default:
            double val = 16.0 * u_m * p[1] * p[2] * ( H - p[2] ) * (H - p[1]) / (H*H*H*H);
            return val ; 
        }
      }
      else
        return 0;
    }
    
    double getMeanVelocity() const
    {
      switch(this->test_case) {
        case 1:
          return 0.0;  
        case 3:
          return 4.0 * u_m * std::sin(get_time()*M_PI/8.0) / 9.0;
        case 2:
        default:
          return (4.0 * u_m )  / (9.0);
      }
    }
    
  protected:
    int test_case;
    double H = 0.41;
    double u_m = 9.0;
  };


  // Constructor.
  NavierStokes( const std::string &mesh_file_name_,
    const unsigned int &degree_velocity_,
    const unsigned int &degree_pressure_,
    const double &T_,
    const double &deltat_,
    const int test_case_ = 2)  // Default to test case 2 if not specified
      : 
      test_case(test_case_),
      mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)), 
      mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)), 
      pcout(std::cout, mpi_rank == 0), 
      inlet_velocity(test_case),
      T(T_), 
      mesh_file_name(mesh_file_name_), 
      degree_velocity(degree_velocity_), 
      degree_pressure(degree_pressure_), 
      deltat(deltat_), 
      mesh(MPI_COMM_WORLD)
      
{}


  // Setup system.
  void
  setup();

  // Solve system.
  void
  solve();

  
  std::vector<double> vec_drag;
  std::vector<double> vec_lift;
  std::vector<double> vec_drag_coeff;
  std::vector<double> vec_lift_coeff;

  std::vector<double> time_prec;
  std::vector<double> time_solve;

protected:
  // Assemble system the first time to create mass-stiffness matrixes 
  void
  assemble(const double &time);
  // Assemble at each time step to only compute the convection matrix that changes overtime
  void
  assemble_time_step(const double &time);

  // Solve the problem for one time step.
  void
  solve_time_step();

  // Output results.
  void
  output(const unsigned int &time_step) const;

  // Compute Lift and Drag coefficients and forces
  std::vector<double>
  compute_forces();

  void
  compute_pressure_difference();

  // MPI parallel. /////////////////////////////////////////////////////////////

  unsigned int test_case;
  // Number of MPI processes.
  const unsigned int mpi_size;

  // This MPI process.
  const unsigned int mpi_rank;

  // Parallel output stream.
  ConditionalOStream pcout;

  // Problem definition. ///////////////////////////////////////////////////////

  // Kinematic viscosity [m2/s].
  const double nu = 1e-3;

  // Density
  const double rho = 1.;

  // Forcing term.
  Functions::ZeroFunction<dim> forcing_term;

  // Inlet velocity.
  InletVelocity inlet_velocity;

  // Final time.
  const double T;

  double drag;
  double lift;

  // Discretization. ///////////////////////////////////////////////////////////

  // Mesh file name.
  const std::string mesh_file_name;

  // Polynomial degree used for velocity.
  const unsigned int degree_velocity;

  // Polynomial degree used for pressure.
  const unsigned int degree_pressure;

  // TIme step.
  const double deltat;

  // g(x).
  Functions::ZeroFunction<dim> function_g;

  // h(x).
  Functions::ZeroFunction<dim> function_h;

  // Initial condition.
  Functions::ZeroFunction<dim> u_0;

  // Mesh.
  parallel::fullydistributed::Triangulation<dim> mesh;

  // Finite element space.
  std::unique_ptr<FiniteElement<dim>> fe;

  // Quadrature formula.
  std::unique_ptr<Quadrature<dim>> quadrature;

  // Quadrature formula used on boundary lines.
  std::unique_ptr<Quadrature<dim - 1>> quadrature_boundary;

  // DoF handler.
  DoFHandler<dim> dof_handler;

  // DoFs owned by current process.
  IndexSet locally_owned_dofs;

  // DoFs owned by current process in the velocity and pressure blocks.
  std::vector<IndexSet> block_owned_dofs;

  // DoFs relevant to the current process (including ghost DoFs).
  IndexSet locally_relevant_dofs;

  // DoFs relevant to current process in the velocity and pressure blocks.
  std::vector<IndexSet> block_relevant_dofs;

  // System matrix.
  TrilinosWrappers::BlockSparseMatrix system_matrix;
  // Stiffness Matrix
  TrilinosWrappers::BlockSparseMatrix stiffness_matrix;
  // Mass Matrix
  TrilinosWrappers::BlockSparseMatrix mass_matrix;
  // Convection matrix at step (u_k)
  TrilinosWrappers::BlockSparseMatrix convection_matrix;

  // Pressure mass matrix, needed for preconditioning, also referred as Mp
  TrilinosWrappers::BlockSparseMatrix pressure_mass;

  // Right-hand side vector in the linear system.
  TrilinosWrappers::MPI::BlockVector system_rhs;

  // System solution (without ghost elements).
  TrilinosWrappers::MPI::BlockVector solution_owned;

  // System solution (including ghost elements).
  TrilinosWrappers::MPI::BlockVector solution;

  TrilinosWrappers::MPI::BlockVector previous_solution;

};

#endif
