#ifndef NAVIER_STOKES_HPP
#define NAVIER_STOKES_HPP

#include "Preconditioners.hpp"
#include "IncludesFile.hpp"


using namespace dealii;

// Class implementing a solver for the NS problem.
class NavierStokes
{
public:
  // Physical dimension (3D)
  static constexpr unsigned int dim = 3;

  // Function for the forcing term.
  class ForcingTerm : public Function<dim>
  {
  public:
    ForcingTerm()
    {
    }

    virtual void
    vector_value(const Point<dim> & /*p*/,
                 Vector<double> &values) const override
    {
      for (unsigned int i = 0; i < dim - 1; ++i)
        values[i] = 0.0;

      values[dim - 1] = -g;
    }

    virtual double
    value(const Point<dim> & /*p*/,
          const unsigned int component = 0) const override
    {
      if (component == dim - 1)
        return -g;
      else
        return 0.0;
    }

  protected:
    const double g = 0.0;
  };


    //Exact sol
    class ExactSolution : public Function<dim>
    {
    public:
      // Constructor.
      ExactSolution() : Function<dim>(dim + 1)
      {
      }
  
      virtual void
      vector_value(const Point<dim> & p, Vector<double> &values) const override
      {
      double factor = -(a * a * std::exp(-2 * nu * b * b * get_time())) / 2.0;
      double term1 = 2.0 * std::sin(a * p[0] + b * p[1]) * std::cos(a * p[2] + b * p[0]) * std::exp(a * (p[1] + p[2]));
      double term2 = 2.0 * std::sin(a * p[1] + b * p[2]) * std::cos(a * p[0] + b * p[1]) * std::exp(a * (p[0] + p[2]));
      double term3 = 2.0 * std::sin(a * p[2] + b * p[0]) * std::cos(a * p[1] + b * p[2]) * std::exp(a * (p[0] + p[1]));
      double term4 = std::exp(2 * a * p[0]) + std::exp(2 * a * p[1]) + std::exp(2 * a * p[2]);
      double pressure = (factor * (term1 + term2 + term3 + term4));

        values[0] = -a * std::exp( -nu * b * b * get_time() ) * ( std::exp(a * p[0]) * std::sin(a * p[1] + b * p[2]) + std::exp(a * p[2]) * std::cos(a * p[0] + b * p[1]) );
        values[1] = -a * std::exp( -nu * b * b * get_time() ) * ( std::exp(a * p[1]) * std::sin(a * p[2] + b * p[0]) + std::exp(a * p[0]) * std::cos(a * p[1] + b * p[2]) );
        values[2] = -a * std::exp( -nu * b * b * get_time() ) * ( std::exp(a * p[2]) * std::sin(a * p[0] + b * p[1]) + std::exp(a * p[1]) * std::cos(a * p[2] + b * p[0]) );
        values[3] = pressure;
      }
      virtual double
      value(const Point<dim> & p, const unsigned int component) const override
      {
      double factor = -(a * a * std::exp(-2 * nu * b * b * get_time())) / 2.0;
      double term1 = 2 * std::sin(a * p[0] + b * p[1]) * std::cos(a * p[2] + b * p[0]) * std::exp(a * (p[1] + p[2]));
      double term2 = 2 * std::sin(a * p[1] + b * p[2]) * std::cos(a * p[0] + b * p[1]) * std::exp(a * (p[0] + p[2]));
      double term3 = 2 * std::sin(a * p[2] + b * p[0]) * std::cos(a * p[1] + b * p[2]) * std::exp(a * (p[0] + p[1]));
      double term4 = std::exp(2 * a * p[0]) + std::exp(2 * a * p[1]) + std::exp(2 * a * p[2]);
      double pressure = factor * (term1 + term2 + term3 + term4);

          if (component == 0)
          {
              return -a * std::exp( -nu * b * b * get_time() ) * ( std::exp(a * p[0]) * std::sin(a * p[1] + b * p[2]) + std::exp(a * p[2]) * std::cos(a * p[0] + b * p[1]) );
          }
          else if (component == 1)
          {
              return -a * std::exp( -nu * b * b * get_time() ) * ( std::exp(a * p[1]) * std::sin(a * p[2] + b * p[0]) + std::exp(a * p[0]) * std::cos(a * p[1] + b * p[2]) );
          }
          else if (component == 2)
          {
              return -a * std::exp( -nu * b * b * get_time() ) * ( std::exp(a * p[2]) * std::sin(a * p[0] + b * p[1]) + std::exp(a * p[1]) * std::cos(a * p[2] + b * p[0]) );
          }
          else
          {
              //return 0;
              return pressure;
          }
        }
      
      private:
      const double nu = 1e-2;
      const double a = M_PI / 4.0;
      const double b = M_PI / 2.0; 

          // Gradient evaluation.
    virtual Tensor<2, dim>
    gradient_tensor(const Point<dim> &p,
             const unsigned int /*component*/ = 0) const
      {
      Tensor<2, dim> values;

      values[0][0] = -a * std::exp( -nu * b * b * get_time() ) * ( a * std::exp(a * p[0]) * std::sin(a * p[1] + b * p[2]) - a * std::exp(a * p[2]) * std::sin(a * p[0] + b * p[1]) );
      values[0][1] = -a * std::exp( -nu * b * b * get_time() ) * ( a * std::exp(a * p[0]) * std::cos(a * p[1] + b * p[2]) - b * std::exp(a * p[2]) * std::sin(a * p[0] + b * p[1]) );
      values[0][2] = -a * std::exp( -nu * b * b * get_time() ) * ( b * std::exp(a * p[0]) * std::cos(a * p[1] + b * p[2]) + a * std::exp(a * p[2]) * std::cos(a * p[0] + b * p[1]) );
        
      values[1][0] = -a * std::exp( -nu * b * b * get_time() ) * ( b * std::exp(a * p[1]) * std::cos(a * p[2] + b * p[0]) + a * std::exp(a * p[0]) * std::cos(a * p[1] + b * p[2]) );
      values[1][1] = -a * std::exp( -nu * b * b * get_time() ) * ( a * std::exp(a * p[1]) * std::sin(a * p[2] + b * p[0]) - a * std::exp(a * p[0]) * std::sin(a * p[1] + b * p[2]) );
      values[1][2] = -a * std::exp( -nu * b * b * get_time() ) * ( a * std::exp(a * p[1]) * std::cos(a * p[2] + b * p[0]) - b * std::exp(a * p[0]) * std::sin(a * p[1] + b * p[2]) );
       
      values[2][0] = -a * std::exp( -nu * b * b * get_time() ) * ( a * std::exp(a * p[2]) * std::cos(a * p[0] + b * p[1]) - b * std::exp(a * p[1]) * std::sin(a * p[2] + b * p[0]) );
      values[2][1] = -a * std::exp( -nu * b * b * get_time() ) * ( b * std::exp(a * p[2]) * std::cos(a * p[0] + b * p[1]) + a * std::exp(a * p[1]) * std::cos(a * p[2] + b * p[0]) );
      values[2][2] = -a * std::exp( -nu * b * b * get_time() ) * ( a * std::exp(a * p[2]) * std::sin(a * p[0] + b * p[1]) - a * std::exp(a * p[1]) * std::sin(a * p[2] + b * p[0]) );
        
      values[3][0] = 0.0;
      values[3][1] = 0.0;
      values[3][2] = 0.0;    

      return values;
      }

    virtual Tensor<1, dim>
     gradient(const Point<dim> &p, const unsigned int component = 0) const override
      {
      
      Tensor<2, dim> grad_tensor = gradient_tensor(p);
      Tensor<1, dim> grad_component;

      for (unsigned int i = 0; i < dim; ++i)
        grad_component[i] = grad_tensor[component][i];

    return grad_component;
    }


    };

// Neumann boundary conditions.
class FunctionH : public Function<dim>
{
public:
  // Constructor.
  FunctionH() : Function<dim>(dim)
  {
  }

  virtual void
  vector_value(const Point<dim> & p, Vector<double> &values) const override
  {
    //expression of H 
    //H=ν∂u/∂n​−pn =  ν∂u/∂n - (0,p(x,-1,z),0)
      double factor = -(a * a * std::exp(-2 * nu * b * b * get_time())) / 2.0;
      double term1 = 2.0 * std::sin(a * p[0] + b * p[1]) * std::cos(a * p[2] + b * p[0]) * std::exp(a * (p[1] + p[2]));   
      double term2 = 2.0 * std::sin(a * p[1] + b * p[2]) * std::cos(a * p[0] + b * p[1]) * std::exp(a * (p[0] + p[2]));
      double term3 = 2.0 * std::sin(a * p[2] + b * p[0]) * std::cos(a * p[1] + b * p[2]) * std::exp(a * (p[0] + p[1]));
      double term4 = std::exp(2.0 * a * p[0]) + std::exp(2.0 * a * p[1]) + std::exp(2.0 * a * p[2]);
      double pressure = factor * (term1 + term2 + term3 + term4);

      values[0] =  - nu * a * std::exp(-nu * b * b * get_time()) * ( a * std::exp(a * p[0]) * std::cos(a * p[1] + b * p[2]) - b * std::exp(a * p[2]) * std::sin(a * p[0] + b * p[1]) );
      values[1] = - nu * a * std::exp(-nu * b * b * get_time()) * (a * std::exp(a * p[1]) * std::sin(a * p[2] + b * p[0]) - a * std::exp(a * p[0]) * std::sin(a * p[1] + b * p[2])) - pressure;
      values[2] = - nu * a * std::exp(-nu * b * b * get_time()) * (b * std::exp(a * p[2]) * std::cos(a * p[0] + b * p[1]) + a * std::exp(a * p[1]) * std::cos(a * p[2] + b * p[0]));

      }
      virtual double value(const Point<dim> &p, const unsigned int component) const override
      {
      double factor = -(a * a * std::exp(-2 * nu * b * b * get_time())) / 2.0;
      double term1 = 2.0 * std::sin(a * p[0] + b * p[1]) * std::cos(a * p[2] + b * p[0]) * std::exp(a * (p[1] + p[2]));   
      double term2 = 2.0 * std::sin(a * p[1] + b * p[2]) * std::cos(a * p[0] + b * p[1]) * std::exp(a * (p[0] + p[2]));
      double term3 = 2.0 * std::sin(a * p[2] + b * p[0]) * std::cos(a * p[1] + b * p[2]) * std::exp(a * (p[0] + p[1]));
      double term4 = std::exp(2.0 * a * p[0]) + std::exp(2.0 * a * p[1]) + std::exp(2.0 * a * p[2]);
      double pressure = factor * (term1 + term2 + term3 + term4);
      if (component == 0)
          {
              return- nu * a * std::exp(-nu * b * b * get_time()) * ( a * std::exp(a * p[0]) * std::cos(a * p[1] + b * p[2]) - b * std::exp(a * p[2]) * std::sin(a * p[0] + b * p[1]) );
          }
          else if (component == 1)
          {
              return - nu * a * std::exp(-nu * b * b * get_time()) * (a * std::exp(a * p[1]) * std::sin(a * p[2] + b * p[0]) - a * std::exp(a * p[0]) * std::sin(a * p[1] + b * p[2])) - pressure;
          }
          else 
              return - nu * a * std::exp(-nu * b * b * get_time()) * (b * std::exp(a * p[2]) * std::cos(a * p[0] + b * p[1]) + a * std::exp(a * p[1]) * std::cos(a * p[2] + b * p[0]));
            }
                
            private:
            const double nu = 1e-2;
            const double a = M_PI / 4.0;
            const double b = M_PI / 2.0;  

    };

  

  // Function for the initial condition.
  class FunctionU0 : public Function<dim>
  {
  public:
    FunctionU0()
    {
    }

    virtual double
    value(const Point<dim> & p, const unsigned int component) const override
    {
    double factor = -(a * a * std::exp(-2 * nu * b * b * 0.0)) / 2.0;
    double term1 = 2 * std::sin(a * p[0] + b * p[1]) * std::cos(a * p[2] + b * p[0]) * std::exp(a * (p[1] + p[2]));
    double term2 = 2 * std::sin(a * p[1] + b * p[2]) * std::cos(a * p[0] + b * p[1]) * std::exp(a * (p[0] + p[2]));
    double term3 = 2 * std::sin(a * p[2] + b * p[0]) * std::cos(a * p[1] + b * p[2]) * std::exp(a * (p[0] + p[1]));
    double term4 = std::exp(2 * a * p[0]) + std::exp(2 * a * p[1]) + std::exp(2 * a * p[2]);
    double pressure = factor * (term1 + term2 + term3 + term4);

        if (component == 0)
        {
            return -a * std::exp(-nu * b * b * 0.0) * 
                  (std::exp(a * p[0]) * std::sin(a * p[1] + b * p[2]) + 
                    std::exp(a * p[2]) * std::cos(a * p[0] + b * p[1]));
        }
        else if (component == 1)
        {
            return -a * std::exp(-nu * b * b * 0.0) * 
                  (std::exp(a * p[1]) * std::sin(a * p[2] + b * p[0]) + 
                    std::exp(a * p[0]) * std::cos(a * p[1] + b * p[2]));
        }
        else if (component == 2)
        {
            return -a * std::exp(-nu * b * b * 0.0) * 
                  (std::exp(a * p[2]) * std::sin(a * p[0] + b * p[1]) + 
                    std::exp(a * p[1]) * std::cos(a * p[2] + b * p[0]));
        }
        else
        {
            return pressure;
        }
      }
    virtual void
    vector_value(const Point<dim> & p, Vector<double> &values) const override
    {
    double factor = -(a * a * std::exp(-2 * nu * b * b * 0.0)) / 2.0;
    double term1 = 2.0 * std::sin(a * p[0] + b * p[1]) * std::cos(a * p[2] + b * p[0]) * std::exp(a * (p[1] + p[2]));
    double term2 = 2.0 * std::sin(a * p[1] + b * p[2]) * std::cos(a * p[0] + b * p[1]) * std::exp(a * (p[0] + p[2]));
    double term3 = 2.0 * std::sin(a * p[2] + b * p[0]) * std::cos(a * p[1] + b * p[2]) * std::exp(a * (p[0] + p[1]));
    double term4 = std::exp(2 * a * p[0]) + std::exp(2 * a * p[1]) + std::exp(2 * a * p[2]);
    double pressure = (factor * (term1 + term2 + term3 + term4));

      values[0] = -a * std::exp( -nu * b * b * 0.0 ) * ( std::exp(a * p[0]) * std::sin(a * p[1] + b * p[2]) + std::exp(a * p[2]) * std::cos(a * p[0] + b * p[1]) );
      values[1] = -a * std::exp( -nu * b * b * 0.0 ) * ( std::exp(a * p[1]) * std::sin(a * p[2] + b * p[0]) + std::exp(a * p[0]) * std::cos(a * p[1] + b * p[2]) );
      values[2] = -a * std::exp( -nu * b * b * 0.0 ) * ( std::exp(a * p[2]) * std::sin(a * p[0] + b * p[1]) + std::exp(a * p[1]) * std::cos(a * p[2] + b * p[0]) );
      values[3] = pressure;
    }
    private:
    const double nu = 1e-2;
    const double a = M_PI / 4.0;
    const double b = M_PI / 2.0; 
  };

  
  // Constructor.
  NavierStokes(const std::string &mesh_file_name_,
               const unsigned int &degree_velocity_,
               const unsigned int &degree_pressure_,
               const double &T_,
               const double &deltat_)
      : mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)), mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)), pcout(std::cout, mpi_rank == 0), T(T_), mesh_file_name(mesh_file_name_), degree_velocity(degree_velocity_), degree_pressure(degree_pressure_), deltat(deltat_), mesh(MPI_COMM_WORLD)
  {
  }

  // Setup system.
  void
  setup();

  // Solve system.
  void
  solve();

  // Compute the error.
  double
  compute_error(const VectorTools::NormType &norm_type);


  std::vector<double> time_prec;
  std::vector<double> time_solve;

protected:
  // Assemble system. We also assemble the pressure mass matrix (needed for the
  // preconditioner).
  void
  assemble(const double &time);

  void
  assemble_time_step(const double &time);

  // Solve the problem for one time step.
  void
  solve_time_step();

  // Output results.
  void
  output(const unsigned int &time_step) const;


  // MPI parallel. /////////////////////////////////////////////////////////////

  // Number of MPI processes.
  const unsigned int mpi_size;

  // This MPI process.
  const unsigned int mpi_rank;

  // Parallel output stream.
  ConditionalOStream pcout;

  // Problem definition. ///////////////////////////////////////////////////////

  // Kinematic viscosity [m2/s].
  const double nu = 1e-2;

  // Density
  const double rho = 1.;

  // Forcing term.
  ForcingTerm forcing_term;



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

  ExactSolution exact_solution;



  // h(x).
  FunctionH function_h;

  // Initial condition.
  FunctionU0 u_0;

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

  std::vector<Tensor<1, dim>> previous_velocity_values;

   std::vector<Tensor<2, dim>> previous_gradient_velocity_values;

  // System matrix.
  TrilinosWrappers::BlockSparseMatrix system_matrix;
  // A
  TrilinosWrappers::BlockSparseMatrix stiffness_matrix;
  // M
  TrilinosWrappers::BlockSparseMatrix mass_matrix;
  // C(u_k)
  TrilinosWrappers::BlockSparseMatrix convection_matrix;

  // Pressure mass matrix, needed for preconditioning. We use a block matrix for
  // convenience, but in practice we only look at the pressure-pressure block.
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
