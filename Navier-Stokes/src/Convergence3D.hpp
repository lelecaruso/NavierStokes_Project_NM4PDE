#ifndef NAVIER_STOKES_HPP
#define NAVIER_STOKES_HPP

#include <deal.II/base/timer.h>

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/grid/tria.h>  
#include <deal.II/fe/fe_tools.h>  

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_in.h>

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_tools.h>

#include <fstream>
#include <iostream>
#include <vector>

using namespace dealii;

// Class implementing a solver for the Stokes problem.
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

  // Homog. Dirichlet boundary conditions.
  class FunctionG : public Function<dim>
  {
  public:
    // Constructor.
    FunctionG() : Function<dim>(dim + 1)
    {
    }

    virtual void
    vector_value(const Point<dim> & /*p*/, Vector<double> &values) const override
    {
      values[0] = 0.;
      values[1] = 0.;
      values[2] = 0.;
      values[3] = 0.;
    }

    virtual double
    value(const Point<dim> & /*p*/, const unsigned int /*component*/) const override
    {
      return 0.;
    }
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
        // values[3] = 0;
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
    //vettore normale al piano y = - 1 è n = (0,1,0)
    //p[1] = -1; //y=-1

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

  
  // Since we're working with block matrices, we need to make our own
  // preconditioner class. A preconditioner class can be any class that exposes
  // a vmult method that applies the inverse of the preconditioner.

  // Identity preconditioner.
  class PreconditionIdentity
  {
  public:
    // Application of the preconditioner: we just copy the input vector (src)
    // into the output vector (dst).
    void
    vmult(TrilinosWrappers::MPI::Vector &dst,
          const TrilinosWrappers::MPI::Vector &src) const
    {
      dst = src;
    }

  protected:
  };
  class PreconditionBlockIdentity
  {
  public:
    // Application of the preconditioner: we just copy the input vector (src)
    // into the output vector (dst).
    void
    vmult(TrilinosWrappers::MPI::BlockVector &dst,
          const TrilinosWrappers::MPI::BlockVector &src) const
    {
      dst = src;
    }

  protected:
  };

  // Block-diagonal preconditioner.
  class PreconditionBlockDiagonal
  {
  public:
    // Initialize the preconditioner, given the velocity stiffness matrix, the
    // pressure mass matrix.
    void
    initialize(const TrilinosWrappers::SparseMatrix &velocity_stiffness_,
               const TrilinosWrappers::SparseMatrix &pressure_mass_)
    {
      velocity_stiffness = &velocity_stiffness_;
      pressure_mass = &pressure_mass_;

      preconditioner_velocity.initialize(velocity_stiffness_);
      preconditioner_pressure.initialize(pressure_mass_);
    }

    // Application of the preconditioner.
    void
    vmult(TrilinosWrappers::MPI::BlockVector &dst,
          const TrilinosWrappers::MPI::BlockVector &src) const
    {
      SolverControl solver_control_velocity(10000,
                                            1e-2 /** src.block(0).l2_norm()*/);
      SolverCG<TrilinosWrappers::MPI::Vector> solver_cg_velocity(
          solver_control_velocity);
      solver_cg_velocity.solve(*velocity_stiffness,
                               dst.block(0),
                               src.block(0),
                               preconditioner_velocity);

      SolverControl solver_control_pressure(10000,
                                            1e-2 /** src.block(1).l2_norm()*/);
      SolverCG<TrilinosWrappers::MPI::Vector> solver_cg_pressure(
          solver_control_pressure);
      solver_cg_pressure.solve(*pressure_mass,
                               dst.block(1),
                               src.block(1),
                               preconditioner_pressure);
    }

  protected:
    // Velocity stiffness matrix.
    const TrilinosWrappers::SparseMatrix *velocity_stiffness;

    // Preconditioner used for the velocity block.
    TrilinosWrappers::PreconditionILU preconditioner_velocity;

    // Pressure mass matrix.
    const TrilinosWrappers::SparseMatrix *pressure_mass;

    // Preconditioner used for the pressure block.
    TrilinosWrappers::PreconditionILU preconditioner_pressure;
  };

  // Block-triangular preconditioner.
  class PreconditionBlockTriangular
  {
  public:
    // Initialize the preconditioner, given the velocity stiffness matrix, the
    // pressure mass matrix.
    void
    initialize(const TrilinosWrappers::SparseMatrix &velocity_stiffness_,
               const TrilinosWrappers::SparseMatrix &pressure_mass_,
               const TrilinosWrappers::SparseMatrix &B_)
    {
      velocity_stiffness = &velocity_stiffness_;
      pressure_mass = &pressure_mass_;
      B = &B_;

      preconditioner_velocity.initialize(velocity_stiffness_);
      preconditioner_pressure.initialize(pressure_mass_);
    }

    // Application of the preconditioner.
    void
    vmult(TrilinosWrappers::MPI::BlockVector &dst,
          const TrilinosWrappers::MPI::BlockVector &src) const
    {
      SolverControl solver_control_velocity(10000,
                                            1e-2 /** src.block(0).l2_norm()*/);
      SolverCG<TrilinosWrappers::MPI::Vector> solver_cg_velocity(
          solver_control_velocity);
      solver_cg_velocity.solve(*velocity_stiffness,
                               dst.block(0),
                               src.block(0),
                               preconditioner_velocity);

      tmp.reinit(src.block(1));
      B->vmult(tmp, dst.block(0));
      tmp.sadd(-1.0, src.block(1));

      SolverControl solver_control_pressure(10000,
                                            1e-2 /* * src.block(1).l2_norm()*/);
      SolverCG<TrilinosWrappers::MPI::Vector> solver_cg_pressure(
          solver_control_pressure);
      solver_cg_pressure.solve(*pressure_mass,
                               dst.block(1),
                               tmp,
                               preconditioner_pressure);
    }

  protected:
    // Velocity stiffness matrix.
    const TrilinosWrappers::SparseMatrix *velocity_stiffness;

    // Preconditioner used for the velocity block.
    TrilinosWrappers::PreconditionILU preconditioner_velocity;

    // Pressure mass matrix.
    const TrilinosWrappers::SparseMatrix *pressure_mass;

    // Preconditioner used for the pressure block.
    TrilinosWrappers::PreconditionILU preconditioner_pressure;

    // B matrix.
    const TrilinosWrappers::SparseMatrix *B;

    // Temporary vector.
    mutable TrilinosWrappers::MPI::Vector tmp;
  };
  class PreconditionSIMPLE
  {
  public:
    void
    initialize(const TrilinosWrappers::SparseMatrix &F_,
               const TrilinosWrappers::SparseMatrix &B_,
               const TrilinosWrappers::SparseMatrix &B_t,
               const TrilinosWrappers::MPI::BlockVector &sol_owned)
    {
      F = &F_;
      B = &B_;
      B_T = &B_t;

      diag_D_inv.reinit(sol_owned.block(0));

      for (unsigned int i : diag_D_inv.locally_owned_elements())
      {
        double temp = F->diag_element(i);
        diag_D_inv[i] = 1.0 / temp;
      }

      // Create S_tilde
      B_.mmult(S_tilde, B_t, diag_D_inv);

      // Initialize the preconditioners
      preconditioner_F.initialize(*F);
      preconditioner_S.initialize(S_tilde);
    }
    void
    vmult(TrilinosWrappers::MPI::BlockVector &dst,
          const TrilinosWrappers::MPI::BlockVector &src) const
    {
      const unsigned int maxiter = 10000;
      const double tol = 1e-2;
      SolverControl solver_F(maxiter, tol * src.block(0).l2_norm());

      SolverGMRES<TrilinosWrappers::MPI::Vector> solver_gmres(solver_F);

      // Store in temporaries the results
      TrilinosWrappers::MPI::Vector y_u = src.block(0);
      TrilinosWrappers::MPI::Vector y_p = src.block(1);

      TrilinosWrappers::MPI::Vector temp_1 = src.block(1);

      solver_gmres.solve(*F, y_u, src.block(0), preconditioner_F);

      B->vmult(temp_1, y_u);
      temp_1 -= src.block(1);

      SolverControl solver_S(maxiter, tol * temp_1.l2_norm());
      SolverCG<TrilinosWrappers::MPI::Vector> solver_cg(solver_S);
      solver_cg.solve(S_tilde, y_p, temp_1, preconditioner_S);

      dst.block(1) = y_p;
      dst.block(1) *= 1. / alpha;
      // temp_1.reinit(dst.block(0));

      B_T->vmult(dst.block(0), dst.block(1));
      // Cannot be same vector
      // D_inv.vmult(dst.block(0), temp_1);
      dst.block(0).scale(diag_D_inv);
      dst.block(0) -= y_u;
      dst.block(0) *= -1.;
    }

  protected:
    const double alpha = 0.5;

    const TrilinosWrappers::SparseMatrix *F;
    const TrilinosWrappers::SparseMatrix *B_T;
    const TrilinosWrappers::SparseMatrix *B;
    TrilinosWrappers::SparseMatrix S_tilde;
    TrilinosWrappers::MPI::Vector diag_D_inv;
    TrilinosWrappers::PreconditionILU preconditioner_F;
    TrilinosWrappers::PreconditionILU preconditioner_S;
  };

  class PreconditionaSIMPLE
  {
  public:
    void
    initialize(const TrilinosWrappers::SparseMatrix &F_,
               const TrilinosWrappers::SparseMatrix &B_,
               const TrilinosWrappers::SparseMatrix &B_t,
               const TrilinosWrappers::MPI::BlockVector &sol_owned)
    {
      F = &F_;
      B = &B_;
      B_T = &B_t;

      diag_D_inv.reinit(sol_owned.block(0));
      diag_D.reinit(sol_owned.block(0));

      for (unsigned int i : diag_D.locally_owned_elements())
      {
        double temp = F->diag_element(i);
        diag_D[i] = -temp;
        diag_D_inv[i] = 1.0 / temp;
      }

      B->mmult(S, *B_T, diag_D_inv);

      preconditionerF.initialize(*F);
      preconditionerS.initialize(S);
    }
    void
    vmult(TrilinosWrappers::MPI::BlockVector &dst,
          const TrilinosWrappers::MPI::BlockVector &src) const
    {

      const unsigned int maxiter = 10000;
      const double tol = 1e-2;
      SolverControl solver_F(maxiter, tol * src.block(0).l2_norm());
      SolverGMRES<TrilinosWrappers::MPI::Vector> solver_gmres(solver_F);

      tmp.reinit(src.block(1));
      // preconditionerF.vmult(dst.block(0), src.block(0));
      solver_gmres.solve(*F, dst.block(0), src.block(0), preconditionerF);

      dst.block(1) = src.block(1);
      B->vmult(dst.block(1), dst.block(0));
      dst.block(1).sadd(-1.0, src.block(1));
      tmp = dst.block(1);

      SolverControl solver_S(maxiter, tol * tmp.l2_norm());
      SolverCG<TrilinosWrappers::MPI::Vector> solver_cg(solver_S);
      solver_cg.solve(S, dst.block(1), tmp, preconditionerS);
      // preconditionerS.vmult(dst.block(1), tmp);

      dst.block(0).scale(diag_D);
      dst.block(1) *= 1.0 / alpha;
      B_T->vmult_add(dst.block(0), dst.block(1));
      dst.block(0).scale(diag_D_inv);
    }

  protected:
    const TrilinosWrappers::SparseMatrix *F;
    const TrilinosWrappers::SparseMatrix *B_T;
    const TrilinosWrappers::SparseMatrix *B;
    TrilinosWrappers::SparseMatrix S;

    TrilinosWrappers::PreconditionILU preconditionerF;
    TrilinosWrappers::PreconditionILU preconditionerS;

    TrilinosWrappers::MPI::Vector diag_D;
    TrilinosWrappers::MPI::Vector diag_D_inv;
    mutable TrilinosWrappers::MPI::Vector tmp;
    mutable TrilinosWrappers::MPI::Vector tmp2;
    const double alpha = 0.5;
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


  // g(x).
  FunctionG function_g;

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

};

#endif
