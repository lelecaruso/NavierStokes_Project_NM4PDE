#ifndef PRECONDITIONERS_HPP
#define PRECONDITIONERS_HPP
#include "IncludesFile.hpp"
using namespace dealii;

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


// Apply the SIMPLE preconditioner.
    //
    // The application of the P_SIMPLE can be divided into two steps:
    // 1. Solve the block lower triangular system:
    //      [ F    0 ] [ sol1_u ] = [ src_u ]
    //      [ B   -S ] [ sol1_p ]   [ src_p ]
    // 2. Solve the subsequent system to correct the intermediate solution:
    //      [ I    D^-1*B^T ] [ dst_u ] = [ sol1_u ]
    //      [ 0      alpha  ] [ dst_p ]   [ sol1_p ]
    //
  class PreconditionSIMPLE 
  {
  public:
    void
    initialize(const TrilinosWrappers::SparseMatrix &F_,
               const TrilinosWrappers::SparseMatrix &B_,
               const TrilinosWrappers::SparseMatrix &B_t,
               const TrilinosWrappers::MPI::BlockVector &sol_owned
               )
    {
      F = &F_;
      B = &B_; 
      B_T = &B_t;
    
      neg_diag_D_inv.reinit(sol_owned.block(0));
      diag_D_inv.reinit(sol_owned.block(0));

      for (unsigned int i : diag_D_inv.locally_owned_elements())
      {
        double temp = F->diag_element(i);
        diag_D_inv[i] = 1.0 / temp;
        neg_diag_D_inv[i] = - 1.0 / temp;     
      }

      // Create S_tilde =B * (D^-1) * B^T,
      // note: Using negative (D^-1) to create - S_tilde 
      B->mmult(negative_S_tilde, *B_T, neg_diag_D_inv); 

      // Initialize the preconditioners
      preconditioner_F.initialize(*F);
      preconditioner_S.initialize(negative_S_tilde);
      
    }
    void
    vmult(TrilinosWrappers::MPI::BlockVector &dst,
          const TrilinosWrappers::MPI::BlockVector &src) const 
    {
      const unsigned int maxiter = 10000;
      const double tol = 1e-2;
      SolverControl solver_F(maxiter, tol * src.block(0).l2_norm());

      SolverGMRES<TrilinosWrappers::MPI::Vector> solver_gmres(solver_F);

      // 1. Solve the block lower triangular system:
    //      [ F    0 ]      [ sol1_u ] =       [ src_u ]
    //      [ B   -S_tilde ] [ sol1_p ]   [ src_p ]

      //Step 1.1 Solve Fsol1_u = src_u

      // Store in temporaries the results
      TrilinosWrappers::MPI::Vector sol1_u = src.block(0);  
      TrilinosWrappers::MPI::Vector sol1_p = src.block(1);
      
      TrilinosWrappers::MPI::Vector temp_1 = src.block(1);

      solver_gmres.solve(*F, sol1_u, src.block(0), preconditioner_F);

      B->vmult(temp_1, sol1_u); //temp_1 = B * sol1_u
      temp_1 -= src.block(1); //temp1 = src_p - B * sol1_u

      //Step 1.2 Solve -S_tilde * sol1_p = src_p - temp1 (RHS)
      SolverControl solver_S(maxiter, tol * temp_1.l2_norm());
      SolverCG<TrilinosWrappers::MPI::Vector> solver_cg(solver_S);
      //Note we have already constructed S-tilde as - S_tilde 
      solver_cg.solve(negative_S_tilde, sol1_p, temp_1, preconditioner_S);

      // temp_1.reinit(dst.block(0));

      //Step 2
      // Step 2: Solve the correction system
        //         [ I   D^-1*B^T ]
        //         [ 0     alphaI  ] dst = sol1
        // =====================================================

        //2.1 scaling alpha

        dst.block(1) = sol1_p;
        dst.block(1) *= 1. / alpha; //scaling 1/alpha * sol1_p
      
        //2.2 dst(0) = sol1_u - inv(D)B.T dst(1)

        dst.block(0) = sol1_u; // Start with sol1_u.
        TrilinosWrappers::MPI::Vector tmp = src.block(0); //to have same dim as sol1_u
        B_T->vmult(tmp, dst.block(1)); ////tmp = BT*dst.block(1)
        tmp.scale(diag_D_inv); //tmp = inv(D)*tmp
        dst.block(0) -= tmp; //sol1_u - tmp
        
    }
  protected:
    const double alpha = 0.5; // parameter (0,1]

    const TrilinosWrappers::SparseMatrix *F;
    const TrilinosWrappers::SparseMatrix *B_T;
    const TrilinosWrappers::SparseMatrix *B;
    TrilinosWrappers::SparseMatrix negative_S_tilde;
    TrilinosWrappers::MPI::Vector diag_D_inv;
    TrilinosWrappers::MPI::Vector neg_diag_D_inv;
    TrilinosWrappers::PreconditionILU preconditioner_F;
    TrilinosWrappers::PreconditionILU preconditioner_S;
  };
//Simple Correct
//Approximate version
  class PreconditionaSIMPLE 
  {
  public:
    void
    initialize(const TrilinosWrappers::SparseMatrix &F_,
               const TrilinosWrappers::SparseMatrix &B_,
               const TrilinosWrappers::SparseMatrix &B_t,
               const TrilinosWrappers::MPI::BlockVector &sol_owned
               )
    {
      F = &F_;
      B = &B_;
      B_T = &B_t;
     

      diag_D_inv.reinit(sol_owned.block(0));
      diag_D.reinit(sol_owned.block(0));
      neg_diag_D_inv.reinit(sol_owned.block(0));

      for (unsigned int i : diag_D.locally_owned_elements())
      {
        double temp = F->diag_element(i);
        diag_D[i] = temp;
        diag_D_inv[i] = 1.0 / temp;
        neg_diag_D_inv[i] = - 1.0 / temp;
      }

      //S_tilde = BD(^-1)B.T
      B->mmult(neg_S, *B_T, neg_diag_D_inv); //Note: we need -S_tilde, so we use - D(^-1)

      preconditioner_F.initialize(*F);
      preconditioner_S.initialize(neg_S); //already assembled neg_S
    }

    void
    vmult(TrilinosWrappers::MPI::BlockVector &dst,
          const TrilinosWrappers::MPI::BlockVector &src) const 
    {

      const unsigned int maxit = 10000;
      const double tol = 1e-2;

      //devo definire due tmp non posso fare block 


       // Prepare 2 temporary vectors to hold intermediate data.
        tmp.reinit(src); 

        // --- Step 1 ---
        // Solve for the primary (first block) variable.
        // This computes an approximate inverse of C applied to the first part of src.
        SolverControl solver_control_F(maxit, tol * src.block(0).l2_norm());
        SolverGMRES<TrilinosWrappers::MPI::Vector> solver_F(solver_control_F);
        solver_F.solve(*F, dst.block(0), src.block(0), preconditioner_F);

        // --- Step 2 ---
        // Copy the secondary part of src into a temporary container.
        // Then update this temporary vector by incorporating the effect of -B acting
        // on the primary variable obtained in Step 1.
  
        B->vmult(dst.block(1), dst.block(0));
        dst.block(1).sadd(-1.0, src.block(1)); // tmp.block(1) = -B * dst.block(0) + src.block(1)
        tmp.block(1) = dst.block(1);

        // --- Step 3 ---
        // Solve the system with the approximate Schur complement.
        // This computes the secondary variable by inverting negS_matrix.
        SolverControl solver_control_S(maxit, tol * tmp.block(1).l2_norm());
        SolverGMRES<TrilinosWrappers::MPI::Vector> solver_S(solver_control_S);
        solver_S.solve(neg_S, dst.block(1), tmp.block(1), preconditioner_S);

        // --- Step 4 ---
        // Scale the primary component by the original diagonal entries.
        // This reintroduces the proper weighting based on C's diagonal.
        dst.block(0).scale(diag_D);

        // --- Step 5 ---
        // Adjust the secondary component by applying the damping factor.
        dst.block(1) /= alpha;

        // --- Step 6 ---
        // Refine the primary component: subtract a correction term computed by
        // multiplying the secondary variable with B^T and then scaling with the inverse
        // of the diagonal entries.
        B_T->vmult(tmp.block(0), dst.block(1));
        dst.block(0) -= tmp.block(0);
        
        // --- Step 7 ---
        // Finalize the update of the primary component by scaling it with Dinv.
        dst.block(0).scale(diag_D_inv);
      
    }

  protected:
    const TrilinosWrappers::SparseMatrix *F;
    const TrilinosWrappers::SparseMatrix *B_T;
    const TrilinosWrappers::SparseMatrix *B;
    TrilinosWrappers::SparseMatrix neg_S;

    TrilinosWrappers::PreconditionILU preconditioner_F;
    TrilinosWrappers::PreconditionILU preconditioner_S;
    

    TrilinosWrappers::MPI::Vector diag_D;
    TrilinosWrappers::MPI::Vector diag_D_inv;
    TrilinosWrappers::MPI::Vector neg_diag_D_inv;
    mutable TrilinosWrappers::MPI::BlockVector tmp;
   
    const double alpha = 1.;
  };
// approximate Simple Correct
  // Yosida preconditioner -- the inverse of Mu is replaced by the inverse of it's diagonal's elements
  class PreconditionYosida 
  {
  public:
    void
    initialize(const TrilinosWrappers::SparseMatrix &F_,
               const TrilinosWrappers::SparseMatrix &B_,
               const TrilinosWrappers::SparseMatrix &B_t,
               const TrilinosWrappers::SparseMatrix &M_,
               const TrilinosWrappers::MPI::BlockVector &sol_owned)
    {
      F = &F_;
      B = &B_;
      B_T = &B_t;
      M = &M_;
      
      diag_D_inv.reinit(sol_owned.block(0));
      neg_diag_D_inv.reinit(sol_owned.block(0));

      for (unsigned int i : diag_D_inv.locally_owned_elements())
      {
        //Note : we have assembled M as M/deltat
        diag_D_inv[i] = ( 1.0 / M->diag_element(i));  //  dt * (Mii)^-1
        neg_diag_D_inv[i] = ( -1.0 / M->diag_element(i));  //  dt * (Mii)^-1
      }

      // Create negative_S_tilde
      B->mmult(negative_S_tilde, *B_T, neg_diag_D_inv);
    
      // Initialize the preconditioners
      preconditioner_F.initialize(*F);
      preconditioner_S.initialize(negative_S_tilde);
    }
    void
    vmult(TrilinosWrappers::MPI::BlockVector &dst,
          const TrilinosWrappers::MPI::BlockVector &src) const 
    {
      const unsigned int maxiter = 100000;
      const double tol = 1e-2;

      SolverControl solver_F(maxiter, tol * src.block(0).l2_norm());
      SolverGMRES<TrilinosWrappers::MPI::Vector> solver_gmres(solver_F);

      // Store in temporaries the results
      TrilinosWrappers::MPI::Vector yu = src.block(0);
      TrilinosWrappers::MPI::Vector yp = src.block(1);
      TrilinosWrappers::MPI::Vector tmp = src.block(1);
      TrilinosWrappers::MPI::Vector tmp2 = src.block(0);

      //Step 1
      // Step 1.1) yu = F^-1 * src.0
      solver_gmres.solve(*F, yu, src.block(0), preconditioner_F);
      
      //Step 1.2) yp = negative_S_tilde^-1(src1-B*yu)
      B->vmult(tmp, yu); //tmp = B*yu
      tmp.add(-1.0, src.block(1)); // tmp = src.block(1) - tmp
      // neg_S*yp = (src(1) - Byu)==tmp(RHS)
      SolverControl solver_S(maxiter, tol * tmp.l2_norm());
      SolverCG<TrilinosWrappers::MPI::Vector> solver_cg(solver_S);
      solver_cg.solve(negative_S_tilde, yp, tmp, preconditioner_S);

      //Step 2) 
      // Step 2.1) dst1 = yp
      dst.block(1) = yp; 

      // Step 2.2) dst0 = yu - F^-1*B_T*yp
        //Step 2.2.1) 
      B_T->vmult(tmp2, dst.block(1)); //tmp2 = B_T*yp (rhs)

      //Solve the linear system 
      res.reinit(src.block(0)); //to store the result of the  lin sys F res = tmp2
      dst.block(0) = yu; //init final velocity dest 
      SolverControl solver_F2(maxiter, tol * tmp2.l2_norm());
      SolverGMRES<TrilinosWrappers::MPI::Vector> solver_gmres2(solver_F2);
      solver_gmres2.solve(*F, res, tmp2, preconditioner_F); // res = F^-1 * tmp2 
      dst.block(0).sadd(-1,res); //update final velocity dest dstu = yu - res

    }

  protected:

    const TrilinosWrappers::SparseMatrix *F;
    const TrilinosWrappers::SparseMatrix *B_T;
    const TrilinosWrappers::SparseMatrix *B;
    const TrilinosWrappers::SparseMatrix *M;
    TrilinosWrappers::SparseMatrix negative_S_tilde;
    TrilinosWrappers::MPI::Vector diag_D_inv;
    TrilinosWrappers::MPI::Vector neg_diag_D_inv;
    TrilinosWrappers::PreconditionILU preconditioner_F;
    TrilinosWrappers::PreconditionILU preconditioner_S;

    mutable TrilinosWrappers::MPI::Vector res;
  };

  // Precondition approximate Yosida: Why it so slow?
  
  class PreconditionaYosida 
 {
  public:
    void
    initialize(const TrilinosWrappers::SparseMatrix &F_,
               const TrilinosWrappers::SparseMatrix &B_,
               const TrilinosWrappers::SparseMatrix &B_t,
               const TrilinosWrappers::SparseMatrix &M_,
               const TrilinosWrappers::MPI::BlockVector &sol_owned)
    {
      F = &F_;
      B = &B_;
      B_T = &B_t;
      M = &M_;

      diag_D_inv.reinit(sol_owned.block(0));
      diag_D.reinit(sol_owned.block(0));
      lump_M.reinit(sol_owned.block(0));


      for (unsigned int i : diag_D.locally_owned_elements())
      {
        double temp = F->diag_element(i);
        diag_D[i] = temp;
        diag_D_inv[i] = 1.0 / temp; // (F_hat)^-1
      }


    // build lumped_M
    for (unsigned int i : lump_M.locally_owned_elements())
      {
        //Note : we have assembled M/deltat
        double temp = 0.0;
        for (auto it = M->begin(i); it != M->end(i); ++it)
            temp += std::abs(it->value());

        
        lump_M[i] = -1.0 / temp; // - deltat * (lump_M)^-1
      }

      //Note: We use -lump_M to create negative S
      B->mmult(negative_S, *B_T, lump_M); // neg_S

      preconditionerF.initialize(*F);
      preconditionerS.initialize(negative_S);
    }

    void
    vmult(TrilinosWrappers::MPI::BlockVector &dst,
          const TrilinosWrappers::MPI::BlockVector &src) const 
    { 
      //Note : diag_D_inv = (F_hat)^-1
      tmp.reinit(src.block(0)); // block 0
      tmp2.reinit(src.block(1)); //block 1 

      const unsigned int maxiter = 100000;
      const double tol = 1e-2;

      // Store in temporaries the results of src updates 
      TrilinosWrappers::MPI::Vector yu = src.block(0);
      TrilinosWrappers::MPI::Vector yp = src.block(1);

      //Step 1) (F_hat)^-1 * src(0) = tmp
      
      tmp = src.block(0);
      tmp.scale(diag_D_inv);
      yu = tmp;      

      //Step 2)   
       B->vmult(tmp2, tmp); //tmp(1) = B*tmp(0)
       yp.sadd(-1.0,tmp2); // src(1) = -tmp(1) + src(1) (RHS)
       
       //Step 3) true solution of neg_S to have better accuracy, instead of neg_S_hat
      SolverControl solver_S(maxiter, tol * yp.l2_norm());
      SolverCG<TrilinosWrappers::MPI::Vector> solver_cg(solver_S);
      solver_cg.solve(negative_S, dst.block(1), yp, preconditionerS); //dst.block(1) updated here 

      yp = dst.block(1); //updating src(1) for next computations

      //Step 4) 
      F->vmult(yu,yu);

      //Step 5)
       B_T->vmult(tmp, yp); //tmp(0) = BT*src(1)
       yu.sadd(-1.0,tmp);

      //Step 6) 
      yu.scale(diag_D_inv);
      dst.block(0) = yu;

    }

  protected:
    const TrilinosWrappers::SparseMatrix *F;
    const TrilinosWrappers::SparseMatrix *B_T;
    const TrilinosWrappers::SparseMatrix *B;
    const TrilinosWrappers::SparseMatrix *M;
    TrilinosWrappers::SparseMatrix negative_S;

    TrilinosWrappers::PreconditionILU preconditionerF;
    TrilinosWrappers::PreconditionILU preconditionerS;

    TrilinosWrappers::MPI::Vector diag_D;
    TrilinosWrappers::MPI::Vector lump_M;
    TrilinosWrappers::MPI::Vector diag_D_inv;
    mutable TrilinosWrappers::MPI::Vector tmp;
    mutable TrilinosWrappers::MPI::Vector tmp2;
  }; 
  
  #endif
  