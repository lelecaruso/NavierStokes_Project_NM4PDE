#include "../include/NavierStokes2D.hpp"
void NavierStokes::setup()
{
  // Create the mesh.
  {
    pcout << "Initializing the mesh" << std::endl;

    Triangulation<dim> mesh_serial;

    GridIn<dim> grid_in;
    grid_in.attach_triangulation(mesh_serial);

    std::ifstream grid_in_file(mesh_file_name);
    grid_in.read_msh(grid_in_file);

    GridTools::partition_triangulation(mpi_size, mesh_serial);
    const auto construction_data = TriangulationDescription::Utilities::
        create_description_from_triangulation(mesh_serial, MPI_COMM_WORLD);
    mesh.create_triangulation(construction_data);

    pcout << "  Number of elements = " << mesh.n_global_active_cells()
          << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the finite element space.
  {
    pcout << "Initializing the finite element space" << std::endl;

    const FE_SimplexP<dim> fe_scalar_velocity(degree_velocity);
    const FE_SimplexP<dim> fe_scalar_pressure(degree_pressure);
    fe = std::make_unique<FESystem<dim>>(fe_scalar_velocity,
                                         dim,
                                         fe_scalar_pressure,
                                         1);

    pcout << "  Velocity degree:           = " << fe_scalar_velocity.degree
          << std::endl;
    pcout << "  Pressure degree:           = " << fe_scalar_pressure.degree
          << std::endl;
    pcout << "  DoFs per cell              = " << fe->dofs_per_cell
          << std::endl;

    quadrature = std::make_unique<QGaussSimplex<dim>>(fe->degree + 1);

    pcout << "  Quadrature points per cell = " << quadrature->size()
          << std::endl;

    quadrature_boundary = std::make_unique<QGaussSimplex<dim - 1>>(fe->degree + 1);

    pcout << "  Quadrature points per boundary cell = " << quadrature_boundary->size()
          << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the DoF handler.
  {
    pcout << "Initializing the DoF handler" << std::endl;

    dof_handler.reinit(mesh);
    dof_handler.distribute_dofs(*fe);

    // We want to reorder DoFs so that all velocity DoFs come first, and then
    // all pressure DoFs.
    std::vector<unsigned int> block_component(dim + 1, 0);
    block_component[dim] = 1;
    DoFRenumbering::component_wise(dof_handler, block_component);

    locally_owned_dofs = dof_handler.locally_owned_dofs();
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

    // Besides the locally owned and locally relevant indices for the whole
    // system (velocity and pressure), we will also need those for the
    // individual velocity and pressure blocks.
    std::vector<types::global_dof_index> dofs_per_block =
        DoFTools::count_dofs_per_fe_block(dof_handler, block_component);
    const unsigned int n_u = dofs_per_block[0];
    const unsigned int n_p = dofs_per_block[1];

    block_owned_dofs.resize(2);
    block_relevant_dofs.resize(2);
    block_owned_dofs[0] = locally_owned_dofs.get_view(0, n_u);
    block_owned_dofs[1] = locally_owned_dofs.get_view(n_u, n_u + n_p);
    block_relevant_dofs[0] = locally_relevant_dofs.get_view(0, n_u);
    block_relevant_dofs[1] = locally_relevant_dofs.get_view(n_u, n_u + n_p);

    pcout << "  Number of DoFs: " << std::endl;
    pcout << "    velocity = " << n_u << std::endl;
    pcout << "    pressure = " << n_p << std::endl;
    pcout << "    total    = " << n_u + n_p << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the linear system.
  {
    pcout << "Initializing the linear system" << std::endl;

    pcout << "  Initializing the sparsity pattern" << std::endl;

    // Velocity DoFs interact with other velocity DoFs (the weak formulation has
    // terms involving u times v), and pressure DoFs interact with velocity DoFs
    // (there are terms involving p times v or u times q). However, pressure
    // DoFs do not interact with other pressure DoFs (there are no terms
    // involving p times q). We build a table to store this information, so that
    // the sparsity pattern can be built accordingly.
    Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1);
    for (unsigned int c = 0; c < dim + 1; ++c)
    {
      for (unsigned int d = 0; d < dim + 1; ++d)
      {
        if (c == dim && d == dim) // pressure-pressure term
          coupling[c][d] = DoFTools::none;
        else // other combinations
          coupling[c][d] = DoFTools::always;
      }
    }

    TrilinosWrappers::BlockSparsityPattern sparsity(block_owned_dofs,
                                                    MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern(dof_handler, coupling, sparsity);
    sparsity.compress();

    // We also build a sparsity pattern for the pressure mass matrix.
    for (unsigned int c = 0; c < dim + 1; ++c)
    {
      for (unsigned int d = 0; d < dim + 1; ++d)
      {
        if (c == dim && d == dim) // pressure-pressure term
          coupling[c][d] = DoFTools::always;
        else // other combinations
          coupling[c][d] = DoFTools::none;
      }
    }
    TrilinosWrappers::BlockSparsityPattern sparsity_pressure_mass(
        block_owned_dofs, MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern(dof_handler,
                                    coupling,
                                    sparsity_pressure_mass);
    sparsity_pressure_mass.compress();

    pcout << "  Initializing the matrices" << std::endl;
    system_matrix.reinit(sparsity);
    mass_matrix.reinit(sparsity);
    convection_matrix.reinit(sparsity);
    stiffness_matrix.reinit(sparsity);
    pressure_mass.reinit(sparsity_pressure_mass);

    pcout << "  Initializing the system right-hand side" << std::endl;
    system_rhs.reinit(block_owned_dofs, MPI_COMM_WORLD);
    pcout << "  Initializing the solution vector" << std::endl;
    solution_owned.reinit(block_owned_dofs, MPI_COMM_WORLD);
    solution.reinit(block_owned_dofs, block_relevant_dofs, MPI_COMM_WORLD);
  }
}




// Function used to assemble the static Matrixes, 
// mass matrix, the stiffness matrix, the pressure matrix
void NavierStokes::assemble(const double &time)
{
  pcout << "===============================================" << std::endl;
  pcout << "Assembling the system" << std::endl;

  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q = quadrature->size();
  const unsigned int n_q_boundary = quadrature_boundary->size();

  FEValues<dim> fe_values(*fe,
                          *quadrature,
                          update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);
  FEFaceValues<dim> fe_boundary_values(*fe,
                                       *quadrature_boundary,
                                       update_values | update_quadrature_points |
                                           update_normal_vectors |
                                           update_JxW_values);

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> cell_mass_matrix(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> cell_stiffness_matrix(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> cell_convection_matrix(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> cell_pressure_mass_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double> cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  system_matrix = 0.0;
  mass_matrix = 0.0;
  stiffness_matrix = 0.0;
  convection_matrix = 0.0;
  system_rhs = 0.0;
  pressure_mass = 0.0;

  FEValuesExtractors::Vector velocity(0);
  FEValuesExtractors::Scalar pressure(dim);

  // Store the current velocity value 
  std::vector<Tensor<1, dim>> current_velocity_values(n_q);
  //Store the current velocity gradient value 
  std::vector<Tensor<2,dim>> current_velocity_gradients(n_q);
  // Store the current velocity divergence value 
  std::vector<double> current_velocity_divergence(n_q);

  for (const auto &cell : dof_handler.active_cell_iterators())
  {

    if (!cell->is_locally_owned())
      continue;

    fe_values.reinit(cell);

    cell_matrix = 0.0;
    cell_mass_matrix = 0.0;
    cell_stiffness_matrix = 0.0;
    cell_convection_matrix = 0.0;
    cell_rhs = 0.0;
    cell_pressure_mass_matrix = 0.0;

    // Retrieve the current solution values.
    fe_values[velocity].get_function_values(solution, current_velocity_values);
    // Retrieve the current solution gradient values
    fe_values[velocity].get_function_gradients(solution, current_velocity_gradients);
    // Retrieve the current solution divergence values
    fe_values[velocity].get_function_divergences(solution, current_velocity_divergence);
    
    for (unsigned int q = 0; q < n_q; ++q)
    {
      
      Vector<double> forcing_term_loc(dim);
      forcing_term.vector_value(fe_values.quadrature_point(q),
                                forcing_term_loc);
      Tensor<1, dim> forcing_term_tensor;
      for (unsigned int d = 0; d < dim; ++d)
        forcing_term_tensor[d] = forcing_term_loc[d];

      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
        {
          
          // Viscosity term.
          cell_stiffness_matrix(i, j) += nu * scalar_product(fe_values[velocity].gradient(i, q), fe_values[velocity].gradient(j, q)) * fe_values.JxW(q);

          // Time derivative discretization.
          cell_mass_matrix(i, j) +=  scalar_product(fe_values[velocity].value(i, q), fe_values[velocity].value(j, q)) / deltat * fe_values.JxW(q);

          // Convective term 
          cell_convection_matrix(i, j) += scalar_product(fe_values[velocity].gradient(j, q) * current_velocity_values[q], fe_values[velocity].value(i, q)) * fe_values.JxW(q);
          
          // Temam Stabilization term
          cell_convection_matrix(i, j) += 0.5 * current_velocity_divergence[q] * scalar_product(fe_values[velocity].value(i, q), fe_values[velocity].value(j, q)) * fe_values.JxW(q);             
          
          // Pressure term in the momentum equation.
          cell_matrix(i, j) -= fe_values[pressure].value(j, q) * fe_values[velocity].divergence(i, q) * fe_values.JxW(q);

          // Pressure term in the continuity equation.
          cell_matrix(i, j) += fe_values[pressure].value(i, q) * fe_values[velocity].divergence(j, q) * fe_values.JxW(q);

          // Pressure mass matrix.
          cell_pressure_mass_matrix(i, j) += fe_values[pressure].value(i, q) * fe_values[pressure].value(j, q) / nu * fe_values.JxW(q);

        }

        // Time derivative discretization on the right hand side
        cell_rhs(i) +=  scalar_product(current_velocity_values[q], fe_values[velocity].value(i, q)) * fe_values.JxW(q) / deltat;

      }
    }

    // Boundary integral for Neumann BCs.
    if (cell->at_boundary())
    {
      for (unsigned int f = 0; f < cell->n_faces(); ++f)
      {
        // NO NEUMANN
        if (cell->face(f)->at_boundary() &&
            (false))
        {
          fe_boundary_values.reinit(cell, f);

          for (unsigned int q = 0; q < n_q_boundary; ++q)
          {
            Vector<double> neumann_loc(dim);
            function_h.vector_value(fe_boundary_values.quadrature_point(q),
                                    neumann_loc);
            Tensor<1, dim> neumann_loc_tensor;
            for (unsigned int d = 0; d < dim; ++d)
              neumann_loc_tensor[d] = neumann_loc[d];

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              cell_rhs(i) +=
                  scalar_product(neumann_loc_tensor, fe_boundary_values[velocity].value(i, q)) * fe_boundary_values.JxW(q);
            }
          }
        }
      }
    }

    cell->get_dof_indices(dof_indices);

    system_matrix.add(dof_indices, cell_matrix);
    mass_matrix.add(dof_indices, cell_mass_matrix);
    convection_matrix.add(dof_indices, cell_convection_matrix);
    stiffness_matrix.add(dof_indices, cell_stiffness_matrix);
    system_rhs.add(dof_indices, cell_rhs);
    pressure_mass.add(dof_indices, cell_pressure_mass_matrix);
  }

  system_matrix.compress(VectorOperation::add);
  mass_matrix.compress(VectorOperation::add);
  convection_matrix.compress(VectorOperation::add);
  stiffness_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);
  pressure_mass.compress(VectorOperation::add);

  // Create the System Matrix F = M + A + C(u_n) + B
  system_matrix.add(1., mass_matrix);
  system_matrix.add(1., convection_matrix);
  system_matrix.add(1., stiffness_matrix);

  // Apply Dirichlet boundary conditions.
  {
    std::map<types::global_dof_index, double> boundary_values;
    std::map<types::boundary_id, const Function<dim> *> boundary_functions;

    // We first impose the Dirichlet boundary conditions on the Inlet.
    inlet_velocity.set_time(time);
    boundary_functions[0] = &inlet_velocity;
    VectorTools::interpolate_boundary_values(dof_handler,
                                            boundary_functions,
                                            boundary_values,
                                            ComponentMask(
                                                {true, true, false}));

    // This ensure the two boundaries do not overlap.
    boundary_functions.clear();
    Functions::ZeroFunction<dim> zero_function(dim + 1);
    
    // We then impose the Dirichlet boundary conditions on Walls and the Obstacle.
    boundary_functions[2] = &zero_function;
    boundary_functions[3] = &zero_function;
    VectorTools::interpolate_boundary_values(dof_handler,
                                            boundary_functions,
                                            boundary_values,
                                            ComponentMask(
                                                {true, true, false}));

    MatrixTools::apply_boundary_values(boundary_values, system_matrix, solution, system_rhs, false);
  }

}

// Function used to assemble at time > deltat to avoid redundant computation of A,M,B
// assemble rhs and convection matrix
void NavierStokes::assemble_time_step(const double &time)
{
  pcout << "===============================================" << std::endl;
  pcout << "Assembling the system" << std::endl;

  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q = quadrature->size();
  const unsigned int n_q_boundary = quadrature_boundary->size();

  FEValues<dim> fe_values(*fe,
                          *quadrature,
                          update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);
  FEFaceValues<dim> fe_boundary_values(*fe,
                                       *quadrature_boundary,
                                       update_values | update_quadrature_points |
                                           update_normal_vectors |
                                           update_JxW_values);


  FullMatrix<double> cell_convection_matrix(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> cell_mass_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double> cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  // We delete the previous Convection Matrix from the system matrix 
  system_matrix.add(-1., convection_matrix);
  convection_matrix = 0.0;
  system_rhs = 0.0;

  FEValuesExtractors::Vector velocity(0);
  FEValuesExtractors::Scalar pressure(dim);
  std::vector<Tensor<1, dim>> boundary_velocity_values(n_q_boundary);
  std::vector<Tensor<1, dim>> prev_boundary_velocity_values(n_q_boundary);

  // Store the current velocity value in a tensor
  std::vector<Tensor<1, dim>> current_velocity_values(n_q);
  //Store the current velocity gradient value in a tensor
  std::vector<Tensor<2,dim>> current_velocity_gradients(n_q);
  // Store the current velocity divergence value in a tensor
  std::vector<double> current_velocity_divergence(n_q);
  // Store the previous velocity value in a tensor
  std::vector<double> prev_velocity_diverg(n_q);

  // Store the current velocity value in a tensor
  std::vector<Tensor<1, dim>> prev_velocity_values(n_q);
  //Store the prev velocity gradient value in a tensor
  std::vector<Tensor<2,dim>> prev_velocity_gradients(n_q);
  // Store the prev velocity divergence value in a tensor
  std::vector<double> prev_velocity_divergence(n_q);
  

  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    if (!cell->is_locally_owned())
      continue;

    fe_values.reinit(cell);

    cell_mass_matrix = 0.0;
    cell_convection_matrix = 0.0;
    cell_rhs = 0.0;

    // Retrieve the previous solution values.
    fe_values[velocity].get_function_values(previous_solution, prev_velocity_values);
    // Retrieve the previous solution gradient values
    fe_values[velocity].get_function_divergences(previous_solution, prev_velocity_diverg);
    // Retrieve the current solution values.
    fe_values[velocity].get_function_values(solution, current_velocity_values);
    //Retrieve the current solution gradient values
    fe_values[velocity].get_function_gradients(solution, current_velocity_gradients);
    // Retrieve the current solution divergence values
    fe_values[velocity].get_function_divergences(solution, current_velocity_divergence);

    for (unsigned int q = 0; q < n_q; ++q)
    {
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
        {

          // Convective term 
          cell_convection_matrix(i, j) += scalar_product(fe_values[velocity].gradient(j, q) * current_velocity_values[q], fe_values[velocity].value(i, q)) * fe_values.JxW(q);
          // Tamam Stabilization term 0.5 = rho / 2
          cell_convection_matrix(i, j) += 0.5 * current_velocity_divergence[q] * scalar_product(fe_values[velocity].value(i, q), fe_values[velocity].value(j, q)) * fe_values.JxW(q);

        }
        // Time derivative discretization on the right hand side BDF2
        cell_rhs(i) +=  scalar_product(current_velocity_values[q], fe_values[velocity].value(i, q)) * fe_values.JxW(q) / deltat;


      }
    }

    // BackFlow Stabilization on open boundary ( mainly for 3D instabilities in turbolent flows )
    if (cell->at_boundary())
    {
      for (unsigned int f = 0; f < cell->n_faces(); ++f)
      {
        
        if (cell->face(f)->at_boundary() &&
            (false))
        {
          fe_boundary_values.reinit(cell, f);
          fe_boundary_values[velocity].get_function_values(solution, boundary_velocity_values);
          fe_boundary_values[velocity].get_function_values(previous_solution, prev_boundary_velocity_values);

          for (unsigned int q = 0; q < n_q_boundary; ++q)
          {
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
              {
                  cell_convection_matrix(i, j) -= 1.5 * std::min((2. * boundary_velocity_values[q] - prev_boundary_velocity_values[q] )* fe_boundary_values.normal_vector(q), 0.0) *
                                                        scalar_product(fe_boundary_values[velocity].value(j, q),fe_boundary_values[velocity].value(i, q)) *
                                                        fe_boundary_values.JxW(q);;
              }
            }
          }
        }
      }
    }

    cell->get_dof_indices(dof_indices);
    convection_matrix.add(dof_indices, cell_convection_matrix);
    system_rhs.add(dof_indices, cell_rhs);
  }
  convection_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);
  pressure_mass.compress(VectorOperation::add);
  system_matrix.add(1., convection_matrix);


  // Dirichlet boundary conditions.

  // Apply Dirichlet boundary conditions.
  {
    std::map<types::global_dof_index, double> boundary_values;
    std::map<types::boundary_id, const Function<dim> *> boundary_functions;

    // We first impose the Dirichlet boundary conditions on the Inlet.
    inlet_velocity.set_time(time);
    boundary_functions[0] = &inlet_velocity;
    VectorTools::interpolate_boundary_values(dof_handler,
                                            boundary_functions,
                                            boundary_values,
                                            ComponentMask(
                                                {true, true, false}));

    // This ensure the two boundaries do not overlap.
    boundary_functions.clear();
    Functions::ZeroFunction<dim> zero_function(dim + 1);
    
    // We then impose the Dirichlet boundary conditions on Walls and the Obstacle.
    boundary_functions[2] = &zero_function;
    boundary_functions[3] = &zero_function;
    VectorTools::interpolate_boundary_values(dof_handler,
                                            boundary_functions,
                                            boundary_values,
                                            ComponentMask(
                                                {true, true, false}));

    MatrixTools::apply_boundary_values(boundary_values, system_matrix, solution, system_rhs, false);
  }

}

// Function used to solve the linear system and assemble the preconditioner
void NavierStokes::solve_time_step(double time)
{
  pcout << "===============================================" << std::endl;

  const unsigned int maxiter = 100000;
  const double tol = 1e-4 /**system_rhs.l2_norm()*/;
  SolverControl solver_control(maxiter, tol, true);
  // solver_control.enable_history_data();
  SolverGMRES<TrilinosWrappers::MPI::BlockVector> solver(solver_control);
  previous_solution = solution;
  // Assemblying the preconditioner
  {

    dealii::Timer timerprec;
    timerprec.restart();
    dealii::Timer timersys;
    
    unsigned int preconditioner_type = 3;
    switch (preconditioner_type)
    {
        // Yosida
        case 0:
        {
            PreconditionYosida yosida;
            yosida.initialize(system_matrix.block(0, 0), system_matrix.block(1, 0), system_matrix.block(0, 1), mass_matrix.block(0, 0), solution_owned);  // Yosida
            timerprec.stop();
            pcout << "Time taken to initialize preconditioner: " << timerprec.wall_time() << " seconds" << std::endl;
            time_prec.push_back(timerprec.wall_time());
            timersys.restart();
            solver.solve(system_matrix, solution_owned, system_rhs, yosida);
            timersys.stop();
            pcout << "Time taken to solve Navier Stokes problem: " << timersys.wall_time() << " seconds" << std::endl;
            time_solve.push_back(timersys.wall_time());
            break;
        }
    
        // SIMPLE
        case 1:
        {
            PreconditionSIMPLE simple;
            simple.initialize(system_matrix.block(0, 0), system_matrix.block(1, 0), system_matrix.block(0, 1), solution_owned);
            timerprec.stop();
            pcout << "Time taken to initialize preconditioner: " << timerprec.wall_time() << " seconds" << std::endl;
            time_prec.push_back(timerprec.wall_time());
            timersys.restart();
            solver.solve(system_matrix, solution_owned, system_rhs, simple);
            timersys.stop();
            pcout << "Time taken to solve Navier Stokes problem: " << timersys.wall_time() << " seconds" << std::endl;
            time_solve.push_back(timersys.wall_time());
            
            break;
        }
    
        // aYosida
        case 2:
        {
            PreconditionaYosida ayosida;
            ayosida.initialize(system_matrix.block(0, 0), system_matrix.block(1, 0), system_matrix.block(0, 1), mass_matrix.block(0, 0), solution_owned);  // Yosida
            timerprec.stop();
            pcout << "Time taken to initialize preconditioner: " << timerprec.wall_time() << " seconds" << std::endl;
            time_prec.push_back(timerprec.wall_time());
            timersys.restart();
            solver.solve(system_matrix, solution_owned, system_rhs, ayosida);
            timersys.stop();
            pcout << "Time taken to solve Navier Stokes problem: " << timersys.wall_time() << " seconds" << std::endl;
            time_solve.push_back(timersys.wall_time());
            
            break;
        }
    
        // aSIMPLE
        case 3:
        {
            PreconditionaSIMPLE asimple;
            asimple.initialize(system_matrix.block(0, 0), system_matrix.block(1, 0), system_matrix.block(0, 1), solution_owned);
            timerprec.stop();
            pcout << "Time taken to initialize preconditioner: " << timerprec.wall_time() << " seconds" << std::endl;
            time_prec.push_back(timerprec.wall_time());
            timersys.restart();
            solver.solve(system_matrix, solution_owned, system_rhs, asimple);
            timersys.stop();
            pcout << "Time taken to solve Navier Stokes problem: " << timersys.wall_time() << " seconds" << std::endl;
            time_solve.push_back(timersys.wall_time());
            
            break;
        }
    
        default:
            throw std::runtime_error("Invalid preconditioner type");
    }
  }
  pcout << "Result:  " << solver_control.last_step() << " GMRES iterations"<< std::endl;
  int Re = int(0.1 * 1.5 * std::sin(time*M_PI/8.0) / .001);
      // Write coefficients to "coeff.csv"
    if (mpi_rank == 0) // Ensure only the root process writes to the file
    {
        std::ofstream coeff_file("gmres.csv", std::ios::app); // Open in append mode
        if (coeff_file.is_open())
        {
            coeff_file << time << ',' << Re << ',' << solver_control.last_step() << "\n";
            coeff_file.close();
        }
        else
        {
            pcout << "Error: Unable to open coeff.csv for writing." << std::endl;
        }
    }
  solution = solution_owned;

}

// Function used to save the output of the simulation
void NavierStokes::output(const unsigned int &time_step, std::vector<double> coeff) const
{
    pcout << "===============================================" << std::endl;

    DataOut<dim> data_out;

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
        data_component_interpretation(
            dim, DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretation.push_back(
        DataComponentInterpretation::component_is_scalar);
    std::vector<std::string> names = {"velocity",
                                      "velocity",
                                      "pressure"};

    data_out.add_data_vector(dof_handler,
                            solution,
                            names,
                            data_component_interpretation);

    std::vector<unsigned int> partition_int(mesh.n_active_cells());
    GridTools::get_subdomain_association(mesh, partition_int);
    const Vector<double> partitioning(partition_int.begin(), partition_int.end());
    data_out.add_data_vector(partitioning, "partitioning");

    data_out.build_patches();

    const std::string output_file_name = "output-navier-stokes-2D";
    data_out.write_vtu_with_pvtu_record("./output2D_1/",
                                        output_file_name,
                                        time_step,
                                        MPI_COMM_WORLD,
                                        numbers::invalid_unsigned_int,
                                        1);

    pcout << "Output written to " << output_file_name << std::endl;

    // Write coefficients to "coeff.csv"
    if (mpi_rank == 0) // Ensure only the root process writes to the file
    {
        std::ofstream coeff_file("coeff_2.csv", std::ios::app); // Open in append mode
        if (coeff_file.is_open())
        {
            coeff_file << time_step << "," << coeff[0] << "," << coeff[1] << "\n";
            coeff_file.close();
        }
        else
        {
            pcout << "Error: Unable to open coeff.csv for writing." << std::endl;
        }
    }

    pcout << "===============================================" << std::endl;    
}


// Function used to update time step and call the solver, compute the forces and output the results
void NavierStokes::solve()
{
  
  pcout << "===============================================" << std::endl;

  // Apply the initial condition.
  {
    pcout << "Applying the initial condition" << std::endl;

    VectorTools::interpolate(dof_handler, u_0, solution_owned);
    solution = solution_owned;

    // Output the initial solution.
    output(0, {0.0, 0.0});
    pcout << "===============================================" << std::endl;
  }
  std::vector<double> coefficients;
  double c_D_max = -999;
  double c_L_min = 999;
  unsigned int time_step = 0;
  double time = 0;
  while (time < T - 0.5 * deltat)
  { 

    time += deltat;
    ++time_step;
    inlet_velocity.set_time(time);

    pcout << "n = " << std::setw(3) << time_step << ", t = " << std::setw(5)
          << time << ":" << std::flush;


    if( time == deltat ) assemble(time);
    else assemble_time_step(time);

    solve_time_step(time);
    if( time == T - deltat )
        compute_pressure_difference();
    {
      coefficients = compute_forces();
      c_D_max = std::max(c_D_max, coefficients[0]);
      c_L_min = std::min(c_L_min, coefficients[1]);
    }
    // Since the starting solution t0 is zero we avoid the initial high forces values
    if( time_step % 1 == 0) output(time_step, coefficients);
  }
  pcout << "===============================================" << std::endl;
  pcout << "Drag Coefficient Max ----->   " << c_D_max << std::endl;
  pcout << std::endl;
  pcout << "Lift Coefficient Min ----->   " << c_L_min << std::endl;
  pcout << "===============================================" << std::endl;
}

std::vector<double> NavierStokes::compute_forces()
{
  pcout << "===============================================" << std::endl;
  pcout << "Computing forces: " << std::endl;

   // Define quadrature for faces
   QGauss<dim - 1> face_quadrature_formula(3);
   const unsigned int n_q_points = face_quadrature_formula.size();

   // Define FE extractors for velocity and pressure
   FEValuesExtractors::Vector velocities(0);
   FEValuesExtractors::Scalar pressure(dim);

   // Containers to store evaluated values
   std::vector<double>         pressure_values(n_q_points);
   std::vector<Tensor<2, dim>> velocity_gradients(n_q_points);

   // Initialize FE face values
   FEFaceValues<dim> fe_face_values(*this->fe,
                                    face_quadrature_formula,
                                    update_values | update_quadrature_points |
                                    update_gradients | update_JxW_values |
                                    update_normal_vectors);

   // Initialize tensors for calculations
   Tensor<1, dim> normal_vector;
   Tensor<2, dim> fluid_stress;
   Tensor<2, dim> fluid_pressure;
   Tensor<1, dim> forces;

   // Initialize drag and lift accumulators
   double local_drag = 0.0;
   double local_lift = 0.0;

   // Iterate over all cells
   for (const auto &cell : this->dof_handler.active_cell_iterators())
   {
       if (!cell->is_locally_owned())
           continue;
           
       if (!cell->at_boundary())
           continue;

       for (unsigned int f = 0; f < cell->n_faces(); ++f)
       {
           if (!cell->face(f)->at_boundary())
               continue;

           unsigned int boundary_id = cell->face(f)->boundary_id();
           bool is_stress_boundary = false;

           // Determine if current face is where stress should be evaluated           
           if (boundary_id == 3) // Obstacle 
               is_stress_boundary = true;

           if (!is_stress_boundary)
               continue;

           // Reinitialize FE face values for the current face
           fe_face_values.reinit(cell, f);

           // Retrieve velocity gradients and pressure values on the face
           fe_face_values[velocities].get_function_gradients(this->solution, velocity_gradients);
           fe_face_values[pressure].get_function_values(this->solution, pressure_values);

           // Iterate over quadrature points on the face
           for (unsigned int q = 0; q < n_q_points; ++q)
           {
               // Get the normal vector
               normal_vector = -fe_face_values.normal_vector(q);

               // Compute fluid pressure tensor (p * I)
               fluid_pressure = 0;
               for (unsigned int d = 0; d < dim; ++d)
                   fluid_pressure[d][d] = pressure_values[q];

               // Compute fluid stress tensor (nu * grad(U) - pI)
               fluid_stress = this->nu * velocity_gradients[q] - fluid_pressure;

               // Compute forces: stress tensor contracted with normal vector and scaled by JxW
               forces = fluid_stress * normal_vector * fe_face_values.JxW(q);

               // Accumulate drag and lift using the scaling factor
               local_drag += forces[0];
               local_lift += forces[1];
           }
       }
   }

  
   // Reduce lift and drag across all processes to rank 0
   double total_lift = 0.0;
   double total_drag = 0.0;

   MPI_Reduce(&local_lift, &total_lift, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(&local_drag, &total_drag, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  const double mean_v = inlet_velocity.getMeanVelocity();
	const double D= 0.1;
	

	const double c_d=(2.*total_drag)/(mean_v*mean_v*D);
	const double c_l=(2.*total_lift)/(mean_v*mean_v*D);
  std::vector<double> coefficients = { c_d , c_l };
	pcout << "Coeff:\t " << c_d << " Coeff:\t " << c_l << std::endl;

  pcout << "===============================================" << std::endl;
  return coefficients;
}


void NavierStokes::compute_pressure_difference()
{
  Point<dim> p_a = { 0.45, 0.2 , 0.205 };
  Point<dim> p_e = { 0.55, 0.2 , 0.205 };

  Vector<double> solution_values1(dim + 1);
  Vector<double> solution_values2(dim + 1);


  bool p1_available = true;
  bool p2_available = true;

  // Attempt to evaluate pressure at p1
  try
  {
      VectorTools::point_value(this->dof_handler, this->solution, p_a,
                               solution_values1);
  }
  catch (const dealii::VectorTools::ExcPointNotAvailableHere &)
  {
      p1_available = false;
  }

  // Attempt to evaluate pressure at p2
  try
  {
      VectorTools::point_value(this->dof_handler, this->solution, p_e,
                               solution_values2);
  }
  catch (const dealii::VectorTools::ExcPointNotAvailableHere &)
  {
      p2_available = false;
  }

    // Initialize pressure variables
    double pres_point1 = 0.0;
    double pres_point2 = 0.0;

    // Assign local pressure values if available
    if (p1_available)
        pres_point1 = solution_values1(dim);
    if (p2_available)
        pres_point2 = solution_values2(dim);

    // Reduce pressure points to rank 0
    double global_pres_point1 = 0.0;
    double global_pres_point2 = 0.0;

    // Assuming only one process has each pressure point, use MPI_MAX to gather the value
    MPI_Reduce(&pres_point1, &global_pres_point1, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&pres_point2, &global_pres_point2, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (this->mpi_rank == 0)
    {
        // Compute pressure difference
        double p_diff = global_pres_point1 - global_pres_point2;
        pcout << "Pressure difference (P(A) - P(B)) = " << p_diff << std::endl;

        /*// Write final aggregated results to CSV
        std::string output_path = this->get_output_directory() + "/lift_drag_output.csv";
        std::ofstream output_file(output_path, std::ios::app);
        if (output_file.is_open())
        {
            output_file << total_drag << ", " << total_lift << ", " << p_diff << "\n";
            output_file.close();
            std::cout << "Wrote aggregated drag/lift data to lift_drag_output.csv" << std::endl;
        }
        else
        {
            std::cerr << "Error: Unable to open lift_drag_output.csv for writing." << std::endl;
        }*/
    }
    // Ensure all processes have completed the reductions
    MPI_Barrier(MPI_COMM_WORLD);
}



