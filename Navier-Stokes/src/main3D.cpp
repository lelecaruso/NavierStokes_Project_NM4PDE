#include "../include/NavierStokes3D.hpp"

// Main function.
int main(int argc, char *argv[])
{ 
   Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
      // Default test case
  int test_case = 2;
 /* if( rank == 0)
  {
    std::cout << "Choose a test case for 2D Navier-Stokes:\n"
          "1 - TEST 1 ( non implementato )\n"
          "2 - TEST 2\n"
          "3 - TEST 3 ( non implementato )\n---> ";
      
    std::cin >> test_case;  // Read user input into test_case
    if (test_case <= 1 || test_case > 3) 
          {
              std::cerr << "Invalid test case number. Using default (2)." << std::endl;
              test_case = 2;
          }    
  }*/

  MPI_Bcast(&test_case, 1, MPI_INT, 0, MPI_COMM_WORLD);


  const std::string mesh_file_name = argc > 1 ? argv[1] : "../mesh/Parallelepiped3D.msh";

  // Taylor-Hood elements
  const unsigned int degree_velocity = 2;
  const unsigned int degree_pressure = 1;

  const double T = 4;  
  const double deltat = 0.0002; 

  dealii::Timer timer;
  // Start the timer for solving the entire problem
  timer.restart();

  NavierStokes problem(mesh_file_name, degree_velocity, degree_pressure, T, deltat, test_case); 

  problem.setup();
  problem.solve();

  // Stop the timer
  timer.stop();

  // Output the elapsed time
  if(rank == 0)
    std::cout << "Time taken to solve ENTIRE Navier Stokes problem: " << timer.wall_time() << " seconds" << std::endl;

  if (rank == 0)
  {
    const std::string output_filename = "forces_results_3D_2case.csv";
    std::ofstream outputFile(output_filename);

    if (!outputFile.is_open())
    {
      std::cerr << "Error opening output file" << std::endl;
      return -1;
    }
    outputFile << "Iteration, Drag, Lift, Coeff Drag, CoeffLift, time prec, time solve" << std::endl;

    for (size_t ite = 0; ite < problem.vec_drag.size(); ite++)
    {
      outputFile << ite * deltat << ", " << problem.vec_drag[ite] << ", " << problem.vec_lift_coeff[ite] << ", " 
                << problem.vec_drag_coeff[ite] << ", " << problem.vec_lift_coeff[ite] << ", "
                << problem.time_prec[ite] << ", " << problem.time_solve[ite]
                << std::endl;
    }
    outputFile.close();
  }

  return 0;
}
