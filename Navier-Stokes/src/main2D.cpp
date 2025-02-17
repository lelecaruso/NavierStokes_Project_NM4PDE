#include "NavierStokes2D.hpp"

// Main function.
int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //const std::string mesh_file_name = argc > 1 ? argv[1] : "../mesh/cilinder_2D_fine.msh";
  const std::string mesh_file_name = argc > 1 ? argv[1] : "../mesh/cilinder_2D_coarse.msh";
  //const std::string mesh_file_name = argc > 1 ? argv[1] : "../mesh/mesh-square-5.msh"; // per facilita di risoluzione del sistema

  //TAYLOR-HOOD 
  const unsigned int degree_velocity = 2;
  const unsigned int degree_pressure = 1;

  //const double T = 24;

  /*
  test 2
  const double Re = 100.0;
  const double A = 0.2175;
  const double B = -5.106;

  const double freq = (A + (B/Re));

  */

  //const double T_test2 = 1.0/freq;

  const double T = 0.0;  //test 3
  const double deltat = 0.05;

  dealii::Timer timer;
  // Start the timer
  timer.restart();

   NavierStokes problem(mesh_file_name, degree_velocity, degree_pressure, T, deltat); //test3
  //NavierStokes problem(mesh_file_name, degree_velocity, degree_pressure, T_test2, deltat);


  if(T>0)
  {
    problem.setup();
    //UNSTEADY CASE
    problem.solve();

  }

  else
  {
    //STEADY CASE
  
  //NavierStokes problem(mesh_file_name, degree_velocity, degree_pressure, T_test2, deltat);
  problem.setupNewton();
  problem.solveNewton();
    
  }

  // Stop the timer
  timer.stop();

  // Output the elapsed time
  if(rank == 0)
    std::cout << "Time taken to solve ENTIRE Navier Stokes problem: " << timer.wall_time() << " seconds" << std::endl;

  if (rank == 0)
  {
    const std::string output_filename = "forces_results_2D_2case.csv";
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
