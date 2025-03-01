#include "Convergence3D.hpp"

// Main function.
int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  const std::string mesh_file_name = argc > 1 ? argv[1] : "../mesh/mesh-cube-5.msh";

  //TAYLOR-HOOD 
  const unsigned int degree_velocity = 2;
  const unsigned int degree_pressure = 1;
  //std::vector<double> errors_L2;

 
<<<<<<< HEAD
  const double T = 1; 
  const double deltat = 0.1;
=======
  const double T = 0.1; 
  const double deltat = 0.001;
>>>>>>> 9da681d37909a019a4f3a0743cbc9f1cc5465b3b

  dealii::Timer timer;
  // Start the timer
  timer.restart();

  NavierStokes problem(mesh_file_name, degree_velocity, degree_pressure, T, deltat); //test3

  problem.setup();
  problem.solve();
  double errors_L2 = problem.compute_error(VectorTools::L2_norm);
  //errors_L2.push_back();
  
  // Stop the timer
  timer.stop();

  // Output the elapsed time
  if(rank == 0)
{
    std::cout << "Time taken to solve ENTIRE Navier Stokes problem: " << timer.wall_time() << " seconds" << std::endl;
    std::cout << std::scientific << "L2-norm of the error at final time step = " << errors_L2;
}


  return 0;
}
