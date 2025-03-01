#include "Convergence3D.hpp"
#include <deal.II/base/convergence_table.h>

// Main function.
int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   ConvergenceTable table;

  const std::vector<std::string> meshes = {"../mesh/mesh-cube-5.msh",
                                          "../mesh/mesh-cube-10.msh"
                                          };
  const std::vector<double>      h_vals = {1.0 / 5.0,
                                           1.0 / 10.0
                                           };
  
  
  //const std::string mesh_file_name = argc > 1 ? argv[1] : "../mesh/mesh-cube-5.msh";

  //TAYLOR-HOOD 
  const unsigned int degree_velocity = 2;
  const unsigned int degree_pressure = 1;
  //std::vector<double> errors_L2;

 
  const double T = 0.005; 
  const double deltat = 0.001;

  dealii::Timer timer;
  // Start the timer
  timer.restart();

  std::ofstream convergence_file("convergence.csv");
  convergence_file << "h,eL2" << std::endl;

  for (unsigned int i = 0; i < meshes.size(); ++i){

  NavierStokes problem(meshes[i], degree_velocity, degree_pressure, T, deltat); //test3

  problem.setup();
  problem.solve();
 
  const double error_L2 = problem.compute_error(VectorTools::L2_norm);
      
  table.add_value("h", h_vals[i]);
  table.add_value("L2", error_L2);
      

  convergence_file << h_vals[i] << "," << error_L2 << ","  << std::endl;
  // Stop the timer
  timer.stop();

}

  if(rank == 0)
  {
  std::cout << "Time taken to solve ENTIRE Navier Stokes problem: " << timer.wall_time() << " seconds" << std::endl;
  table.evaluate_all_convergence_rates(ConvergenceTable::reduction_rate_log2);
  table.set_scientific("L2", true);
  table.write_text(std::cout);
  }

  return 0;
}
