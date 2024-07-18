/*
ebtel++
A code for computing the evolution of dynamically heated, spatially-averaged solar coronal loops.
*/

#include <time.h>
#include "boost/program_options.hpp"
#include "boost/numeric/odeint.hpp"
#include "loop.h"
#include "dem.h"
#include "observer.h"

int main(int argc, char *argv[])
{
  //Declarations
  int num_steps;
  state_type state;
  char config[256];
  LOOP loop;
  DEM dem;
  OBSERVER obs;

  //Parse command line options with boost
  namespace po = boost::program_options;
  po::options_description description("\nebtel++\nA code for efficiently computing the evolution of dynamically-heated coronal loops. Based on the Enthalpy-based Thermal Evolution of Loops (EBTEL) model of Klimchuk et al. (2008) and Cargill et al. (2012). For more information, consult the documentation.\n\nOptional command line arguments");
  description.add_options()
    ("help,h","This help message")
    ("config,c",po::value<std::string>()->default_value("config/ebtel.example.cfg.xml"),"Configuration file for EBTEL.");
  po::variables_map vm;
  po::store(po::command_line_parser(argc,argv).options(description).run(), vm);
  if(vm.count("help"))
  {
  	std::cout << description;
  	return 0;
  }
  po::notify(vm);

  //Copy parameter to char array
  std::strcpy(config,vm["config"].as<std::string>().c_str());

  // Create loop object
  loop = new Loop(config);
  // Create DEM object
  if(loop->parameters.calculate_dem)
  {
    dem = new Dem(loop);
  }
  else
  {
    dem = new Dem();
  }
  // Configure observer
  obs = new Observer(loop,dem);

  // Set initional conditions of the loop
  state = loop->CalculateInitialConditions();
  // Set initial state for loop and dem
  obs->Observe(state, 0.0);

  // Set up Runge-Kutta integrator
  typedef boost::numeric::odeint::runge_kutta_cash_karp54< state_type > stepper_type;
  auto controlled_stepper = boost::numeric::odeint::make_controlled(loop->parameters.adaptive_solver_error, loop->parameters.adaptive_solver_error, stepper_type());

  // Integrate
  num_steps = 0;
  if(loop->parameters.use_adaptive_solver)
  {
    // Set maximum number of allowed failures
    int max_failures = 1000;
    // Initialize time and timestep
    double tau = loop->parameters.tau;
    double t = loop->parameters.tau;
    double old_tau,old_t;
    // Start integration loop
    while(t<loop->parameters.total_time)
    {
      int fail = 1;
      int num_failures = 0;
      while(fail>0)
      {
        // Throw error if exceeded max number of failures to avoid infinite loop
        if(num_failures>max_failures)
        {
          throw std::runtime_error("Adaptive solver exceeded maximum number of allowed failures.");
        }
        old_tau = tau;
        old_t = t;
        fail = controlled_stepper.try_step(loop->CalculateDerivs,state,t,tau);
        // Force NaNs to fail
        if(!fail) fail = obs->CheckNan(state,t,tau,old_t,old_tau);
        num_failures++;
      }
      // Enforce thermal conduction timescale limit
      double tau_tc = 4e-10*state[2]*pow(loop->parameters.loop_length,2)*pow(std::fmax(state[3],state[4]),-2.5);
      // Limit abrupt changes in the timestep with safety factor
      tau = std::fmax(std::fmin(tau,0.5*tau_tc),loop->parameters.adaptive_solver_safety*tau);
      // Control maximum timestep
      tau = std::fmin(tau,loop->parameters.tau_max);
      tau = std::fmin(tau,loop->CalculateTimeNextHeating(t));
      // Save the state
      obs->Observe(state,t);
      num_steps += 1;
    }
  }
  else
  {
    // Constant timestep integration
    num_steps = boost::numeric::odeint::integrate_const( controlled_stepper, loop->CalculateDerivs, state, loop->parameters.tau, loop->parameters.total_time, loop->parameters.tau, obs->Observe);
    if(obs->CheckNan(state))
    {
        throw std::runtime_error("NaNs were detected in the output.  Check the input configuration.");
    }
  }

  //Print results to file
  if(!loop->parameters.use_adaptive_solver)
  {
    num_steps = std::fmin(loop->parameters.N,num_steps);
  }
  loop->PrintToFile(num_steps);
  if(loop->parameters.calculate_dem)
  {
    dem->PrintToFile(num_steps);
  }

  //Cleanup
  delete obs;
  delete loop;
  delete dem;

  return 0;
}
