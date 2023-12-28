/*
ebtel++
A code for computing the evolution of dynamically heated, spatially-averaged solar coronal loops.
*/

#include "boost/program_options.hpp"
#include "ebtel.h"

int main(int argc, char *argv[])
{
  //Declarations
  char config[256];

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

  //Run ebtel
  run(config);

  return 0;
}
