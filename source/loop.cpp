/*
loop.cpp
Loop object that will hold all information about the loop and be evolved in time.
*/

#include "loop.h"

Parameters Loop::parameters;
Terms Loop::terms;
HEATER Loop::heater;

Loop::Loop(char *config)
{
  tinyxml2::XMLElement *root;

  //Open file
  tinyxml2::XMLError load_ok = doc.LoadFile(config);
  if(load_ok != 0)
  {
    std::string filename(config);
    std::string error_message = "Failed to load XML configuration file " + filename;
    throw std::runtime_error(error_message);
  }
  //Parse file and read into data structure
  root = doc.FirstChildElement();
  //Numeric parameters
  parameters.total_time = std::stod(get_element_text(root,"total_time"));
  parameters.tau = std::stod(get_element_text(root,"tau"));
  parameters.tau_max = std::stod(get_element_text(root,"tau_max"));
  parameters.loop_length = std::stod(get_element_text(root,"loop_length"));
  parameters.area_ratio_tr_corona = std::stod(get_element_text(root, "area_ratio_tr_corona"));
  parameters.area_ratio_0_corona = std::stod(get_element_text(root, "area_ratio_0_corona"));
  parameters.adaptive_solver_error = std::stod(get_element_text(root,"adaptive_solver_error"));
  parameters.adaptive_solver_safety = std::stod(get_element_text(root,"adaptive_solver_safety"));
  parameters.saturation_limit = std::stod(get_element_text(root,"saturation_limit"));
  parameters.c1_cond0 = std::stod(get_element_text(root,"c1_cond0"));
  parameters.c1_rad0 = std::stod(get_element_text(root,"c1_rad0"));
  parameters.helium_to_hydrogen_ratio = std::stod(get_element_text(root,"helium_to_hydrogen_ratio"));
  parameters.surface_gravity = std::stod(get_element_text(root,"surface_gravity"));
  //Boolean parameters
  parameters.force_single_fluid = string2bool(get_element_text(root,"force_single_fluid"));
  parameters.use_c1_loss_correction = string2bool(get_element_text(root,"use_c1_loss_correction"));
  parameters.use_c1_grav_correction = string2bool(get_element_text(root,"use_c1_grav_correction"));
  parameters.use_flux_limiting = string2bool(get_element_text(root,"use_flux_limiting"));
  parameters.calculate_dem = string2bool(get_element_text(root,"calculate_dem"));
  parameters.use_adaptive_solver = string2bool(get_element_text(root,"use_adaptive_solver"));
  parameters.save_terms = string2bool(get_element_text(root,"save_terms"));
  //String parameters
  parameters.output_filename = get_element_text(root,"output_filename");

  //Estimate results array length
  parameters.N = int(std::ceil(parameters.total_time/parameters.tau));

  //Compute components of loop length
  double loop_length_ratio_tr_total = std::stod(get_element_text(root, "loop_length_ratio_tr_total"));
  parameters.loop_length_tr = parameters.loop_length * loop_length_ratio_tr_total;
  parameters.loop_length_corona = parameters.loop_length - parameters.loop_length_tr;

  //Initialize heating object
  heater = new Heater(get_element(root,"heating"));

  //Initialize DEM object
  if(parameters.calculate_dem)
  {
    parameters.dem_options = get_element(root,"dem");
  }

  // Call the setup function
  Setup();
}

Loop::Loop(void)
{
  heater = new Heater();
}

Loop::~Loop(void)
{
  //Destructor--free some stuff here
  doc.Clear();
  delete heater;
}

void Loop::Setup(void)
{
  // Calculate needed He abundance corrections
  CalculateAbundanceCorrection(parameters.helium_to_hydrogen_ratio);

  //Reserve memory for results
  results.time.resize(parameters.N);
  results.heat.resize(parameters.N);
  results.pressure_e.resize(parameters.N);
  results.pressure_i.resize(parameters.N);
  results.temperature_e.resize(parameters.N);
  results.temperature_i.resize(parameters.N);
  results.density.resize(parameters.N);
  results.velocity.resize(parameters.N);
}

state_type Loop::GetState(void)
{
  return __state;
}

void Loop::SetState(state_type state)
{
  __state = state;
}

state_type Loop::CalculateInitialConditions(void)
{
  // TODO: Modify initial conditions calculation to account for area expansion
  int i = 0;
  int i_max = 100;
  double tol = 1e-2;
  double temperature_old = (double)LARGEST_DOUBLE;
  double density_old = (double)LARGEST_DOUBLE;
  double temperature,density;
  double pe_initial,pi_initial;
  double radiative_loss;
  double error_temperature,error_density;
  double c1 = 2.0;
  double c2 = CalculateC2();
  double heat = heater->Get_Heating(0.0);
  state_type state;

  while(i<i_max)
  {
    if(i > 0)
    {
      c1 = CalculateC1(temperature_old, temperature_old, density_old);
    }
    temperature = c2*std::pow(3.5*c1/(1.0 + c1)*std::pow(parameters.loop_length_corona,2)*heat/(SPITZER_ELECTRON_CONDUCTIVITY + SPITZER_ION_CONDUCTIVITY),2.0/7.0);
    radiative_loss = CalculateRadiativeLoss(temperature);
    density = std::sqrt(heat/(radiative_loss*(1.0 + c1)));
    error_temperature = std::abs(temperature - temperature_old)/temperature;
    error_density = std::abs(density - density_old)/density;
    if(std::fmax(error_density,error_temperature) < tol)
    {
      break;
    }
    i++;
    temperature_old = temperature;
    density_old = density;
  }

  // Set current state in order pressure_e, pressure_i, density
  pi_initial = parameters.boltzmann_correction*BOLTZMANN_CONSTANT*density*temperature;
  pe_initial = BOLTZMANN_CONSTANT*density*temperature;
  if(parameters.force_single_fluid)
  {
    double p_initial = (pe_initial + pi_initial)/2.;
    pe_initial = p_initial;
    pi_initial = p_initial;
  }

  state = {{ pe_initial, pi_initial, density, temperature, temperature }};

  return state;
}

void Loop::PrintToFile(int num_steps)
{
  std::ofstream f;
  f.open(parameters.output_filename);
  for(int i=0;i<num_steps;i++)
  {
    //Use 10 decimal places when printing the time
    f << std::fixed << std::setprecision(std::numeric_limits<double>::digits10)
    << results.time[i] << "\t"
    << std::setprecision(6) << std::scientific
    << results.temperature_e[i] << "\t"
    << results.temperature_i[i] << "\t"
    << results.density[i] << "\t"
    << results.pressure_e[i] << "\t"
    << results.pressure_i[i] << "\t"
    << results.velocity[i] << "\t"
    << results.heat[i] << "\n";
  }
  f.close();

  if(parameters.save_terms)
  {
    f.open(parameters.output_filename+".terms");
    for(int i=0;i<num_steps;i++)
    {
      f << terms.f_e[i] << "\t"
      << terms.f_i[i] << "\t"
      << terms.c1[i] << "\t"
      << terms.radiative_loss[i] << "\n";
    }
    f.close();
  }
}

void Loop::CalculateDerivs(const state_type &state, state_type &derivs, double time)
{
  // TODO: Modify derivatives to account for area expansion
  double dpe_dt,dpi_dt,dn_dt,dTe_dt,dTi_dt;
  double psi_tr,psi_c,xi,R_tr,enthalpy_flux;

  double f_e = CalculateThermalConduction(state[3],state[2],"electron");
  double f_i = CalculateThermalConduction(state[4],state[2],"ion");
  double radiative_loss = CalculateRadiativeLoss(state[3]);
  double heat = heater->Get_Heating(time);
  double c1 = CalculateC1(state[3],state[4],state[2]);
  double c2 = CalculateC2();
  double c3 = CalculateC3();
  double collision_frequency = CalculateCollisionFrequency(state[3],state[2]);

  xi = state[0]/state[1];
  // R_tr = c1*std::pow(state[2],2)*radiative_loss*parameters.loop_length;
  R_c = std::pow(state[2],2)*radiative_loss*parameters.loop_length_corona;
  psi_tr = (f_e + R_tr - xi*f_i)/(1.0 + xi);
  psi_c = BOLTZMANN_CONSTANT*state[2]*collision_frequency*(state[4] - state[3]);
  enthalpy_flux = GAMMA_MINUS_ONE/GAMMA*(-f_e - R_tr + psi_tr);
  L_a = paramter.loop_length_corona + (A_tr/A_c) * parameter.loop_length_tr

  // dpe_dt = GAMMA_MINUS_ONE*(heat*heater->partition + 1.0/parameters.loop_length*(psi_tr - R_tr*(1.0 + 1.0/c1))) + psi_c;
  dpe_dt = GAMMA_MINUS_ONE*(heat*heater->partition + psi_c/L_a*(1 + parameters.area_ratio_tr_corona + (psi_tr/psi_c)) - R_c/L_a*(1 + parameters.area_ratio_tr_corona*c1));
  // dpi_dt = GAMMA_MINUS_ONE*(heat*(1.0 - heater->partition) - 1.0/parameters.loop_length*psi_tr) - psi_c;
  dpi_dt = GAMMA_MINUS_ONE*(heat*(1.0 - heater->partition) - psi_c/L_a*(1 + parameters.area_ratio_tr_corona + (psi_tr/psi_c)));
  // Divide pressure equally if single-fluid case
parameters.area_t
  if(parameters.force_single_fluid)
  {
    double tmp_dpe_dt = dpe_dt;
    dpe_dt = 0.5*(tmp_dpe_dt + dpi_dt);
    dpi_dt = 0.5*(tmp_dpe_dt + dpi_dt);
  }
 // dn_dt = c2/(c3*parameters.loop_length*BOLTZMANN_CONSTANT*state[3])*enthalpy_flux;
 dn_dt = ((c2 * (GAMMA_MINUS_ONE))/(c3*parameters.loop_length_corona*BOLTZMANN_CONSTANT*state[3]))*(parameters.area_ratio_tr_corona*(parameters.loop_length_corona/L_a)(c1-parameters.loop_length_tr/parameters.loop_length_corona)*R_c-(parameters.area_ratio_tr_corona*(paramters.loop_length_corona/L_a))*(psi_tr/psi_c - parameters.loop_length_tr/parameters.loop_length_corona)*psi_c+parameters.area_ratio_0_corona*f_e)

  dTe_dt = state[3]*(1/state[0]*dpe_dt - 1/state[2]*dn_dt);
  dTi_dt = state[4]*(1/state[1]*dpi_dt - 1/state[2]*dn_dt);

  derivs[0] = dpe_dt;
  derivs[1] = dpi_dt;
  derivs[2] = dn_dt;
  derivs[3] = dTe_dt;
  derivs[4] = dTi_dt;
}

void Loop::SaveResults(int i,double time)
{
  // Get heating profile and velocity
  double heat = heater->Get_Heating(time);
  double velocity = CalculateVelocity(__state[3], __state[4], __state[0]);

  // Save results to results structure
  if(i >= parameters.N)
  {
    results.time.push_back(time);
    results.heat.push_back(heat);
    results.temperature_e.push_back(__state[3]);
    results.temperature_i.push_back(__state[4]);
    results.pressure_e.push_back(__state[0]);
    results.pressure_i.push_back(__state[1]);
    results.density.push_back(__state[2]);
    results.velocity.push_back(velocity);
  }
  else
  {
    results.time[i] = time;
    results.heat[i] = heat;
    results.temperature_e[i] = __state[3];
    results.temperature_i[i] = __state[4];
    results.pressure_e[i] = __state[0];
    results.pressure_i[i] = __state[1];
    results.density[i] = __state[2];
    results.velocity[i] = velocity;
  }
}

void Loop::SaveTerms(void)
{
  // Calculate terms
  double f_e = CalculateThermalConduction(__state[3], __state[2], "electron");
  double f_i = CalculateThermalConduction(__state[4], __state[2], "ion");
  double c1 = CalculateC1(__state[3], __state[4], __state[2]);
  double radiative_loss = CalculateRadiativeLoss(__state[3]);

  // Save terms
  terms.f_e.push_back(f_e);
  terms.f_i.push_back(f_i);
  terms.c1.push_back(c1);
  terms.radiative_loss.push_back(radiative_loss);
}

double Loop::CalculateThermalConduction(double temperature, double density, std::string species)
{
  double kappa,mass,k_B;
  double f_c,f;
  double c2 = CalculateC2();

  if(species.compare("electron")==0)
  {
    kappa = SPITZER_ELECTRON_CONDUCTIVITY;
    mass = ELECTRON_MASS;
    k_B = BOLTZMANN_CONSTANT;
  }
  else
  {
    kappa = SPITZER_ION_CONDUCTIVITY;
    mass = parameters.ion_mass_correction*PROTON_MASS;
    k_B = parameters.boltzmann_correction*BOLTZMANN_CONSTANT;
  }

  f_c = -2.0/7.0*kappa*std::pow(temperature/c2,3.5)/parameters.loop_length_corona;

  if(parameters.use_flux_limiting)
  {
    double f_s = -parameters.saturation_limit*1.5/std::sqrt(mass)*density*std::pow(k_B*temperature,1.5);
    f = -f_c*f_s/std::sqrt(std::pow(f_c,2) + std::pow(f_s,2));
  }
  else
  {
    f = f_c;
  }

  return f;
}

double Loop::CalculateRadiativeLoss(double temperature)
{
  double chi, alpha;
  double log_temperature = std::log10(temperature);

  if( log_temperature <= 4.97 )
	{
	    chi = 1.09e-31;
	    alpha = 2.0;
	}
	else if( log_temperature <= 5.67 )
	{
	    chi = 8.87e-17;
	    alpha = -1.0;
	}
	else if( log_temperature <= 6.18 )
	{
	    chi = 1.90e-22;
	    alpha = 0.0;
	}
	else if( log_temperature <= 6.55 )
	{
	    chi = 3.53e-13;
	    alpha = -3.0/2.0;
	}
	else if( log_temperature <= 6.90 )
	{
	    chi = 3.46e-25;
	    alpha = 1.0/3.0;
	}
	else if( log_temperature <= 7.63 )
	{
	    chi = 5.49e-16;
	    alpha = -1.0;
	}
	else // NOTE: free-free radiation is included in the parameter values for log_10 T > 7.63
	{
	    chi = 1.96e-27;
	    alpha = 1.0/2.0;
	}

	return chi * std::pow( 10.0, (alpha*log_temperature) );
}

double Loop::CalculateCollisionFrequency(double temperature_e,double density)
{
  // TODO: find a reference for this formula
  double coulomb_logarithm = 23.0 - std::log(std::sqrt(density/1.0e+13)*std::pow(BOLTZMANN_CONSTANT*temperature_e/(1.602e-9),-1.5));
  return 16.0*SQRT_PI/3.0*ELECTRON_CHARGE_POWER_4/(parameters.ion_mass_correction*PROTON_MASS*ELECTRON_MASS)*std::pow(2.0*BOLTZMANN_CONSTANT*temperature_e/ELECTRON_MASS,-1.5)*density*coulomb_logarithm;
}

double Loop::CalculateC1(double temperature_e, double temperature_i, double density)
{
  // TODO: Modify c1 calculation for area expansion
  double c1;
  double density_eqm_2,density_ratio;

  double c1_eqm0 = 2.0;
  double c2 = CalculateC2();
  double grav_correction = 1.0;
  double loss_correction = 1.0;
  double scale_height = CalculateScaleHeight(temperature_e,temperature_i);
  double radiative_loss = CalculateRadiativeLoss(temperature_e);

  if(parameters.use_c1_grav_correction)
  {
    grav_correction = std::exp(4.0*std::sin(_PI_/5.0)*parameters.loop_length_corona/(_PI_*scale_height));
  }
  if(parameters.use_c1_loss_correction)
  {
    loss_correction = 1.95e-18/std::pow(temperature_e,2.0/3.0)/radiative_loss;
  }

  density_eqm_2 = (SPITZER_ELECTRON_CONDUCTIVITY + SPITZER_ION_CONDUCTIVITY)*std::pow(temperature_e/c2,3.5)/(3.5*std::pow(parameters.loop_length_corona,2)*c1_eqm0*loss_correction*grav_correction*radiative_loss);
  density_ratio = std::pow(density,2)/density_eqm_2;

  if(density_ratio<1.0)
  {
    c1 = (2.0*c1_eqm0 + parameters.c1_cond0*(1.0/density_ratio - 1.0))/(1.0 + 1.0/density_ratio);
  }
  else
  {
    c1 = (2.0*c1_eqm0 + parameters.c1_rad0*(density_ratio - 1.0))/(1.0 + density_ratio);
  }

  return c1*loss_correction*grav_correction;
}

double Loop::CalculateC2(void)
{
  return 0.9;
}

double Loop::CalculateC3(void)
{
  return 0.6;
}

double Loop::CalculateC4(void)
{
  return 1.0;
}

double Loop::CalculateScaleHeight(double temperature_e,double temperature_i)
{
  return BOLTZMANN_CONSTANT*(temperature_e + parameters.boltzmann_correction*temperature_i)/(parameters.ion_mass_correction*PROTON_MASS)/(parameters.surface_gravity * (double)SOLAR_SURFACE_GRAVITY);
}

void Loop::CalculateAbundanceCorrection(double helium_to_hydrogen_ratio)
{
  double z_avg = (1.0 + 2.0*helium_to_hydrogen_ratio)/(1.0 + helium_to_hydrogen_ratio);
  parameters.boltzmann_correction = 1.0/z_avg;
  parameters.ion_mass_correction = (1.0 + 4.0*helium_to_hydrogen_ratio)/(2.0 + 3.0*helium_to_hydrogen_ratio)*(1.0 + z_avg)/z_avg;
}

double Loop::CalculateVelocity(double temperature_e, double temperature_i, double pressure_e)
{
  //TODO: Modify enthalpy calculation for expansion
  double c4 = CalculateC4();
  double density = pressure_e/(BOLTZMANN_CONSTANT*temperature_e);
  double c1 = CalculateC1(temperature_e,temperature_i,density);
  double R_tr = c1*std::pow(density,2)*CalculateRadiativeLoss(temperature_e)*parameters.loop_length_corona;
  double fe = CalculateThermalConduction(temperature_e,density,"electron");
  double fi = CalculateThermalConduction(temperature_i,density,"ion");
  double sc = CalculateScaleHeight(temperature_e,temperature_i);
  double xi = temperature_e/temperature_i/parameters.boltzmann_correction;

  double coefficient = c4*xi*GAMMA_MINUS_ONE/(GAMMA*(xi+1));
  double pressure_e_0 = pressure_e*std::exp(2.0*parameters.loop_length_corona*std::sin(_PI_/5.0)/(_PI_*sc));
  double enthalpy_flux = -(fe + fi + R_tr);

  return coefficient*enthalpy_flux/pressure_e_0;
}

double Loop::ControlTimeStep(const state_type &state, double time, double tau)
{
  double tau_tc = 4e-10*state[2]*pow(parameters.loop_length_corona, 2)*pow(std::fmax(state[3], state[4]), -2.5);
  // Limit abrupt changes in the timestep with safety factor
  tau = std::fmax(std::fmin(tau, 0.5*tau_tc), parameters.adaptive_solver_safety*tau);
  // Control maximum timestep
  tau = std::fmin(tau, parameters.tau_max);
  tau = std::fmin(tau, heater->Get_Time_To_Next_Heating_Change(time));
  return tau;
}
