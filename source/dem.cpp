/*
dem.cpp
Methods for DEM class
*/

#include "dem.h"


Dem::Dem(void)
{
  // Default constructor
}

Dem::Dem(LOOP loop_object)
{
  loop = loop_object;
  // Set some parameters for later calculations
  use_new_method = string2bool(get_element_text(loop->parameters.dem_options,"use_new_method"));
  // Configure temperature vector from inputs
  tinyxml2::XMLElement * temperature_node = get_element(loop->parameters.dem_options,"temperature");
  int nbins = std::stoi(temperature_node->Attribute("bins"));
  double temperature_min = pow(10.0,std::stod(temperature_node->Attribute("log_min")));
  double temperature_max = pow(10.0,std::stod(temperature_node->Attribute("log_max")));
  double delta_temperature = (log10(temperature_max) - log10(temperature_min))/(nbins-1);
  // Set temperature bins and associated radiative losses
  __temperature.resize(nbins);
  __radiative_loss.resize(nbins);
  for(int i=0;i<nbins;i++)
  {
    __temperature[i] = temperature_min*pow(10.0,i*delta_temperature);
    __radiative_loss[i] = loop->CalculateRadiativeLoss(__temperature[i]);
  }
  // Resize DEM arrays
  dem_TR.resize(loop->parameters.N);
  dem_corona.resize(loop->parameters.N);
}

Dem::~Dem(void)
{
  // Destructor--free some stuff here if needed
}

void Dem::CalculateDEM(int i)
{
  // TODO: check whether DEM calculation needs to be modified for expansion
  state_type loop_state = loop->GetState();
  double velocity = loop->CalculateVelocity(loop_state[3],loop_state[4],loop_state[0]);
  double scale_height = loop->CalculateScaleHeight(loop_state[3],loop_state[4]);
  double f_e = loop->CalculateThermalConduction(loop_state[3],loop_state[2],"electron");
  double R_tr = loop->CalculateC1(loop_state[3],loop_state[4],loop_state[2])*pow(loop_state[2],2)*loop->CalculateRadiativeLoss(loop_state[3])*loop->parameters.loop_length_corona;
  // Calculate coronal temperature range
  double temperature_corona_max = fmax(loop_state[3]/loop->CalculateC2(),1.1e+4);
  double temperature_corona_min = fmax(loop_state[3]*(2.0 - 1.0/loop->CalculateC2()),1.0e+4);
  // Calculate coronal emission
  double delta_temperature = pow(10.0,0.5/100.0)*temperature_corona_max - pow(10.0,-0.5/100.0)*temperature_corona_min;
  double coronal_emission = 2.0*pow(loop_state[2],2)*loop->parameters.loop_length_corona/delta_temperature;

  bool dem_tr_negative = false;
  std::vector<double> tmp_dem_corona(__temperature.size()), tmp_dem_tr(__temperature.size());

  for(int j=0;j<__temperature.size();j++)
  {
    // Coronal DEM
    tmp_dem_corona[j] = 0.0;
    if(__temperature[j]<=temperature_corona_max && __temperature[j]>=temperature_corona_min)
    {
      tmp_dem_corona[j] = coronal_emission;
    }
    // Transition Region DEM
    tmp_dem_tr[j] = 0.0;
    if(__temperature[j]<loop->CalculateC3()/loop->CalculateC2()*loop_state[3])
    {
      if(dem_tr_negative && i>0)
      {
        tmp_dem_tr[j] = dem_TR[i-1][j];
      }
      else
      {
        tmp_dem_tr[j] = CalculateDEMTR(j,loop_state[2],velocity,loop_state[0],scale_height,R_tr,f_e);
        if(tmp_dem_tr[j] < 0.0)
        {
          dem_tr_negative=true;
          std::cout << "Negative DEM at timestep " << i << std::endl;
          j = -1;
          continue;
        }
      }
    }
  }
  if(i>=loop->parameters.N)
  {
    dem_TR.push_back(tmp_dem_tr);
    dem_corona.push_back(tmp_dem_corona);
  }
  else
  {
    dem_TR[i] = tmp_dem_tr;
    dem_corona[i] = tmp_dem_corona;
  }
}

void Dem::PrintToFile(int num_steps)
{
  // Open file streams
  std::ofstream f_corona;
  std::ofstream f_tr;
  f_corona.open(loop->parameters.output_filename+".dem_corona");
  f_tr.open(loop->parameters.output_filename+".dem_tr");
  // First row of each file is the temperature array
  for(int j=0;j<__temperature.size();j++)
  {
    f_corona << __temperature[j] << "\t";
    f_tr << __temperature[j] << "\t";
  }
  f_corona << "\n";
  f_tr << "\n";
  // Print TR and corona DEM at each timestep
  for(int i=0;i<num_steps;i++)
  {
    for(int j=0;j<dem_TR[i].size();j++)
    {
      f_corona << dem_corona[i][j] << "\t";
      f_tr << dem_TR[i][j] << "\t";
    }
    f_corona << "\n";
    f_tr << "\n";
  }

  // Close the file
  f_corona.close();
  f_tr.close();
}

double Dem::CalculateDEMTR(int j,double density,double velocity,double pressure,double scale_height,double R_tr,double f_e)
{
  double dem_tr;

  if(use_new_method)
  {
    double a = (SPITZER_ELECTRON_CONDUCTIVITY + SPITZER_ION_CONDUCTIVITY)*pow(__temperature[j],1.5);
    double b = -GAMMA*(1.0+ loop->parameters.boltzmann_correction)*BOLTZMANN_CONSTANT/GAMMA_MINUS_ONE*density*velocity;
    double density_squared = pow(pressure/BOLTZMANN_CONSTANT/__temperature[j],2)*exp(4.0*loop->parameters.loop_length_corona*sin(_PI_/5.0)/scale_height/_PI_);
    double c = -density_squared*__radiative_loss[j];
    double dTds_plus = (-b + sqrt(pow(b,2) - 4.0*a*c))/(2.0*a);
    double dTds_minus = (-b - sqrt(pow(b,2) - 4.0*a*c))/(2.0*a);
    double dTds = fmax(dTds_plus,dTds_minus);
    dem_tr = 2.0*density_squared/dTds;
  }
  else
  {
    double dem_evap = GAMMA_MINUS_ONE*(SPITZER_ELECTRON_CONDUCTIVITY+SPITZER_ION_CONDUCTIVITY)/GAMMA/(1.0+loop->parameters.boltzmann_correction)/pow(BOLTZMANN_CONSTANT,3)*pow(pressure,2)/(density*velocity*sqrt(__temperature[j]));
    double dem_condense = -GAMMA*(1.0 + loop->parameters.boltzmann_correction)*BOLTZMANN_CONSTANT*density*velocity/GAMMA_MINUS_ONE/__radiative_loss[j];
    double dem_eqm = sqrt(2.0*(SPITZER_ELECTRON_CONDUCTIVITY+SPITZER_ION_CONDUCTIVITY)/7.0/__radiative_loss[j])*pressure/BOLTZMANN_CONSTANT*pow(__temperature[j],-0.25);

    double f_e_plus_R_tr = f_e + R_tr;
    if(f_e_plus_R_tr==0.0)
    {
      f_e_plus_R_tr = 1.0e+10*f_e;
    }

    dem_tr = 2.0*(f_e*dem_evap - f_e*R_tr/f_e_plus_R_tr*dem_eqm + R_tr*dem_condense)/(f_e - f_e*R_tr/f_e_plus_R_tr + R_tr);
  }

  return dem_tr;
}
