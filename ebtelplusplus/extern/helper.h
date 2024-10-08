/*
helper.h
General purpose includes to be used everywhere
*/

#ifndef HELPER_H
#define HELPER_H

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <limits>
#include <filesystem>
#include <set>
#include "boost/array.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace fs = std::filesystem;
namespace py = pybind11;


// Structure to hold all input parameters
struct Parameters {
  /* Total simulation time (in s) */
  double total_time;
  /* Timestep (in s); when using adaptive solver, the initial timestep*/
  double tau;
  /* Maximum allowed timestep (in s) when using adaptive solver */
  double tau_max;
  /* Coronal portion of loop half length (in cm) */
  double loop_length_corona;
  /* Ratio of transition region portion of loop to coronal portion of loop */
  double loop_length_ratio_tr_corona;
  /* Ratio of average cross-sectional area in transition region to average cross-sectional area in corona */
  double area_ratio_tr_corona;
  /* Ratio of cross-sectional area at TR-corona interface to average cross-sectional area in corona */
  double area_ratio_0_corona;
  /* Truncation error tolerance for adaptive solver */
  double adaptive_solver_error;
  /* Safety factor on allowed timestep for adaptive solver */
  double adaptive_solver_safety;
  /* Heat flux saturation limit; 1/6 is a typical value */
  double saturation_limit;
  /* Nominal conductive c1 value; 6.0 is recommended*/
  double c1_conduction;
  /* Nominal radiative c1 value; 0.6 is recommended */
  double c1_radiation;
  /* Force collision frequency to 0.9 s^-1 to emulate single-species fluid */
  bool force_single_fluid;
  /* Switch for radiative loss correction to c1 */
  bool use_c1_radiation_correction;
  /* Switch for gravitational correction to c1 */
  bool use_c1_gravity_correction;
  /* Switch for classical Spitzer conductivity */
  bool use_flux_limiting;
  /* Switch for calculating DEM; if True, runtimes will be much longer */
  bool calculate_dem;
  /* Switch for using the adaptive solver option */
  bool use_adaptive_solver;
  /* What radiative losses to assume:
        power_law: use the default power-law fit to radiative losses (Klimchuk et al 2008)
        variable: use a look-up table with time-variable abundance factor f
        photospheric: use a look-up table with abundance factor f = 1
        coronal: use a look-up table with abundance factor f = 4
  */
  std::string radiative_loss;
  std::string radiation_data_dir;
  /* Switch for using look-up tables for the radiative loss calculation */
  bool use_lookup_table_losses;
  /* DEM calculation parameters */
  bool dem_use_new_tr_method;
  int dem_temperature_bins;
  double dem_temperature_min;
  double dem_temperature_max;
  /* Correction to ion equation of state */
  double boltzmann_correction;
  /* Ion mass correction to account for He abundance */
  double ion_mass_correction;
  /* Ratio of helium to hydrogen  */
  double helium_to_hydrogen_ratio;
  /* Gravitational acceleration at stellar surface */
  double surface_gravity;
  /* Number of grid points */
  size_t N;
  
  /* Variables and arrays used for variable abundance radiative losses */
  /* The temperature values in the look-up table for radiative losses */
  //double log10_temperature_array[101];
  std::vector<double> log10_temperature_array;
  std::vector<double> log10_density_array;
  std::vector<double> abundance_array;
  /* The look-up table's radiative loss rate as a 
   * function of [abundance][temperature][density] */
  std::vector<std::vector<std::vector<double> > > log10_loss_rate_array;
  /* The density in the corona before upflows begin, used to calculate
   * the change in abundance factor */
  double initial_density;
  /* The density at the previous time step*/
  double previous_density;
  /* The abundance factor in the corona before upflows begin */
  double initial_abundance_factor;
  /* The abundance factor at the previous time step */
  double previous_abundance_factor;
  /* Whether the flows are upflowing (into the corona) or not, 
   * which is used to determine whether the abundance factor changes
   * with time. That is, flows out of the corona do not affect 
   * the abundance factor. */
  bool upflowing;
  /* The power law losses are used to calculate the initial conditions, so 
   * this bool is used to tell the code that. */
  bool initial_radiation;
};

// Structure to hold all results
struct Results {
  /* Time (in s) */
  std::vector<double> time;
  /* Electron temperature (in K) */
  std::vector<double> temperature_e;
  /* Ion temperature (in K) */
  std::vector<double> temperature_i;
  /* Electron pressure (in dyne cm^-2 s^-1) */
  std::vector<double> pressure_e;
  /* Ion pressure (in dyne cm^-2 s^-1) */
  std::vector<double> pressure_i;
  /* Number density (in cm^-3) */
  std::vector<double> density;
  /* Velocity (in cm s^-1) */
  std::vector<double> velocity;
  /* Heating rate (in erg cm^-3 s^-1) */
  std::vector<double> heat;
};

// Structure to hold equation terms
struct Terms {
  /* Electron heat flux (in erg cm^-2 s^-1)*/
  std::vector<double> f_e;
  /* Ion heat flux (in erg cm^-2 s^-1)*/
  std::vector<double> f_i;
  /* Radiative_loss (in erg cm^3 s^-1) */
  std::vector<double> radiative_loss;
  /* c1 coefficient */
  std::vector<double> c1;
};

// Generic type for state vectors and derivatives
typedef boost::array<double, 5> state_type;

#endif
