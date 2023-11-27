/*
loop.h
Loop class definition
*/

#ifndef LOOP_H
#define LOOP_H

#include "helper.h"
#include "heater.h"
#include "util/constants.h"

// Loop object
//
// Class for holding all of the information about the loop. It can
// be passed a configuration and after the initial conditions are
// set up, it is evolved through time.
//
class Loop {
private:

  /* Results structure */
  Results results;

  /* Pointer to doc tree */
  tinyxml2::XMLDocument doc;

  /* Current state of the system */
  state_type __state;

  // Calculate c4
  // @return ratio of average to base velocity
  //
  // Calculate the ratio of average to base velocity. Set to 1 for now
  //
  double CalculateC4(void);

  // Calculate coulomb collision frequency
  // @temperature_e electron temperature (in K)
  // @density number density (in cm${^-3}$)
  //
  // Calculate the coulomb collision frequency for binary collisions
  // between electrons and ions according to Eq. 2.5e and Section 3 of
  // [Braginskii (1965)](http://adsabs.harvard.edu/abs/1965RvPP....1..205B).
  //
  // @return coulomb collision frequency (in s^-1)
  //
  static double CalculateCollisionFrequency(double temperature_e,double density);

  // Calculate correction for He abundance
  //
  void CalculateAbundanceCorrection(double helium_to_hydrogen_ratio);

public:

  /* Instance of the <Heater> object */
  static HEATER heater;

  /* Parameter structure*/
  static Parameters parameters;

  /* Terms structure */
  static Terms terms;

  // Constructor
  // @config main configuration file
  //
  // Setup the loop object by reading in parameters from the configuration
  // file <ebtel_config> into the <parameters> structure. The constructor also creates the <heater> object for calculating
  // the heating profile.
  //
  Loop(char * config);

  // Default constructor
  //
  // Create object without any configuration. Useful if
  // parameters are going to be read in from memory rather
  // than from a configuration file.
  //
  Loop(void);

  /* Destructor */
  ~Loop(void);

  // Setup object
  //
  // Allocate space for results and set some parameters. If you
  // create the object with the empty constructor, you need to
  // call this later on. If you use the default config file approach,
  // this is called automatically.
  //
  void Setup(void);

  // Set initial conditions
  //
  // Calculate the equilibrium values of pressure, temperature, and
  // density based on the supplied loop half-length and initial heating
  // according to the equilibrium solutions of the EBTEL equations.
  //
  // @return the initial state of the loop
  //
  state_type CalculateInitialConditions(void);

  // Print results to file
  // @num_steps number of steps taken by the integration routine
  //
  // Print results of EBTEL simulation to filename supplied in configuration file.
  // See documentation for the structure of the file itself and
  // instructions on how to parse it.
  //
  void PrintToFile(int num_steps);

  // Save results to structure
  // @i Current timestep
  // @time Current time (in s)
  //
  void SaveResults(int i, double time);

  // Save equation terms to structure
  //
  void SaveTerms(void);

  // Return current state publicly
  //
  // @return vector holding the current state of the loop
  //
  state_type GetState(void);

  // Set current state
  // @state electron pressure, ion pressure, and density to set as the current loop state
  //
  void SetState(state_type state);

  // Calculate $c_1$
  // @temperature_e electron temperature (in K)
  // @temperature_i ion temperature (in K)
  // @density number density (in cm$^{-3}$)
  //
  // Calculate the $c_1$ parameter, the ratio between the transition and coronal radiative losses
  //
  // @return $c_1$ parameter
  //
  static double CalculateC1(double temperature_e,double temperature_i,double density);

  // Calculate $c_2$
  //
  // Calculate the ratio of the average to apex temperature. Fixed at 0.9 for now.
  //
  // @return $c_2$ parameter
  //
  static double CalculateC2(void);

  // Calculate $c_3$
  //
  // Calculate the ratio of the base (corona/interface point)to apex temperature.
  // Fixed at 0.6 for now.
  //
  // @return $c_3$ parameter
  //
  static double CalculateC3(void);

  // Calculate velocity
  //
  // Calculate the velocity using the base total pressure and the enthalpy
  // flux as determined by our EBTEL equations.
  //
  // @return velocity averaged over the loop half-length (in cm s$^{-1}$)
  //
  double CalculateVelocity(void);

  // Calculate temperature scale height
  // @temperature_e electron temperature (in K)
  // @temperature_i ion temperature (in K)
  //
  // Calculate the temperature scale height of the loop. This parameter is used when
  // accounting for gravitational stratification in the model.
  //
  // @return the temperature scale height (in cm)
  //
  static double CalculateScaleHeight(double temperature_e,double temperature_i);

  // Calculate thermal conduction
  // @temperature temperature (in K)
  // @density density (in cm$^{-3}$)
  // @species either "electron" or "ion"
  //
  // Calculate the heat flux for either the electrons or ions, depending on the value of <species>.
  // The classical Spitzer formula is used. If <Parameters.use_flux_limiting> is set to true
  // in the configuration file, then a flux limiter is used to prevent runaway cooling.
  //
  // @return electron or ion heat flux (in erg cm$^{-2}$ s$^{-1}$)
  //
  static double CalculateThermalConduction(double temperature,double density,std::string species);

  // Calculate radiative losses
  // @temperature electron temperature (in K)
  //
  // Calculate the radiative loss at a particular temperature using the power-law approximation.
  // The formulation used here is based on the calculations of John Raymond (1994, private
  // communication) and twice the coronal abundances of Meyer (1985). This is the same power-law
  // radiative loss function as is implemented in the HYDRAD code and the EBTEL IDL code.
  //
  // @return radiative loss function (in erg cm$^3$ s$^{-1}$)
  //
  static double CalculateRadiativeLoss(double temperature);

  // Calculate derivatives of EBTEL equations
  // @state current state of the loop
  // @time current time (in s)
  //
  // Calculate the rate of change of the electron pressure, ion pressure, and
  // density according to the two-fluid EBTEL equations. A full derivation of these
  // equations can be found in Appendix B of [Barnes et al. (2016)](http://adsabs.harvard.edu/abs/2016ApJ...829...31B).
  //
  // @return the time derivatives of the electron pressure, ion pressure, and density
  //
  static void CalculateDerivatives(const state_type &state, state_type &derivs, double time);

  // Control the integration time step
  // @tau time step (in s)
  // @state current state of the loop
  // @time current time (in s)
  //
  // Calculate the time step to account for thermal conduction as well as the time
  // until the next heating event such that a heating event does not get skipped.
  //
  // @return updated time step (in s)
  static double ControlTimeStep(const state_type &state, double time, double tau);
};
// Pointer to the <Loop> class
typedef Loop* LOOP;

#endif
