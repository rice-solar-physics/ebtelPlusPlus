/*
dem.h
Class definition for DEM object
*/

#ifndef DEM_H
#define DEM_H

#include "helper.h"
#include "loop.h"
#include "../Radiation_Model/source/radiation.h"
#include "../rsp_toolkit/source/xmlreader.h"
#include "../rsp_toolkit/source/file.h"
#include "../rsp_toolkit/source/constants.h"

// DEM object
//
// Class for holding all of the methods needed to calculate
// the differential emission measure in the transition region
// and the corona. Requires the <loop> object for knowledge about
// the evolution of the coronal loop.
//
class Dem{
private:
  /* Loop object */
  LOOP loop;

  /* Method option for DEM TR calculation */
  bool use_new_method;

  /* Temperature range */
  std::vector<double> __temperature;

  /* Radiative loss */
  std::vector<double> __radiative_loss;

  /*Transition region DEM*/
  std::vector<std::vector<double> > dem_TR;

  /*Coronal DEM*/
  std::vector<std::vector<double> > dem_corona;

  // Calculate TR DEM
  //
  double CalculateDEMTR(int j,double density,double velocity,double pressure,double scale_height,double R_tr,double f_e);

public:
  // Default constructor
  // @loop <Loop> object that provides needed parameters and methods
  //
  // Setup Dem object to calculate differential emission measure in both the
  // transition region and the corona.
  //
  Dem(LOOP loop);

  // Destructor
  //
  ~Dem(void);

  // Calculate DEM
  // @i Timestep index
  //
  // Front end for DEM calculations. Calls methods to calculate both
  // the transition region and coronal DEM across the entire specified
  // temperature range.
  //
  void CalculateDEM(int i);

  // Print results to file
  // @excess number of timesteps to clip from the end of the array; only nonzero for adaptive solver
  //
  // Print coronal and transition region DEM arrays to separate files.
  // The filenames are the output filename as given in <loop>,
  // suffixed by `.dem_corona` and `.dem_tr`, respectively. The first
  // row of each file is the temperature vector, <__temperature>.
  //
  void PrintToFile(int excess);
};
// Pointer to the <Dem> class
typedef Dem* DEM;

#endif
