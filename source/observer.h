/* observer.h
Class definition for observer class
*/

#ifndef OBSERVER_H
#define OBSERVER_H

#include "helper.h"
#include "loop.h"
#include "dem.h"

// Observer object
//
// The Observer watches the integrator routine and is called at every timestep
// to do something. Here, the observer updates the state and saves any necessary
// results.
//
class Observer {
private:
  /* Timestep counter for printing results */
  static int i;
  /* <Loop> object, used for calling save method */
  static LOOP loop;
  /* <Dem> object, used for calling save method */
  static DEM dem;
public:
  // Default constructor
  // @loop <Loop> instance used for saving loop results
  // @dem <Dem> instance used for saving emission measure results
  //
  // Class for monitoring the integration routine. This object includes methods
  // for watching the integration and saving any needed parameters at each timestep.
  //
  Observer(LOOP loop,DEM dem);

  // Destructor
  ~Observer(void);

  // Observer for the integrator
  // @state current state of the loop system
  // @time current time
  //
  // Method called at each step in the integration. It calls methods
  // from <Loop> and <Dem> to save relevant results. 
  //
  static void Observe(const state_type &state, const double time);
};
// Pointer to the <Observer> class
typedef Observer* OBSERVER;

#endif
