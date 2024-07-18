/* observer.cpp
Function definitions for Observer methods
*/

#include "observer.h"

int Observer::i;
LOOP Observer::loop;
DEM Observer::dem;

Observer::Observer(LOOP loop_object,DEM dem_object)
{
  // Initialize counter
  i = 0;
  // Set needed objects
  loop = loop_object;
  dem = dem_object;
}

Observer::~Observer(void)
{
  // Destroy things here
}

void Observer::Observe(const state_type &state, const double time)
{
  // Store state
  loop->SetState(state);
  // Save terms
  if(loop->parameters.save_terms)
  {
    loop->SaveTerms();
  }
  // Calculate DEM
  if(loop->parameters.calculate_dem)
  {
    dem->CalculateDEM(i);
  }
  // Save results
  loop->SaveResults(i,time);
  // Increment counter
  i++;
}

int Observer::CheckNan(state_type &state, double &time, double &tau, double old_time, double old_tau)
{
  // Check for NaNs in the state
  for(int j=0; j<state.size(); j++)
  {
    if(std::isnan(state[j]))
    {
      // Reset and fail if NaNs found
      time = old_time;
      tau = old_tau*loop->parameters.adaptive_solver_safety;
      state = loop->GetState();
      return 1;
    }
  }

  // Pass otherwise
  return 0;
}

int Observer::CheckNan(state_type &state)
{
  // Check for NaNs in the state
  for(int j=0; j<state.size(); j++)
  {
    if(std::isnan(state[j]))
    {
      return 1;
    }
  }

  // Pass otherwise
  return 0;
}
