#ifndef __BOUNDARY_CONDITIONS_INCLUDED__
#define __BOUNDARY_CONDITIONS_INCLUDED__
#include <string>
#include "AnalyticSolution.h"
class BoundaryConditions
{
 public:
  BoundaryConditions(const bool analytic_=true,
		     const double time_=-1000.,
		     const double experimental_air_temperature_=-1000.,
		     const double experimental_solar_radiation_=-1000.,
		     const double experimental_wind_speed_=-1000.,
		     const double experimental_relative_humidity_=-1000.,
		     const double experimental_precipitation_=-1000.,
		     const double experimental_surface_temperature_=-1000.,
		     const double experimental_surface_pressure_=-75.2025,
		     const bool moisture_movement_=false);
  //------------------------------------------------------------
  double get_solar_flux        (const std::string surface_type,
				const std::string author,
				const double canopy_density);
  double get_infrared_flux     (const std::string surface_type,
				const std::string author,
				const std::string direction,
				const double old_surface_temperature,
				const double canopy_temperature,
				const double canopy_density);
  double get_convective_flux   (const std::string surface_type,
				const std::string author,
				const std::string direction,
				const double canopy_density,
				const double old_surface_temperature,
				const double new_surface_temperature,
				const bool new_estimate);
  double get_evaporative_flux  (const std::string surface_type,
				const std::string author,
				const double canopy_density,
				const double old_surface_temperature,
				const double new_surface_temperature,
				const bool new_estimate
				/*const bool linear_evaporative*/);
  double get_canopy_temperature(/*const std::string surface_type,*/
				const std::string author,
				const double canopy_temperature);
  double get_inbound_heat_flux (const std::string surface_type,
				const std::string author,
				const double shading_factor,
				const double canopy_temperature,
				const double canopy_density,
				const double old_surface_temperature,
				const double new_surface_temperature,
				const bool new_estimate
				/*const bool linear_evaporative*/);
  void print_inbound_heat_fluxes (double &inbound_solar_flux,
				  double &inbound_convective_flux,
				  double &inbound_evaporative_flux,
				  double &inbound_infrared_flux,
				  double &outbound_convection_coefficient,
				  double &outbound_infrared_coefficient,
				  double &outbound_evaporative_coefficient);
  double get_outbound_coefficient (const std::string surface_type,
				   const std::string author,
				   const double canopy_density,
				   const double old_surface_temperature,
				   const double new_surface_temperature,
				   const double new_estimate
				   /*const bool linear_evaporative*/);

  double get_evaporative_mass_flux(const std::string surface_type,
				   const std::string author,
				   const double canopy_density,
				   const double old_surface_temperature,
				   const double new_surface_temperature,
				   const bool new_estimate);
  
private:
  const bool   analytic;
  const double time;
  double air_temperature;
  double solar_radiation;
  double wind_speed;
  double relative_humidity;
  double precipitation;
  double surface_temperature;
  double surface_pressure;
  const bool moisture_movement;
 
  double pi;
  double phi;
  double theta;

  double print_solar_radiation;
  double print_inbound_convective;
  double print_inbound_infrared;
  double print_inbound_evaporative;
  double print_outbound_convection_coefficient;
  double print_outbound_infrared_coefficient;
  double print_outbound_evaporative_coefficient;
  //  double print_outbound_evaporative_coefficient;
};
#endif
