#ifndef __SURFACE_COEFFICIENTS_INCLUDED__
#define __SURFACE_COEFFICIENTS_INCLUDED__
//*****************************************************************************************************************
/*
  This function return the heat transfer coefficients for four heat transfer processes:
  solar radiation, air convection, infrarred radiation and evaporation. The return value
  is in W/m2K.
  
  Some of these coefficients depend only on the type of surface (e.g. solar albedo) but
  other depende on the particular meteorological conditions the surface is exposed to
  (e.g evaporative). For this reason the functions take meteorological data as input.

  Of course, different fields and authors offer different heat transfer formulations.
  In this program I'm including three of formulations that I found useful during my
  research. These formulations can be found in the following publications:

  [1] - Jansson, C., Almkvist, E., Jansson, P., 2006. Heat balance of an asphalt 
  surface: observations and physically-based simulations. Meteorological Applications
  13, 203–212. doi:10.1017/S1350482706002179

  [2] - Herb, W.R., Janke, B., Mohseni, O., Stefan, H.G., 2008. Ground surface 
  temperature simulation for different land covers. Journal of Hydrology 356, 327–343. 
  doi:10.1016/j.jhydrol.2008.04.020

  [3] - Best, M.J., 1998. A Model to Predict Surface Temperatures. Boundary-Layer
  Meteorology 88, 279–306. doi:10.1023/A:1001151927113

  An important note. At the time of writing this code, I'm not modelling heat moisture
  diffusion. This makes a bit difficult to accurately give a heat transfer coefficient
  and the heat flux associated to the evaporative process. So far a constant value
  based on experimental mesuarements and literature review is being returned.
  
  The last function retunrs an average heat flux value (W/m2) for evaporation. 
  This is because there is no analytical expression for the behaviour of relative humudity.
*/
#include <iostream>

class SurfaceCoefficients
{
public:
  SurfaceCoefficients (double surface_water_tension_=-75.2025,
		       double precipitation_=0.,
		       bool moisture_movement_=false);
  double get_sky_emissivity (const std::string author,
			     const double air_temperature,
			     const double relative_humidity);
  double Clasius_Clapeyron_saturated_vapour_pressure  (const double temperature);
  double Philip_1957_surface_vapour_pressure  (const double surface_temperature);
  double vapour_pressure_surface_new_estimate (const double old_surface_temperature,
					       const double new_surface_temperature);
  // double linear_inbound_evaporative_flux (const std::string surface_type,
  // 					  const std::string author,
  // 					  const double new_air_temperature,
  // 					  const double new_relative_humidity,
  // 					  const double new_wind_speed,
  // 					  const double old_surface_temperature);
  // double linear_outbound_evaporative_coefficient (const std::string surface_type,
  // 						  const std::string author,
  // 						  const double new_air_temperature,
  // 						  const double new_relative_humidity,
  // 						  const double new_wind_speed,
  // 						  const double old_surface_temperature);
  /* Jansson et al. (2006) */
  double get_absortivity_Jansson                   (const std::string surface_type);
  double get_infrared_inbound_coefficient_Jansson  (const std::string surface_type,
						    const double air_temperature,
						    const double relative_humidity);
  double get_infrared_outbound_coefficient_Jansson (const std::string surface_type);
  double get_convective_coefficient_Jansson        (const std::string surface_type,
						    const double wind_speed);
  double get_evaporative_coefficient_Jansson       (const std::string surface_type,
						    const double wind_speed);
  double get_evaporative_flux_Jansson              (const std::string surface_type,
						    const double air_temperature,
						    const double relative_humidity,
						    const double wind_speed,
						    double old_surface_temperature,
						    double new_surface_temperature,
						    const bool new_estimate);
  double get_evaporative_mass_flux_Jansson         (const std::string surface_type,
						    const double air_temperature,
						    const double relative_humidity,
						    const double wind_speed,
						    const double old_surface_temperature,
						    const double new_surface_temperature,
						    const bool new_estimate);
  /* Herb et al. (2008) */
  double get_absortivity_Herb                  (const std::string surface_type);
  double get_infrared_inbound_coefficient_Herb (const std::string surface_type,
						const double air_temperature,
						const double relative_humidity);
  double get_infrared_outbound_coefficient_Herb(const std::string surface_type);
  double get_convective_coefficient_Herb       (/*const std::string surface_type,*/
						const double air_temperature,
						const double relative_humidity,
						const double wind_speed,
						const double old_surface_temperature,
						const double new_surface_temperature,
						const bool new_estimate);
  double get_evaporative_coefficient_Herb      (/*const std::string surface_type,*/
						const double air_temperature,
						const double relative_humidity,
						const double wind_speed,
						const double old_surface_temperature,
						const double new_surface_temperature,
						const bool new_estimate);
  double get_evaporative_flux_Herb             (/*const std::string surface_type,*/
						const double air_temperature,
						const double relative_humidity,
						const double wind_speed,
						double old_surface_temperature,
						double new_surface_temperature,
						const bool new_estimate);
  double get_evaporative_mass_flux_Herb        (/*const std::string surface_type,*/
						const double air_temperature,
						const double relative_humidity,
						const double wind_speed,
						const double old_surface_temperature,
						const double new_surface_temperature,
						const bool new_estimate);
  /* Best 1998 */
  double get_absortivity_Best (const std::string surface_type);
  double get_infrared_inbound_coefficient_Best(const std::string surface_type,
					       const double air_temperature,
					       const double relative_humidity);
  double get_infrared_inbound_coefficient_canopy_sky_Best(const std::string surface_type,
							  const double air_temperature,
							  const double relative_humidity);
  double get_infrared_inbound_coefficient_canopy_soil_Best(const std::string surface_type);
  double get_infrared_outbound_coefficient_Best(const std::string surface_type);
  double get_convective_coefficient_Best (/*const std::string surface_type,*/
					  const double air_temperature,
					  const double relative_humidity,
					  const double wind_speed,
					  const double old_surface_temperature,
					  const double new_surface_temperature,
					  const bool new_estimate);
  double get_evaporative_coefficient_Best(/*const std::string surface_type,*/
					  const double air_temperature,
					  const double relative_humidity,
					  const double wind_speed,
					  const double old_surface_temperature,
					  const double new_surface_temperature,
					  const bool new_estimate);
  double get_evaporative_flux_Best       (/*const std::string surface_type,*/
					  const double air_temperature,
					  const double relative_humidity,
					  const double wind_speed,
					  double old_surface_temperature,
					  double new_surface_temperature,
					  const bool new_estimate);
  double get_evaporative_mass_flux_Best  (/*const std::string surface_type,*/
					  const double air_temperature,
					  const double relative_humidity,
					  const double wind_speed,
					  const double old_surface_temperature,
					  const double new_surface_temperature,
					  const bool new_estimate);
  double get_convective_coefficient_canopy_Best(/*const std::string surface_type,*/
						const double wind_speed);
  double get_evaporative_flux_canopy_Best(/*const std::string surface_type,*/
					  const double air_temperature,
					  const double relative_humidity,
					  const double wind_speed,
					  const double solar_radiation,
					  const double surface_temperature,
					  const bool derivative);
private:
  double surface_water_tension;
  double precipitation;
  bool moisture_movement;
  double molar_mass_water;
  double molar_mass_air;
  double gravity_constant;
  double gas_constant;
  double steffan_boltzmann;
  double von_karman; 
  
  double air_density;
  double air_specific_heat_capacity;
  double atmospheric_pressure;
  
  double water_latent_heat_of_vaporization;
  
  double equilibrium_temperature ;  // estimated equilibrium temperature at the bottom of the soil
  
  double absortivity_soil_Jansson;
  double absortivity_road_Jansson;
  
  double absortivity_soil_Herb;
  double absortivity_road_Herb;

  double absortivity_soil_Best;
  double absortivity_road_Best;

  double absortivity_snow;

  double emissivity_soil_Jansson;
  double emissivity_road_Jansson;
  
  double emissivity_soil_Herb;
  double emissivity_road_Herb;

  double emissivity_soil_Best;
  double emissivity_road_Best;
  double emissivity_canopy_Best;

  double emissivity_snow;
  
  double a_Jansson;
  double b_Jansson;

  double z_road_momentum_roughness;
  double z_road_heat_roughness;
  
  double z_soil_momentum_roughness;
  double z_soil_heat_roughness;
  
  double z_ref_wind;
  double z_ref_temperature;
  
  double coefficient_forced_convection;
  double coefficient_natural_convection;
  double coefficient_sheltering;

  double psychrometric_constant;
};

#endif
