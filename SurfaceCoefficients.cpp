#include "SurfaceCoefficients.h"
#include <math.h>

SurfaceCoefficients::SurfaceCoefficients (double surface_water_tension_,
					  double precipitation_,
					  bool moisture_movement_)
{
  //surface_water_tension = 430 * 0.1;     // constant, (m) EXPERIMENTAL
  //surface_water_tension = -75.2025;      // constant, (m) Garrat 1994
  if (moisture_movement==false)
    surface_water_tension=-75.2025;
  else
    surface_water_tension=surface_water_tension_;
  precipitation         = precipitation_;
  moisture_movement     = moisture_movement_;
  molar_mass_water      = 0.0180153;     // constant, (kg/mol)
  molar_mass_air        = 0.02897;       // constant, (kg/mol)
  gravity_constant      = 9.81;          // constant, (m/s2)
  gas_constant          = 8.3144621;     // constant, (J/molK)
  steffan_boltzmann     = 5.67E-8;       // constant, (W/m2K4)
  von_karman            = 0.41;          // constant, (dimensionless)

  air_density = 1.2041;                        // (kg/m3)
  air_specific_heat_capacity = 1012;           // (J/kgK)
  water_latent_heat_of_vaporization = 2.26E6;  // (J/kg)

  atmospheric_pressure              = 101325.; // (Pa)  
  
  absortivity_soil_Jansson = 0.85;         // Garrat 1994,  (dimensionless)
  absortivity_road_Jansson = 0.90;         // Jansson 2006, (dimensionless)
  
  absortivity_soil_Herb = 0.85;            // Herb 2008, (dimensionless)
  absortivity_road_Herb = 0.88;            // Herb 2008, (dimensionless) //0.88

  absortivity_soil_Best = 0.85;            // Herb 2008, (dimensionless)
  absortivity_road_Best = 0.88;            // Herb 2008, (dimensionless)

  absortivity_snow = 0.15;                 // Garrat 1994, fresh snow 0.05 - 0.35

  emissivity_soil_Jansson = 0.97;          // Garrat 1994,  (dimensionless)
  emissivity_road_Jansson = 1.00;          // Jansson 2006, (dimensionless)

  emissivity_soil_Herb = 0.95;             // Herb 2008, (dimensionless)
  emissivity_road_Herb = 0.94;             // Herb 2008, (dimensionless) //0.94
  
  emissivity_soil_Best = 0.95;             // Herb 2008, (dimensionless)
  emissivity_road_Best = 0.94;             // Herb 2008, (dimensionless)
  emissivity_canopy_Best = 0.95;           // Herb 2008, (dimensionless)

  emissivity_snow=0.95;                    // Garrat 1994
  
  a_Jansson = 0.540;                       // Monteith 1961, Brunt 1932, (dimensionless)
  b_Jansson = 0.065;                       // Monteith 1961, Brunt 1932, (1/hPa)
  
  z_road_momentum_roughness = 5E-4; // (m)
  z_road_heat_roughness     = 1E-4; // (m)
  
  z_soil_momentum_roughness = 1E-3; // (m)
  z_soil_heat_roughness     = 1E-3; // (m)
  
  z_ref_wind        = 3.0;          // (m)
  z_ref_temperature = 3.0;          // (m)
  
  coefficient_forced_convection  = 0.0015; // (dimensionless) nominal value according to Herb et al 2008 0.0015
  coefficient_natural_convection = 0.0015; // (m/s  1/K^0.33) nominal value according to Herb et al 2008 0.0015
  coefficient_sheltering         = 1.0000; // (dimensionless)

  psychrometric_constant =
    (air_specific_heat_capacity * atmospheric_pressure /
     ((molar_mass_water/molar_mass_air) * water_latent_heat_of_vaporization)); // (Pa/K)
}

double SurfaceCoefficients::
get_sky_emissivity (const std::string author,
		    const double air_temperature,
		    const double relative_humidity)
{
  double vapour_pressure_air=
    (relative_humidity/100.)
    *Clasius_Clapeyron_saturated_vapour_pressure(air_temperature);
  double sky_emissivity=0;
  
  if (author=="Jansson")//Brunt 1932 (vapor pressure must be in hPa)
    {
      sky_emissivity=
	a_Jansson +
	b_Jansson*pow((vapour_pressure_air/100.),0.5); 
    }
  else if (author=="Herb" ||
	   author=="Best")
    {
      double coefficient_cloud_cover=0.5;
      /*
	Edinger 1956? Apparently the vapor pressure 
	must be in hPa (requires confirmation, I
	haven't been able to locate the reference.
	It is cited in Herb et al 2008.
      */
      sky_emissivity=
	coefficient_cloud_cover
	+ 0.67*(1-coefficient_cloud_cover)*
	pow(vapour_pressure_air/100.,0.08); 
    }
  else
    {
      std::cout << "Error: author not implemented." << std::endl
		<< "Error in get_sky_emissivity."   << std::endl;
      throw 100;
    }
  return (sky_emissivity);
}

double SurfaceCoefficients::
Clasius_Clapeyron_saturated_vapour_pressure(const double temperature)
{
  /*
   * Clausius-Clapeyron equation
   * Relates equilibrium or saturation vapor pressure and temperature to the
   * latent heat of the phase change, without requiring specific volume data
   */
  return(611.*
	 exp((water_latent_heat_of_vaporization*molar_mass_water/gas_constant)*//[Pa]
	     ((1./273.15)-(1./(temperature+273.15)))));
}

double SurfaceCoefficients::
Philip_1957_surface_vapour_pressure(const double surface_temperature)
{
  double vapour_pressure_saturation_surface=
    Clasius_Clapeyron_saturated_vapour_pressure(surface_temperature);
  double vapour_pressure_surface=
    vapour_pressure_saturation_surface*
    exp(surface_water_tension*molar_mass_water*gravity_constant
	/(gas_constant*(surface_temperature+273.15)));//(Pa)
  return (vapour_pressure_surface);
}

double SurfaceCoefficients::
vapour_pressure_surface_new_estimate(const double old_surface_temperature,
				     const double new_surface_temperature)
{
  double old_vapour_pressure_surface=
    Philip_1957_surface_vapour_pressure(old_surface_temperature);
  double constant_1=
    (surface_water_tension*molar_mass_water*gravity_constant/gas_constant)
    -(water_latent_heat_of_vaporization*molar_mass_water/gas_constant);//[K]
  
  double average_temperature=
    0.5*old_surface_temperature+0.5*new_surface_temperature;
  double old_vapur_pressure_surface_derivative=
    -1.*Philip_1957_surface_vapour_pressure(average_temperature)
    *(constant_1/((average_temperature+273.15)*(average_temperature+273.15)));
  
  return (old_vapour_pressure_surface
	  + (old_vapur_pressure_surface_derivative
	     *(new_surface_temperature-old_surface_temperature)));
}

// double SurfaceCoefficients::
// linear_inbound_evaporative_flux (const std::string surface_type,
// 				 const std::string author,
// 				 const double new_air_temperature,
// 				 const double new_relative_humidity,
// 				 const double new_wind_speed,
// 				 const double old_surface_temperature)
// {
//   double inbound_evaporative_coefficient;
//   if (author=="Jansson")
//     {
//       inbound_evaporative_coefficient = 
// 	get_evaporative_coefficient_Jansson (surface_type,
// 					     new_wind_speed);
//     }
//   else if (author=="Herb")
//     {
//       inbound_evaporative_coefficient = 
// 	get_evaporative_coefficient_Herb (surface_type,
// 					  new_air_temperature,
// 					  new_relative_humidity,
// 					  new_wind_speed,
// 					  old_surface_temperature);
//     }
//   else if (author=="Best")
//     {
//       inbound_evaporative_coefficient = 
// 	get_evaporative_coefficient_Herb (surface_type,
// 					  new_air_temperature,
// 					  new_relative_humidity,
// 					  new_wind_speed,
// 					  old_surface_temperature);
//     }
//   else
//     {
//       std::cout << "Error in linear_inbound_evaporative_flux\n"
// 		<< "Author not implemented\n";
//     }
  
//   double vapor_pressure_air =
//     (new_relative_humidity / 100.)
//     * Clasius_Clapeyron_saturated_vapour_pressure (new_air_temperature) ; // [Pa]
//   double vapor_pressure_surface 
//     = Philip_1957_surface_vapour_pressure (old_surface_temperature);
//   double constant_1 
//     = (surface_water_tension*molar_mass_water*gravity_constant/gas_constant)
//     + (water_latent_heat_of_vaporization*molar_mass_water/gas_constant);  // [K]
  
//   return (inbound_evaporative_coefficient*vapor_pressure_air
// 	  - (inbound_evaporative_coefficient*vapor_pressure_surface*
// 	     (1.-(constant_1*old_surface_temperature
// 		  /((old_surface_temperature+273.15)*(old_surface_temperature+273.15))))));
// }

// double SurfaceCoefficients::
// linear_outbound_evaporative_coefficient (const std::string surface_type,
// 					 const std::string author,
// 					 const double new_air_temperature,
// 					 const double new_relative_humidity,
// 					 const double new_wind_speed,
// 					 const double old_surface_temperature)
// {
//   double inbound_evaporative_coefficient;
//   if (author=="Jansson")
//     {
//       inbound_evaporative_coefficient = 
// 	get_evaporative_coefficient_Jansson (surface_type,
// 					     new_wind_speed);
//     }
//   else if (author=="Herb")
//     {
//       inbound_evaporative_coefficient = 
// 	get_evaporative_coefficient_Herb (surface_type,
// 					  new_air_temperature,
// 					  new_relative_humidity,
// 					  new_wind_speed,
// 					  old_surface_temperature);
//     }
//   else if (author=="Best")
//     {
//       inbound_evaporative_coefficient = 
// 	get_evaporative_coefficient_Herb (surface_type,
// 					  new_air_temperature,
// 					  new_relative_humidity,
// 					  new_wind_speed,
// 					  old_surface_temperature);
//     }
//   else
//     {
//       std::cout << "Error in linear_inbound_evaporative_flux\n"
// 		<< "Author not implemented\n";
//     }
  
//   double vapor_pressure_surface 
//     = Philip_1957_surface_vapour_pressure (old_surface_temperature);
//   double constant_1 
//     = (surface_water_tension*molar_mass_water*gravity_constant/gas_constant)
//     + (water_latent_heat_of_vaporization*molar_mass_water/gas_constant);  // [K]
  
//   return (inbound_evaporative_coefficient*constant_1*vapor_pressure_surface
// 	  /((old_surface_temperature+273.15)*(old_surface_temperature+273.15)));
// }
/*********************************************************************************/
/***************                                                   ***************/
/***************      Jansson et al. (2006) Surface Coefficients   ***************/
/***************                                                   ***************/
/*********************************************************************************/
double SurfaceCoefficients::
get_absortivity_Jansson(const std::string surface)
{
  double coefficient=0;
  if (surface=="road")
    coefficient=absortivity_road_Jansson;
  else if (surface=="soil")
    coefficient=absortivity_soil_Jansson;
  else
    {
      std::cout << "Error: surface not implemented."   << std::endl
		<< "Error in get_absortivity_Jansson." << std::endl;
      throw 3;
    }
  return(coefficient);
}

double SurfaceCoefficients::
get_infrared_inbound_coefficient_Jansson(const std::string surface_type,
					 const double air_temperature,   // (C)
					 const double relative_humidity) // (%)
{
  double sky_emissivity
    =get_sky_emissivity ("jansson",
			 air_temperature,
			 relative_humidity);
  double surface_emissivity=0.;
  if (surface_type=="road")
    surface_emissivity=emissivity_road_Jansson;
  else if (surface_type=="soil")
    surface_emissivity=emissivity_soil_Jansson;
  else
    {
      std::cout << "Error: surface not implemented."                    << std::endl
		<< "Error in get_infrared_inbound_coefficient_Jansson." << std::endl;
      throw 3;
    }
  return (steffan_boltzmann*sky_emissivity*surface_emissivity); // (W/m2K4)
}

double SurfaceCoefficients::
get_infrared_outbound_coefficient_Jansson(const std::string surface_type)
{
  double surface_emissivity=0.;
  if (surface_type=="road")
    surface_emissivity=emissivity_road_Jansson;
  else if (surface_type=="soil")
    surface_emissivity=emissivity_soil_Jansson;
  else
    {
      std::cout << "Error: surface not implemented."                     << std::endl
		<< "Error in get_infrared_outbound_coefficient_Jansson." << std::endl;
      throw 3;
    }
  return (steffan_boltzmann*surface_emissivity); // (W/m2K4)
}

double SurfaceCoefficients::
get_convective_coefficient_Jansson(const std::string surface_type,
				   const double wind_speed)
{
  double resistance=0;
  if (surface_type=="road")
    resistance = ((1./(pow(von_karman,2)*wind_speed)) *
		  log10(z_ref_wind/z_road_momentum_roughness) *
		  log10(z_ref_temperature/z_road_heat_roughness));
  else if (surface_type=="soil")
    resistance = ((1./(pow(von_karman,2)*wind_speed)) *
		  log10(z_ref_wind/z_soil_momentum_roughness) *
		  log10(z_ref_temperature/z_soil_heat_roughness));
  else
    {
      std::cout << "Error: surface not implemented."      << std::endl
		<< "Error in get_convective_coefficient." << std::endl;
      throw 3;
    }
  return(air_density*air_specific_heat_capacity/resistance);
}

double SurfaceCoefficients::
get_evaporative_coefficient_Jansson(const std::string surface_type,
				    const double wind_speed)
{
  /*
    The evaporative coeffiecient in Jansson's formulation (turbulent
    formulation) is identical to the convective coefficient divided
    by the psychrometric constant
  */
  return (get_convective_coefficient_Jansson(surface_type,
					     wind_speed)/psychrometric_constant);
}

double SurfaceCoefficients::
get_evaporative_flux_Jansson(const std::string surface_type,
			     const double air_temperature  ,
			     const double relative_humidity,
			     const double wind_speed       ,
			     double old_surface_temperature,
			     double new_surface_temperature,
			     const bool new_estimate)
{
  if (moisture_movement==false)
    {
      old_surface_temperature=10.9;
      new_surface_temperature=10.9;
    }

  if (surface_type=="soil")
    {   
      /*
	The equations used for (Jansson et al 2006) for evaporation
	from a bare soil assume a saturated soil. This is the same
	assumption as in the equations used by (Herb et al 2008).
	This equation does not include a term to take into account
	the contribution of natural convection. For this reason, for
	low wind speeds (more usual that one would expect according
	to the experimental data) the convective and evaporative heat
	fluxes are negligible, and this in turn bring the soil to very
	high temperatures. However, the advantage of this formulation
	is that it has been modified by (Philip et al 1957) to be 
	applied to soils by including a 'soil moisture availability
	factor'
	From thermodynamics laws Philip (1957) dereive an
	expression for the relative humidity of air in equilibrium
	with the water in the soil pore. A prerequisite for the
	equation to be valid is that the air close to the pore
	water surface is in equilibrium with the pore water.
	However, in a drying soil, where vapor is continously
	transported to the atmosphere, the vapor pressure of air
	adjacent to the evaporating water surface will not be in
	equilibrium with the liquid water in the pore, i.e. the
	calculated relative humidity won't be representative at the
	soil surface.
      */
      double vapour_pressure_surface=0.;
      if (new_estimate==false)
	vapour_pressure_surface
	  =Philip_1957_surface_vapour_pressure (old_surface_temperature);
      else
	vapour_pressure_surface
	  =vapour_pressure_surface_new_estimate (old_surface_temperature,
						 new_surface_temperature);
      /*
	An obvious drawback of this formulation is that we need
	to know the surface water tension at the surface (or near
	the surface) in order ot use it and unless we have
	experimental measurements or we are dealing with a fully
	coupled heat-moisture diffusion model, it might prove to
	be difficult to give a reasonable number and even more
	difficult to propose its variation in time.
      */
      double evaporative_coefficient 
	=get_evaporative_coefficient_Jansson(surface_type,
					     wind_speed);
      double vapour_pressure_air
	=(relative_humidity/100.)
	* Clasius_Clapeyron_saturated_vapour_pressure (air_temperature) ; // (Pa)
      
      return(evaporative_coefficient 
	     *(vapour_pressure_air-vapour_pressure_surface)); // (W/m2)
    }
  else if (surface_type=="road")
    {
      /*
	Assumption: No evaporation from the pavement surface.
	This means that we are considering it as a nonporous
	surface. In principle, if there is rainfall, the
	surface will get wet and some evaporation will occur.
	So, we are also assuming that the rainfall events are
	scarce. This assumption, however, might be wrong as,
	in my experience, heavy rainfalls in the uk are uncommon,
	drizle is a lot more usual.
      */
      return 0.; // (W/m2)
    }
  else
    {
      std::cout << "Error: surface not implemented." << std::endl
		<< "Error in get_evaporative."       << std::endl;
      throw 3;
      return -1;
    }
}

double SurfaceCoefficients::
get_evaporative_mass_flux_Jansson(const std::string surface_type,
				  const double air_temperature  ,
				  const double relative_humidity,
				  const double wind_speed       ,
				  const double old_surface_temperature,
				  const double new_surface_temperature,
				  const bool new_estimate)
{
  return ((1./(/*air_density*/water_latent_heat_of_vaporization))*
	  get_evaporative_flux_Jansson (surface_type,air_temperature,
					relative_humidity,wind_speed,
					old_surface_temperature,new_surface_temperature,
					new_estimate));
}
/******************************************************************************************************/
/***********************                                                    ***************************/
/***********************      Herb et al. (2008) Surface Coefficients       ***************************/
/***********************                                                    ***************************/
/******************************************************************************************************/
double SurfaceCoefficients::
get_absortivity_Herb(const std::string surface)
{
  double coefficient=0;
  if (surface=="road")
    coefficient=absortivity_road_Herb;
  else if (surface=="soil")
    coefficient=absortivity_soil_Herb;
  else if (surface=="snow")
    coefficient=absortivity_snow;
  else
    {
      std::cout << "Error: surface not implemented." << std::endl
		<< "Error in get_absortivity_Herb."  << std::endl;
      throw 3;
    }
  return coefficient;
}

double SurfaceCoefficients::
get_infrared_inbound_coefficient_Herb(const std::string surface_type,
				      const double air_temperature,   // (C)
				      const double relative_humidity) // (%)
{
  double sky_emissivity
    =get_sky_emissivity ("Herb",
			 air_temperature,
			 relative_humidity);
  double surface_emissivity=0.;
  if (surface_type=="road")
    surface_emissivity=emissivity_road_Herb;
  else if (surface_type=="soil")
    surface_emissivity=emissivity_soil_Herb;
  else
    {
      std::cout << "Error: surface not implemented."         << std::endl
		<< "Error in get_infrared_coefficient_Herb." << std::endl;
      throw 3;
    }
  return(steffan_boltzmann*sky_emissivity*surface_emissivity); // (W/m2K4)
}

double SurfaceCoefficients::
get_infrared_outbound_coefficient_Herb(const std::string surface_type)
{
  double surface_emissivity=0.;
  if (surface_type=="road")
    surface_emissivity=emissivity_road_Herb;
  else if (surface_type=="soil")
    surface_emissivity=emissivity_soil_Herb;
  else if (surface_type=="snow")
    surface_emissivity=emissivity_snow;
  else
    {
      std::cout << "Error: surface not implemented."                  << std::endl
		<< "Error in get_infrared_outbound_coefficient_Herb." << std::endl;
      throw 3;
    }
  return (steffan_boltzmann*surface_emissivity);//[W/m2K4]
}

double SurfaceCoefficients::
get_convective_coefficient_Herb(/*const std::string surface_type,*/
				const double air_temperature,     // (C)
				const double relative_humidity,   // (%)
				const double wind_speed,          // (m/s)
				const double old_surface_temperature, // (C)
				const double new_surface_temperature, // (C)
				const bool new_estimate) 
{
  /*
    "The advantage of this formulation is that both mechanical induced turbulence
    and atmospheric convective transport are included"
    [Henderson 1981, "Modelling evaporation from reservoirs"] 
    ...
    The main assumption in this formulation is that the air near the surface (soil
    or water reservoir) is saturated at all times. However:
    ...
    "When an early saturated soil dries, the free water in larger pores evaporates
    first. The remaining water is retained mainly in smaller soil pores by strong
    capillary forces, and the vapour pressure of the air in equilibrium with the
    pore water will therefore be lower than in the air close to a free water surface"
    [Alvenas and Janssen 1997, "Model for evaporation, moisture and temperature"]
    ...
    This means that the assumption of a soil saturated at all times should be revised
    if we are in weather conditions that are likely to dry the soil.
    For the road case the situation is different. We assume the road surface is dry
    except when a rain event occurs. This should be decided by the function calling
    this member.
    
    Clausius-Clapeyron equation
    Relates equilibrium or saturation vapour pressure and temperature to the
    latent heat of the phase change, without requiring specific volume data
  */
  /*
    the 'old_surface_temperature' passed to this function is assumed to be an
    'anchor temperature', this is, a temperature around which we expect the 
    surface temperature to be around.
  */
  double vapour_pressure_surface=0.;
  double mixing_ratio_surface=0.;
  double virtual_surface_temperature=0.;
  double surface_temperature=//[K]
    0.5*new_surface_temperature+0.5*old_surface_temperature;
  if (new_estimate==false)
    {
      vapour_pressure_surface=
	Philip_1957_surface_vapour_pressure(surface_temperature);
      mixing_ratio_surface=
	(molar_mass_water/molar_mass_air)*
	vapour_pressure_surface/atmospheric_pressure;//[Pa/Pa]
      virtual_surface_temperature=
	(1.+0.6*mixing_ratio_surface)*(surface_temperature+273.15);//[K]
    }
  else
    {
      vapour_pressure_surface=
	vapour_pressure_surface_new_estimate(surface_temperature,
					     surface_temperature);
      mixing_ratio_surface=
	(molar_mass_water/molar_mass_air)
	*vapour_pressure_surface/atmospheric_pressure;//[Pa/Pa]
      virtual_surface_temperature=
	(1.+0.6*mixing_ratio_surface)*(surface_temperature+273.15);//[K]
    }
  double saturated_vapour_pressure_air=
    Clasius_Clapeyron_saturated_vapour_pressure(air_temperature);
  double mixing_ratio_air=
    ((relative_humidity/100.)*
     (molar_mass_water/molar_mass_air)
     *saturated_vapour_pressure_air/atmospheric_pressure);//[Pa/Pa]
  double virtual_air_temperature=
    (1.+0.6*mixing_ratio_air)*(air_temperature+273.15);//[K]
  /*
    Ryan et al. equation for flux evaporation includes a term for natural convection.
    This term is proportional to the difference in virtual temperatures between the
    surface and air. This term is raised to a fractional power (1/3). For negative
    numbers (i.e. when the virtual temperature of the air is higher than the virtual
    temperature of the surface) would need the computation of complex numbers. However
    for that case, it doesn't make sense to include the term for natural convection
    because it would correspond to a stable temperature profile that prevents natural
    convection. In this case, I simply set to zero this term in order to neglect its
    contribution.
  */
  double virtual_temperature_difference=0.;
  if (virtual_surface_temperature>virtual_air_temperature)
    virtual_temperature_difference
      =virtual_surface_temperature-virtual_air_temperature;
  else
    virtual_temperature_difference=0.;
  /*
    Specific humidity (q) is the number of grams of water vapour per unit mass of 
    air plus water vapour. In applications, q is numerically close enough to the
    mixing ratio that one seldom needs to distinguish between them. [North and 
    Erukhimova 2009, "Atmospheric Thermodynamics Elementary Physics and Chemistry"
  */
  return (air_density*air_specific_heat_capacity
	  *(coefficient_forced_convection*coefficient_sheltering*wind_speed + 
	    coefficient_natural_convection*pow(virtual_temperature_difference,1./3.))); // (W/m^2K)
}

double SurfaceCoefficients::
get_evaporative_coefficient_Herb(/*const std::string surface_type,*/
				 const double air_temperature,     // (C)
				 const double relative_humidity,   // (%)
				 const double wind_speed,          // (m/s)
				 const double old_surface_temperature, // (C)
				 const double new_surface_temperature, // (C)
				 const bool new_estimate) // (C)
{
  /*
    The evaporative coeffiecient in Herb's and Best's formulations
    (non-turbulent formulations) is identical to the convective 
    coefficient except that the convective coefficient is multiplied
    by the air's specific heat capacity while the evaporative
    coefficient is multiplied by the water's latent heat of vaporization
  */
  return ((water_latent_heat_of_vaporization/air_specific_heat_capacity)
	  *get_convective_coefficient_Herb(/*surface_type,*/air_temperature,
					   relative_humidity,wind_speed,
					   old_surface_temperature,
					   new_surface_temperature,
					   new_estimate)
	  *((molar_mass_water/molar_mass_air)/atmospheric_pressure));
}

double SurfaceCoefficients::
get_evaporative_flux_Herb(const double air_temperature,  // (C)
			  const double relative_humidity,// (%)
			  const double wind_speed,       // (m/s)
			  double old_surface_temperature,// (C)
			  double new_surface_temperature,// (C)
			  const bool new_estimate)// (C)
{
  /*
    the 'old_surface_temperature' passed to this function is assumed to be an
    'anchor temperature', this is, a temperature around which we expect the 
    surface temperature to be around.
  */
  double surface_temperature=//[K]
    0.5*new_surface_temperature+0.5*old_surface_temperature;
  double vapour_pressure_surface=0.;
  if (new_estimate==false)
    vapour_pressure_surface=
      Philip_1957_surface_vapour_pressure(surface_temperature);
  
  else
    vapour_pressure_surface=
      vapour_pressure_surface_new_estimate(surface_temperature,
					    surface_temperature);
  double vapour_pressure_air=
    (relative_humidity/100.)
    *Clasius_Clapeyron_saturated_vapour_pressure(air_temperature) ;//[Pa]
  
  return(get_evaporative_coefficient_Herb(/*surface_type,*/air_temperature,
					  relative_humidity,wind_speed,
					  old_surface_temperature,
					  new_surface_temperature,
					  new_estimate)
	 *(vapour_pressure_air-vapour_pressure_surface));//[W/m^2]
}

double SurfaceCoefficients::
get_evaporative_mass_flux_Herb(/*const std::string surface_type,*/
			       const double air_temperature,     // (C)
			       const double relative_humidity,   // (%)
			       const double wind_speed,          // (m/s)
			       const double old_surface_temperature, // (C)
			       const double new_surface_temperature, // (C)
			       const bool new_estimate) // (C)
{
  return ((1./(water_latent_heat_of_vaporization/*/air_specific_heat_capacity*/))*
	  get_evaporative_flux_Herb(/*surface_type,*/air_temperature,
				    relative_humidity,wind_speed,
				    old_surface_temperature,new_surface_temperature,
				    new_estimate));
}
/******************************************************************************************************/
/**************************                                                   *************************/
/**************************      Best et al. (2008) Surface Coefficients      *************************/
/**************************                                                   *************************/
/******************************************************************************************************/
double SurfaceCoefficients::
get_absortivity_Best (const std::string surface)
{
  double coefficient=0;
  if (surface=="road")
    coefficient=absortivity_road_Best;
  else if (surface=="soil")
    coefficient=absortivity_soil_Best;
  else
    {
      std::cout << "Error: surface not implemented." << std::endl
		<< "Error in get_absortivity_Best."  << std::endl;
      throw 3;
    }
  return coefficient;
}

double SurfaceCoefficients::
get_infrared_inbound_coefficient_Best(const std::string surface_type,
				      const double air_temperature,
				      const double relative_humidity)
{
  double sky_emissivity = 
    get_sky_emissivity ("Best",
			air_temperature,
			relative_humidity);
  double surface_emissivity = 0.;
  if (surface_type == "road")
    surface_emissivity = emissivity_road_Best;
  else if (surface_type == "soil")
    surface_emissivity = emissivity_soil_Best;
  else if (surface_type == "canopy")
    surface_emissivity = emissivity_canopy_Best;
  else
    {
      std::cout << "Error: surface not implemented."          << std::endl
		<< "Error in get_infrared_coefficient_Best."  << std::endl;
      throw 3;
    }
  return (steffan_boltzmann*sky_emissivity*surface_emissivity); // (W/m2K4)
}

double SurfaceCoefficients::
get_infrared_outbound_coefficient_Best(const std::string surface_type)
{
  double surface_emissivity = 0.;
  if (surface_type == "road")
    surface_emissivity = emissivity_road_Best;
  else if (surface_type == "soil")
    surface_emissivity = emissivity_soil_Best;
  else if (surface_type == "canopy")
    surface_emissivity = emissivity_canopy_Best;
  else
    {
      std::cout << "Error: surface not implemented."                  << std::endl
		<< "Error in get_infrared_outbound_coefficient_Best." << std::endl;
      throw 3;
    }
  return (steffan_boltzmann*surface_emissivity); // (W/m2K4)
}

double SurfaceCoefficients::
get_convective_coefficient_Best(const double air_temperature,     // (C)
				const double relative_humidity,   // (%)
				const double wind_speed,          // (m/s)
				const double old_surface_temperature, // (C)
				const double new_surface_temperature, // (C)
				const bool new_estimate)
{
  /*
    The convective coeffiecients in Herb's and Best's formulations
    (non-turbulent formulations) are identical, so there is no
    necessity to repeat code here, just call Herb's convective
    coefficient function
  */
  return (get_convective_coefficient_Herb(air_temperature,
					  relative_humidity,wind_speed,
					  old_surface_temperature,
					  new_surface_temperature,
					  new_estimate));
}

double SurfaceCoefficients::
get_evaporative_coefficient_Best(const double air_temperature,     // (C)
				 const double relative_humidity,   // (%)
				 const double wind_speed,          // (m/s)
				 const double old_surface_temperature, // (C)
				 const double new_surface_temperature, // (C)
				 const bool new_estimate)
{
  /*
    The evaporative coeffiecients in Herb's and Best's formulations
    (non-turbulent formulations) are identical, so there is no
    necessity to repeat code here, just call Herb's evaporative
    coefficient function
  */ 
  return (get_evaporative_coefficient_Herb(air_temperature,
					   relative_humidity,wind_speed,
					   old_surface_temperature,
					   new_surface_temperature,
					   new_estimate));
}

double SurfaceCoefficients::
get_evaporative_flux_Best(const double air_temperature,     // (C)
			  const double relative_humidity,   // (%)
			  const double wind_speed,          // (m/s)
			  double old_surface_temperature, // (C)
			  double new_surface_temperature, // (C)
			  const bool new_estimate) // (C)
{
  /*
    the 'old_surface_temperature' passed to this function is assumed to be an
    'anchor temperature', this is, a temperature around which we expect the 
    surface temperature to be around.
  */
  double surface_temperature=//[K]
    0.5*new_surface_temperature+0.5*old_surface_temperature;
  double vapour_pressure_surface=0.;
  if (new_estimate==false)
    vapour_pressure_surface=
      Philip_1957_surface_vapour_pressure(surface_temperature);
  else
    vapour_pressure_surface=
      vapour_pressure_surface_new_estimate(surface_temperature,
					   surface_temperature);
  double vapour_pressure_air=
    (relative_humidity/100.)
    *Clasius_Clapeyron_saturated_vapour_pressure(air_temperature);//[Pa]
  
  return(get_evaporative_coefficient_Best(air_temperature,
					  relative_humidity,wind_speed,
					  old_surface_temperature,
					  new_surface_temperature,
					  new_estimate)
	  *(vapour_pressure_air-vapour_pressure_surface));//[W/m^2]
}

double SurfaceCoefficients::
get_convective_coefficient_canopy_Best(const double wind_speed)//(m/s)
{ 
  return (0.01*(wind_speed+0.3)*air_density*air_specific_heat_capacity); // (W/m^2K)
}

double SurfaceCoefficients::
get_evaporative_flux_canopy_Best(const double air_temperature,     // (C)
				 const double relative_humidity,   // (%)
				 const double wind_speed,          // (m/s)
				 const double solar_radiation,     // (W/m2)
				 const double surface_temperature, // (C)
				 const bool derivative) 
{
  double maximum_solar_radiation = 1000.;
  double moisture_content = 0.24;
  double wilting_point = 0.5*moisture_content;
  double stomata_resistance = 200.*(maximum_solar_radiation/(solar_radiation+0.3*maximum_solar_radiation) 
				    + pow(wilting_point/moisture_content,2));
  
  double saturated_vapor_partial_pressure_surface=0.;

  //  if (moisture_movement==true)
    saturated_vapor_partial_pressure_surface
      =Clasius_Clapeyron_saturated_vapour_pressure(surface_temperature);
  // else
  //   saturated_vapor_partial_pressure_surface
  //     =Clasius_Clapeyron_saturated_vapour_pressure(/*surface_temperature*/10.9);

  double saturated_vapor_partial_pressure_air 
    = Clasius_Clapeyron_saturated_vapour_pressure (air_temperature);
  double mixing_ratio_surface =
    (molar_mass_water/molar_mass_air)
    * saturated_vapor_partial_pressure_surface/atmospheric_pressure; // (Pa/Pa)
  double mixing_ratio_air =
    ((relative_humidity/100.) *
     (molar_mass_water/molar_mass_air)
     * saturated_vapor_partial_pressure_air/atmospheric_pressure);     // (Pa/Pa)
  
  if (derivative == false)
    return (air_density * water_latent_heat_of_vaporization 
	    * (0.01*(wind_speed+0.3)/(1.+stomata_resistance*0.01*(wind_speed+0.3)))
	    * (mixing_ratio_air - mixing_ratio_surface) ); // (W/m^2)
  else 
    return(-1. * air_density * water_latent_heat_of_vaporization 
	   * (0.01*(wind_speed+0.3)/(1.+stomata_resistance*0.01*(wind_speed+0.3)))
	   * mixing_ratio_surface * ((water_latent_heat_of_vaporization*molar_mass_water/gas_constant)
				     /pow(surface_temperature+273.15,2)));
  
}

double SurfaceCoefficients::
get_evaporative_mass_flux_Best(const double air_temperature,     // (C)
			       const double relative_humidity,   // (%)
			       const double wind_speed,          // (m/s)
			       const double old_surface_temperature, // (C)
			       const double new_surface_temperature, // (C)
			       const bool new_estimate) // (C)
{
  return ((1./(water_latent_heat_of_vaporization/*/air_specific_heat_capacity*/))*
	  get_evaporative_flux_Best(/*surface_type,*/air_temperature,
				    relative_humidity,wind_speed,
				    old_surface_temperature,new_surface_temperature,
				    new_estimate));
}
