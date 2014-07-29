#ifndef __SURFACE_COEFFICIENTS_INCLUDED__
#define __SURFACE_COEFFICIENTS_INCLUDED__

//*****************************************************************************************************************

// This function return the heat transfer coefficients for three heat transfer processes: solar radiation, air
// convection and infrarred radiation. The return value is in W/m2K. The calculations are based on average mete-
// orological values from experimental data recorded by trl for 2005 and 2006.
// The last function retunrs an average heat flux value (W/m2) for evaporation. This is because there is no ana-
// lytical expression for the behaviour of relative humudity.

#include <iostream>
#include <math.h>

class SurfaceCoefficients
{
public:
  SurfaceCoefficients ();
  /* Jansson et al. (2006) */
  double get_absortivity_Jansson (const std::string surface);
  
  double get_infrared_inbound_coefficient_Jansson  (const std::string surface,
						    const double relative_humidity,
						    const double air_temperature);
  
  double get_infrared_outbound_coefficient_Jansson (const std::string surface);
  
  double get_infrared_coefficient_Jansson (const std::string surface,
					   const double air_temperature,
					   const double relative_humidity,
					   const double previous_air_temperature,
					   const double previous_temperature_surface);
  
  double get_convective_coefficient_Jansson (const std::string surface,
					     const double wind_speed);
  
  double get_evaporative_flux_Jansson (const std::string surface,
				       const double air_temperature,
				       const double wind_speed,
				       const double relative_humidity,
				       const double surface_temperature);
  
  double evaporative_inbound_coefficient  (const std::string surface_type,
					   const double wind_speed,
					   const double relative_humidity,
					   const double air_temperature);
  double evaporative_outbound_coefficient (const std::string surface_type,
					   const double wind_speed,
					   const double surface_temperature);  
  
  /* Herb et al. (2008) */
  double get_absortivity_Herb (const std::string surface);

  double get_infrared_inbound_coefficient_Herb(const std::string surface_type,
					       const double relative_humidity,
					       const double air_temperature);
  
  double get_infrared_outbound_coefficient_Herb(const std::string surface_type);
  
  double get_infrared_coefficient_Herb (const std::string surface,
					const double air_temperature,
					const double relative_humidity,
					const double previous_air_temperature,
					const double previous_temperature_surface);
  
  double get_convective_coefficient_Herb(const std::string surface_type,
					 const double wind_speed,
					 const double relative_humidity,
					 const double air_temperature,
					 const double temperature_surface);

  double get_evaporative_flux_Herb(const std::string surface_type,
				   const double air_temperature,
				   const double wind_speed,
				   const double relative_humidity,
				   const double temperature_surface);
  /* Best 1998 */
  double get_absortivity_Best (const std::string surface);

  double get_infrared_inbound_coefficient_Best(const std::string surface_type,
					       const double relative_humidity,
					       const double air_temperature);
  
  

  double get_infrared_inbound_coefficient_canopy_sky_Best(const std::string surface_type,
							  const double relative_humidity,
							  const double air_temperature);
  
  double get_infrared_inbound_coefficient_canopy_soil_Best(const std::string surface_type);
  
  double get_infrared_outbound_coefficient_Best(const std::string surface_type);
  
  double get_convective_coefficient_Best(const std::string surface_type,
					 const double wind_speed,
					 const double relative_humidity,
					 const double air_temperature,
					 const double temperature_surface);
  
  double get_evaporative_flux_Best(const std::string surface_type,
				   const double air_temperature,
				   const double wind_speed,
				   const double relative_humidity,
				   const double new_temperature_surface);

  double get_convective_coefficient_canopy_Best(const std::string surface_type,
						const double wind_speed);

  double get_evaporative_flux_canopy_Best(const std::string surface_type,
					  const double air_temperature,
					  const double wind_speed,
					  const double relative_humidity,
					  const double solar_radiation,
					  const double temperature_surface,
					  const bool derivative);

private:
  double surface_water_tension;
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

  double emissivity_soil_Jansson;
  double emissivity_road_Jansson;
  
  double emissivity_soil_Herb;
  double emissivity_road_Herb;

  double emissivity_soil_Best;
  double emissivity_road_Best;
  double emissivity_canopy_Best;

  
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
};

SurfaceCoefficients::SurfaceCoefficients ()
{
  //surface_water_tension = 430 * 0.1;     // constant, (m) EXPERIMENTAL
  surface_water_tension = 75.2025;       // constant, (m) Garrat 1994
  molar_mass_water      = 0.0180153;     // constant, (kg/mol)
  molar_mass_air        = 0.02897;       // constant, (kg/mol)
  gravity_constant      = 9.81;          // constant, (m/s2)
  gas_constant          = 8.3144621;     // constant, (J/molK)
  steffan_boltzmann     = 5.67E-8;       // constant, (W/m2K4)
  von_karman            = 0.41;          // constant, (dimensionless)

  air_density = 1.2041;                        // (kg/m3)
  air_specific_heat_capacity = 1012;           // (J/kgK)
  water_latent_heat_of_vaporization = 2.45E6;  // (J/kg) 

  atmospheric_pressure              = 101325.; // (Pa)  
  
  absortivity_soil_Jansson = 0.85;         // Garrat 1994,  (dimensionless)
  absortivity_road_Jansson = 0.90;         // Jansson 2006, (dimensionless)
  
  absortivity_soil_Herb = 0.85;            // Herb 2008, (dimensionless)
  absortivity_road_Herb = 0.88;            // Herb 2008, (dimensionless) //0.88

  absortivity_soil_Best = 0.85;            // Herb 2008, (dimensionless)
  absortivity_road_Best = 0.88;            // Herb 2008, (dimensionless)

  emissivity_soil_Jansson = 0.97;          // Garrat 1994,  (dimensionless)
  emissivity_road_Jansson = 1.00;          // Jansson 2006, (dimensionless)

  emissivity_soil_Herb = 0.95;             // Herb 2008, (dimensionless)
  emissivity_road_Herb = 0.94;             // Herb 2008, (dimensionless) //0.94
  
  emissivity_soil_Best = 0.95;             // Herb 2008, (dimensionless)
  emissivity_road_Best = 0.94;             // Herb 2008, (dimensionless)
  emissivity_canopy_Best = 0.95;           // Herb 2008, (dimensionless)
  
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

}
/*************************************************************************************************************************************************************************/
/***************************************************                                                   *******************************************************************/
/***************************************************      Jansson et al. (2006) Surface Coefficients   *******************************************************************/
/***************************************************                                                   *******************************************************************/
/*************************************************************************************************************************************************************************/
double SurfaceCoefficients::get_absortivity_Jansson (const std::string surface)
{
  double coefficient = 0;
  if (surface == "road")
    coefficient = absortivity_road_Jansson;
  else if (surface == "soil")
    coefficient = absortivity_soil_Jansson;
  else
    {
      std::cout << "Error: not implemented"             << std::endl
		<< "error in get_absortivity_Jansson"   << std::endl;
      throw 3;
    }
  return coefficient;
}


double SurfaceCoefficients::get_infrared_coefficient_Jansson (const std::string surface_type,
							      const double air_temperature,              // (C)
							      const double relative_humidity,            // (%)
							      const double previous_air_temperature,     // (C)
							      const double previous_temperature_surface) // (C)
{
  double vapor_partial_pressure_air = ((relative_humidity/100.) * 611. *
				       exp( (water_latent_heat_of_vaporization*molar_mass_water/gas_constant) * // (Pa)
					    ( (1./273.15) - (1./(air_temperature+273.15)) ) ));
  
  //double sky_emissivity = a_Jansson + b_Jansson*pow((vapor_partial_pressure_air/100.),0.5); // Brunt 1932 (vapor pressure must be in hPa)

  double coefficient_cloud_cover = 0.59;
  double sky_emissivity = (coefficient_cloud_cover + 
			   0.67*(1-coefficient_cloud_cover)*pow(vapor_partial_pressure_air/100.,0.08)); /* Edinger 1956? Apparently the vapor pressure 
													   must be in hPa (requires confirmation, I
													   haven't been able to locate the reference.
													   It is cited in Herb et al 2008.
													*/
  double sky_temperature = pow(sky_emissivity,0.25) * (previous_air_temperature + 273.15);
  
  double surface_emissivity = 0.;
  if (surface_type == "road")
    surface_emissivity = emissivity_road_Jansson;
  else if (surface_type == "soil")
    surface_emissivity = emissivity_soil_Jansson;
  else
    {
      std::cout << "Error: not implemented"                    << std::endl
		<< "error in get_infrared_coefficient_Jansson" << std::endl;
      throw 3;
    }
  
  // std::cout << sky_temperature << std::endl;
  // std::cout << 	pow(0.25*(pow(sky_temperature,2) + pow((previous_temperature_surface + 273.15),2)) *
  // 		    (    sky_temperature    +    ( previous_temperature_surface + 273.15)    ),(1./3.)) << std::endl;
  
  return (steffan_boltzmann*surface_emissivity *
	  (pow(sky_temperature,2) + pow((previous_temperature_surface + 273.15),2)) *
	  (    sky_temperature    +    ( previous_temperature_surface + 273.15)    )); // (W/m2K)
  
  // return (4. * steffan_boltzmann * pow(0.5 * ((sky_temperature) + (previous_temperature_surface + 273.15)),3));
}

double SurfaceCoefficients::get_infrared_inbound_coefficient_Jansson (const std::string surface_type,
								      const double relative_humidity, // (%)
								      const double air_temperature)   // (C)
{
  double vapor_partial_pressure_air = ((relative_humidity/100.) * 611. *
				       exp( (water_latent_heat_of_vaporization*molar_mass_water/gas_constant) * // (Pa)
					    ( (1./273.15) - (1./(air_temperature+273.15)) ) ));
  //  double sky_emissivity = a_Jansson + b_Jansson*pow(vapor_partial_pressure_air/100.,0.5); // Brunt 1932 (vapor pressure must be in hPa)

  double coefficient_cloud_cover = 0.59;
  double sky_emissivity = (coefficient_cloud_cover + 
			   0.67*(1-coefficient_cloud_cover)*pow(vapor_partial_pressure_air/100.,0.08)); /* Edinger 1956? Apparently the vapor pressure 
													   must be in hPa (requires confirmation, I
													   haven't been able to locate the reference.
													   It is cited in Herb et al 2008.
													*/
  double surface_emissivity = 0.;
  if (surface_type == "road")
    surface_emissivity = emissivity_road_Jansson;
  else if (surface_type == "soil")
    surface_emissivity = emissivity_soil_Jansson;
  else
    {
      std::cout << "Error: not implemented"                            << std::endl
		<< "error in get_infrared_inbound_coefficient_Jansson" << std::endl;
      throw 3;
    }
  return (steffan_boltzmann * sky_emissivity * surface_emissivity); // (W/m2K4)
}

double SurfaceCoefficients::get_infrared_outbound_coefficient_Jansson (const std::string surface_type)
{
  double surface_emissivity;
  
  if (surface_type == "road")
    surface_emissivity = emissivity_road_Jansson;
  else if (surface_type == "soil")
    surface_emissivity = emissivity_soil_Jansson;
  else
    {
      std::cout << "Error: not implemented"                             << std::endl
		<< "error in get_infrared_outbound_coefficient_Jansson" << std::endl;
      throw 3;
    }
  return (steffan_boltzmann*surface_emissivity); // (W/m2K4)
}

double SurfaceCoefficients::get_convective_coefficient_Jansson(const std::string surface_type,
							       const double wind_speed)
{
  double resistance = 0;
  if (surface_type == "road")
    resistance = ((1./(pow(von_karman,2)*wind_speed)) *
		  log10(z_ref_wind/z_road_momentum_roughness) *
		  log10(z_ref_temperature/z_road_heat_roughness));
  else if (surface_type == "soil")
    resistance = ((1./(pow(von_karman,2)*wind_speed)) *
		  log10(z_ref_wind/z_soil_momentum_roughness) *
		  log10(z_ref_temperature/z_soil_heat_roughness));
  else
    {
      std::cout << "Error: not implemented"              << std::endl
		<< "error in get_convective_coefficient" << std::endl
		<< "surface type unknown"                << std::endl
		<< "options: 'road' 'soil'"              << std::endl;
      throw 3;
      return 1;
    }
  return (air_density*air_specific_heat_capacity/resistance);
}

double SurfaceCoefficients::get_evaporative_flux_Jansson (const std::string surface_type,
							  const double air_temperature  ,
							  const double wind_speed       ,
							  const double relative_humidity,
							  const double surface_temperature)
{
  if (surface_type == "soil")
    {   
      /*
	The equations used for (Jansson et al 2006) for evaporation from a bare soil assume a saturated soil.
	This is the same assumption as in the equations used by (Herb et al 2008). This equation does not
	include a term to take into account the contribution of natural convection. For this reason, for low
	wind speeds (more usual that one would expect according to the experimental data) the convective
	and evaporative heat fluxes are negligible, and this in turn bring the soil to very high temperatures.
	However, the advantage of this formulation is that it has been modified by (Philip et al 1957) to be
	applied to soils by including a 'soil moisture availability factor'
      */
      double vapor_pressure_saturation_surface = 611. * exp((water_latent_heat_of_vaporization*molar_mass_water/gas_constant) * // (Pa)
							    ((1./273.15) - (1./(surface_temperature + 273.15)) ) ); 
      /*
	from thermodynamics laws Philip (1957) dereive an expression for the relative humidity of air in 
	equilibrium with the water in the soil pore. A prerequisite for the equation to be valid is that 
	the air close to the pore water surface is in equilibrium with the pore water. However, in a drying
	soil, where vapor is continously transported to the atmosphere, the vapor pressure of air adjacent
	to the evaporating water surface will not be in equilibrium with the liquid water in the pore, i.e.
	the calculated relative humidity won't be representative at the soil surface.
      */
      double vapor_pressure_surface            = (vapor_pressure_saturation_surface * 
						  exp (-1. * surface_water_tension * molar_mass_water * gravity_constant/(gas_constant * (surface_temperature + 273.15)) )); // (Pa)
      /*
	An obvious drawback of this formulation is that we need to know the surface water tension at the
	surface (or near the surface) in order ot use it and unless we have experimental measurements 
	or we are dealing with a fully coupled heat-moisture diffusion model, it might prove to be 
	difficult to give a reasonable number and even more difficult to propose its variation in time.
      */
      double vapor_pressure_saturation_air = 611. * exp((water_latent_heat_of_vaporization*molar_mass_water/gas_constant) * // (Pa)
							((1./273.15) - (1./(air_temperature + 273.15)) ) );      
      double vapor_pressure_air            = vapor_pressure_saturation_air * relative_humidity / 100. ; // (Pa)
      
      double psychrometric_constant = (air_specific_heat_capacity * atmospheric_pressure /
				       ((molar_mass_water/molar_mass_air) * water_latent_heat_of_vaporization)); // (Pa/K)
      
      double aerodynamic_resistance = ((1./(pow(von_karman,2)*wind_speed)) *
				       log10(z_ref_wind/z_soil_momentum_roughness) *
				       log10(z_ref_temperature/z_soil_heat_roughness)); // (s/m)

      // std::cout << air_density * air_specific_heat_capacity / (psychrometric_constant * aerodynamic_resistance) << std::endl;
      // std::cout << vapor_pressure_air << std::endl;
      // std::cout << vapor_pressure_surface << std::endl;
      // std::cout << vapor_pressure_saturation_surface << std::endl;
      // std::cout << exp (-1. * surface_water_tension * molar_mass_water * gravity_constant/(gas_constant * (surface_temperature + 273.15)) ) << std::endl;

      // std::cout << 611. * exp((water_latent_heat_of_vaporization*molar_mass_water/gas_constant) * // (Pa)
      // 			      ((1./273.15) ) ) << std::endl;
      
      return (air_density * air_specific_heat_capacity * (vapor_pressure_air - vapor_pressure_surface) / (psychrometric_constant * aerodynamic_resistance) ); // (W/m2)
    }
  else if (surface_type == "road")
    {
      /*
	Assumption: No evaporation from the pavement surface. This means that we are considering it as
	a nonporous surface. In principle, if there is rainfall, the surface will get wet and some
	evaporation will occur. So, we are also assuming that the rainfall events are scarce. This
	assumption, however, might be wrong as, in my experience, heavy rainfalls in the uk are uncommon,
	drizle is a lot more usual.
       */
      return 0.; // (W/m2)
    }
  else
    {
      std::cout << "Error: not implemented"   << std::endl
		<< "error in get_evaporative" << std::endl;
      throw 3;
      return 1;
    }
  
}

//*****************************************************************************************************************



double SurfaceCoefficients::evaporative_inbound_coefficient (const std::string surface_type = "soil",
							     const double wind_speed        = 1.1413  /* average_wind_speed (m/s) based on TRL's meteorolgical data       */,
							     const double relative_humidity = 80.66   /* average_relative_humidity (%) based on TRL's meteorological data */,
							     const double air_temperature   = 10.24   /* average_air_temperature (K) based on TRL's meteorological data   */)
{
  if (surface_type == "soil")
    {
      double vapor_pressure_saturation_air = exp(20.386 - 5132/(air_temperature+273.15)) * (0.133322368); // (kPa)
      double vapor_pressure_air            = vapor_pressure_saturation_air * relative_humidity / 100. ; // (kPa)
      
      //double atmospheric_pressure        = 101.325;  // (kPa)
      double latent_heat_of_vaporization = 2.45E6; // (J/kg) 
      double ratio_water_molecular_weight_to_dry_air = 0.622; // (dimensionless)
      
      double psychrometric_constant = air_specific_heat_capacity * (atmospheric_pressure/1000.) / 
	(ratio_water_molecular_weight_to_dry_air * latent_heat_of_vaporization); // (kPa/K)
      
      double aerodynamic_resistance = (1./(pow(von_karman,2)*wind_speed)) *
	log10(z_ref_wind/z_soil_momentum_roughness) *
	log10(z_ref_temperature/z_soil_heat_roughness); // (s/m)
      
      // std::cout << wind_speed << "\t"
      // 		<< relative_humidity << "\t"
      // 		<< air_temperature << "\t"
      // 		<< vapor_pressure_saturation_air << "\t"
      // 		<< vapor_pressure_air << "\t"
      // 		<< psychrometric_constant << "\t"
      // 		<< aerodynamic_resistance << "\t" 
      // 		<< air_density * air_specific_heat_capacity * vapor_pressure_air / (psychrometric_constant * aerodynamic_resistance) 
      // 		<< std::endl;
      
      return (air_density * air_specific_heat_capacity * vapor_pressure_air / (psychrometric_constant * aerodynamic_resistance) ); // (W/m2)
    }
  else if (surface_type == "road")
    {
      return 0.; // (W/m2)
    }
  else
    {
      std::cout << "Error: not implemented"   << std::endl
		<< "error in get_evaporative" << std::endl;
      throw 3;
      return 1;
    }
  
}

double SurfaceCoefficients::evaporative_outbound_coefficient (const std::string surface_type = "soil",
							      const double wind_speed        = 1.1413  /* average_wind_speed (m/s) based on TRL's meteorolgical data       */,
							      const double surface_temperature = -100000.)
{
  if (surface_type == "soil")
    {
      double vapor_pressure_saturation_surface = exp(20.386 - 5132/(surface_temperature+273.15)) * (0.133322368); // (kPa)
      double vapor_pressure_surface            = (vapor_pressure_saturation_surface 
						  * exp (-1.*surface_water_tension*molar_mass_water*gravity_constant/(gas_constant*(surface_temperature+273.15)) )); // (kPa)
      
      //double atmospheric_pressure        = 101.325;  // (kPa)
      double latent_heat_of_vaporization = 2.45E6; // (J/kg) 
      double ratio_water_molecular_weight_to_dry_air = 0.622; // (dimensionless)
      
      double psychrometric_constant = air_specific_heat_capacity * (atmospheric_pressure/1000.) / 
	(ratio_water_molecular_weight_to_dry_air * latent_heat_of_vaporization); // (kPa/K)
      
      double aerodynamic_resistance = (1./(pow(von_karman,2)*wind_speed)) *
	log10(z_ref_wind/z_soil_momentum_roughness) *
	log10(z_ref_temperature/z_soil_heat_roughness); // (s/m)
      
      // std::cout << wind_speed << "\t"
      // 		<< surface_temperature << "\t"
      // 		<< vapor_pressure_saturation_surface << "\t"
      // 		<< vapor_pressure_surface << "\t"
      // 		<< psychrometric_constant << "\t"
      // 		<< aerodynamic_resistance << "\t" 
      // 		<< air_density * air_specific_heat_capacity * vapor_pressure_surface / (psychrometric_constant * aerodynamic_resistance) 
      // 		<< std::endl;


      return (air_density * air_specific_heat_capacity * vapor_pressure_surface / (psychrometric_constant * aerodynamic_resistance) ); // (W/m2)
    }
  else if (surface_type == "road")
    {
      return 0.; // (W/m2)
    }
  else
    {
      std::cout << "Error: not implemented"   << std::endl
		<< "error in get_evaporative" << std::endl;
      throw 3;
      return 1;
    }
  
}

/*************************************************************************************************************************************************************************/
/***************************************************                                                   *******************************************************************/
/***************************************************      Herb et al. (2008) Surface Coefficients      *******************************************************************/
/***************************************************                                                   *******************************************************************/
/*************************************************************************************************************************************************************************/
double SurfaceCoefficients::get_absortivity_Herb (const std::string surface)
{
  double coefficient = 0;
  if (surface == "road")
    coefficient = absortivity_road_Herb;
  else if (surface == "soil")
    coefficient = absortivity_soil_Herb;
  else
    {
      std::cout << "Error: not implemented"          << std::endl
		<< "error in get_absortivity_Herb"   << std::endl;
      throw 3;
    }
  return coefficient;
}

double SurfaceCoefficients::get_infrared_inbound_coefficient_Herb(const std::string surface_type,
								  const double relative_humidity, // (%)
								  const double air_temperature)   // (C)
{
  double coefficient_cloud_cover = .0;
  double vapor_partial_pressure_air = ((relative_humidity/100.) * 
				       611. * exp( (water_latent_heat_of_vaporization*molar_mass_water/gas_constant)* // (Pa)
						   ( (1./273.15) - (1./(air_temperature+273.15)) ) ) ); 
  
  double sky_emissivity = (coefficient_cloud_cover + 
			   0.67*(1-coefficient_cloud_cover)*pow(vapor_partial_pressure_air/100.,0.08)); /* Edinger 1956? Apparently the vapor pressure 
													   must be in hPa (requires confirmation, I
													   haven't been able to locate the reference.
													   It is cited in Herb et al 2008.
													*/
  
  double surface_emissivity = 0.;
  if (surface_type == "road")
    surface_emissivity = emissivity_road_Herb;
  else if (surface_type == "soil")
    surface_emissivity = emissivity_soil_Herb;
  else
    {
      std::cout << "Error: not implemented"                   << std::endl
		<< "error in get_infrared_coefficient_Herb"   << std::endl;
      throw 3;
    }
  return ( steffan_boltzmann * sky_emissivity * surface_emissivity); // (W/m2K4)
}

double SurfaceCoefficients::get_infrared_outbound_coefficient_Herb(const std::string surface_type)
{
  double surface_emissivity = 0.;
  if (surface_type == "road")
    surface_emissivity = emissivity_road_Herb;
  else if (surface_type == "soil")
    surface_emissivity = emissivity_soil_Herb;
  else
    {
      std::cout << "Error: not implemented"                          << std::endl
		<< "error in get_infrared_outbound_coefficient_Herb" << std::endl;
      throw 3;
    }
  return ( steffan_boltzmann * surface_emissivity); // (W/m2K4)
}

double SurfaceCoefficients::get_infrared_coefficient_Herb (const std::string surface_type,
							   const double air_temperature,              // (C)
							   const double relative_humidity,            // (%)
							   const double previous_air_temperature,     // (C)
							   const double previous_temperature_surface) // (C)
{
  double coefficient_cloud_cover = .0;
  double vapor_partial_pressure_air = ((relative_humidity/100.) * 611. *
				       exp( (water_latent_heat_of_vaporization*molar_mass_water/gas_constant) * // (Pa)
					    ( (1./273.15) - (1./(air_temperature+273.15)) ) ) );
  
  double sky_emissivity  = (coefficient_cloud_cover + 
			    0.67 * (1-coefficient_cloud_cover) * pow(vapor_partial_pressure_air/100.,0.08)); /* Edinger 1956? Apparently the vapor pressure 
														must be in hPa (requires confirmation, I
														haven't been able to locate the reference.
														It is cited in Herb et al 2008.
													     */  
  double sky_temperature = pow(sky_emissivity,0.25) * (previous_air_temperature + 273.15); 
  
  double surface_emissivity = 0.;
  if (surface_type == "road")
    surface_emissivity = emissivity_road_Herb;
  else if (surface_type == "soil")
    surface_emissivity = emissivity_soil_Herb;
  else
    {
      std::cout << "Error: not implemented"                 << std::endl
		<< "error in get_infrared_coefficient_Herb" << std::endl;
      throw 3;
    }
  return (steffan_boltzmann * surface_emissivity *
	  (pow(sky_temperature,2) + pow((previous_temperature_surface + 273.15),2)) *
	  (    sky_temperature    +    ( previous_temperature_surface + 273.15)    )); // (W/m2K)
}


double SurfaceCoefficients::get_convective_coefficient_Herb(const std::string /*surface_type*/,
							    const double wind_speed,          // (m/s)
							    const double relative_humidity,   // (%)
							    const double air_temperature,     // (C)
							    const double temperature_surface) // (C)
{
  // double coefficient_forced_convection  = 0.0015; // (dimensionless) nominal value according to Herb et al 2008
  // double coefficient_natural_convection = 0.0015; // (m/s  1/K^0.33) nominal value according to Herb et al 2008
  // double coefficient_sheltering         = 1.0000; // (dimensionless)
  
  double saturated_vapor_partial_pressure_surface = 611. * exp( (water_latent_heat_of_vaporization*molar_mass_water/gas_constant)*  // (Pa)
								( (1./273.15) - (1./(temperature_surface+273.15)) ) );  
  
  double saturated_vapor_partial_pressure_air     = 611. * exp( (water_latent_heat_of_vaporization*molar_mass_water/gas_constant)* // (Pa)
								( (1./273.15) - (1./(air_temperature+273.15)) ) );
  
  double mixing_ratio_surface = (molar_mass_water/molar_mass_air) * saturated_vapor_partial_pressure_surface/atmospheric_pressure; // (Pa/Pa)
  double mixing_ratio_air     = ((relative_humidity/100.) *
				 (molar_mass_water/molar_mass_air) * saturated_vapor_partial_pressure_air/atmospheric_pressure);     // (Pa/Pa)
  
  double virtual_temperature_surface    = (1. + 0.6*mixing_ratio_surface) * temperature_surface + 273.15; // (K)
  double virtual_air_temperature        = (1. + 0.6*mixing_ratio_air    ) * air_temperature     + 273.15; // (K)
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
  double virtual_temperature_difference = 0.;
  if (virtual_temperature_surface > virtual_air_temperature)
    virtual_temperature_difference = virtual_temperature_surface - virtual_air_temperature;
  else
    virtual_temperature_difference = 0.;
  /* 
     Specific humidity (q) is the number of grams of water vapor per unit mass
     of air plus water vapor. In applications, q is numerically close enough to
     the mixing ratio that one seldom needs to distinguish between them.
     [North and Erukhimova 2009, "Atmospheric Thermodynamics Elementary Physics
     and Chemistry"
  */
  return (air_density*air_specific_heat_capacity*(coefficient_forced_convection*coefficient_sheltering*wind_speed + 
						  coefficient_natural_convection*pow(virtual_temperature_difference,1./3.))); // (W/m^2K)
}

double SurfaceCoefficients::get_evaporative_flux_Herb(const std::string /*surface_type*/,
						      const double air_temperature,     // (C)
						      const double wind_speed,          // (m/s)
						      const double relative_humidity,   // (%)
						      const double temperature_surface) // (C)
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
    capillary forces, and the vapor pressure of the air in equilibrium with the
    pore water will therefore be lower than in the air close to a free water surface"
    [Alvenas and Janssen 1997, "Model for evaporation, moisture and temperature"]
    ...
    This means that the assumption of a soil saturated at all times should be revised
    if we are in weather conditions that are likely to dry the soil.
    For the road case the situation is different. We assume the road surface is dry
    except when a rain event occurs. This should be decided by the function calling
    this member.
  */
  // double coefficient_forced_convection  = 0.0015; // (dimensionless) nominal value according to Herb et al 2008
  // double coefficient_natural_convection = 0.0015; // (m/s  1/K^0.33) nominal value according to Herb et al 2008
  // double coefficient_sheltering         = 1.0000; // (dimensionless)
  /* 
     Clausius-Clapeyron equation
     Relates equilibrium or saturation vapor pressure and temperature to the
     latent heat of the phase change, without requiring specific volume data
  */
  double saturated_vapor_partial_pressure_surface = 611. * exp( (water_latent_heat_of_vaporization*molar_mass_water/gas_constant)*  // (Pa)
								( (1./273.15) - (1./(temperature_surface+273.15)) ) );  
  
  double saturated_vapor_partial_pressure_air     = 611. * exp( (water_latent_heat_of_vaporization*molar_mass_water/gas_constant)* // (Pa)
								( (1./273.15) - (1./(air_temperature+273.15)) ) );
 
  double mixing_ratio_surface = (molar_mass_water/molar_mass_air) * saturated_vapor_partial_pressure_surface/atmospheric_pressure; // (Pa/Pa)
  double mixing_ratio_air     = ((relative_humidity/100.) *
				 (molar_mass_water/molar_mass_air) * saturated_vapor_partial_pressure_air/atmospheric_pressure);     // (Pa/Pa)
  
  double virtual_temperature_surface    = (1. + 0.6*mixing_ratio_surface) * temperature_surface + 273.15; // (K)
  double virtual_air_temperature        = (1. + 0.6*mixing_ratio_air    ) * air_temperature     + 273.15; // (K)
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
  double virtual_temperature_difference = 0.;
  if (virtual_temperature_surface > virtual_air_temperature)
    virtual_temperature_difference = virtual_temperature_surface - virtual_air_temperature;
  else
    virtual_temperature_difference = 0.;
  /* 
     Specific humidity (q) is the number of grams of water vapor per unit mass
     of air plus water vapor. In applications, q is numerically close enough to
     the mixing ratio that one seldom needs to distinguish between them.
     [North and Erukhimova 2009, "Atmospheric Thermodynamics Elementary Physics
     and Chemistry"
  */
  // std::cout << air_density * water_latent_heat_of_vaporization << std::endl 
  // 	    << coefficient_forced_convection*wind_speed + coefficient_natural_convection*pow(virtual_temperature_difference,1./3.) << "\t"
  // 	    << (mixing_ratio_surface - mixing_ratio_air) << std::endl
  // 	    << saturated_vapor_partial_pressure_surface << "\t" << mixing_ratio_surface << "\t" << temperature_surface + 273.15 << "\t" 
  // 	    << virtual_temperature_surface << std::endl
  // 	    << saturated_vapor_partial_pressure_air     << "\t" << mixing_ratio_air     << "\t" << air_temperature + 273.15     << "\t" 
  // 	    << virtual_air_temperature << std::endl
  // 	    << std::endl;
  // if (surface_type=="road")
  // std::cout << "Atmospheric pressure:" << atmospheric_pressure << std::endl
  // 	    << "Air Density: " << air_density                       << std::endl
  // 	    << "Lv         : " << water_latent_heat_of_vaporization << std::endl
  // 	    << "Cfc        : " << coefficient_forced_convection     << std::endl
  // 	    << "Cnc        : " << coefficient_natural_convection    << std::endl
  // 	    << "Csh        : " << coefficient_sheltering            << std::endl
  // 	    << "Wind Speed         : " << wind_speed                << std::endl
  // 	    << "Air temperature    : " << air_temperature           << std::endl
  // 	    << "Surface temperature: " << temperature_surface       << std::endl
  // 	    << "Saturated vapor pressure surface: " << saturated_vapor_partial_pressure_surface << std::endl
  // 	    << "Saturated vapor pressure air    : " << saturated_vapor_partial_pressure_air     << std::endl
  // 	    << "Mixing ratio air    : " << mixing_ratio_air     << std::endl
  // 	    << "Mixing ratio surface: " << mixing_ratio_surface << std::endl
  // 	    << "Virtual temperature air       : " << virtual_air_temperature        << std::endl
  // 	    << "Virtual temperature surface   : " << virtual_temperature_surface    << std::endl
  // 	    << "Virtual temperature difference: " << virtual_temperature_difference << std::endl
  // 	    << "Evaporative Heat Flux: " << (air_density * water_latent_heat_of_vaporization * (coefficient_forced_convection  * coefficient_sheltering * wind_speed + 
  // 												coefficient_natural_convection * pow(virtual_temperature_difference,1./3.))
  // 					     * (mixing_ratio_air - mixing_ratio_surface) ) << std::endl
  // 	    << std::endl;

  return (air_density * water_latent_heat_of_vaporization * (coefficient_forced_convection  * coefficient_sheltering * wind_speed + 
							     coefficient_natural_convection * pow(virtual_temperature_difference,1./3.))
	  * (mixing_ratio_air - mixing_ratio_surface) ); // (W/m^2)
}

/*************************************************************************************************************************************************************************/
/***************************************************                                                   *******************************************************************/
/***************************************************      Best et al. (2008) Surface Coefficients      *******************************************************************/
/***************************************************                                                   *******************************************************************/
/*************************************************************************************************************************************************************************/
double SurfaceCoefficients::get_absortivity_Best (const std::string surface)
{
  double coefficient = 0;
  if (surface == "road")
    coefficient = absortivity_road_Best;
  else if (surface == "soil")
    coefficient = absortivity_soil_Best;
  else
    {
      std::cout << "Error: not implemented"          << std::endl
		<< "error in get_absortivity_Best"   << std::endl;
      throw 3;
    }
  return coefficient;
}

double SurfaceCoefficients::get_infrared_inbound_coefficient_Best(const std::string surface_type,
								  const double relative_humidity,
								  const double air_temperature)
{
  double coefficient_cloud_cover = .0;
  double vapor_partial_pressure_air = ((relative_humidity/100.) * 
				       611. * exp( (water_latent_heat_of_vaporization*molar_mass_water/gas_constant)* // (Pa)
						   ( (1./273.15) - (1./(air_temperature+273.15)) ) ) ); 
  double sky_emissivity = (coefficient_cloud_cover + 
			   0.67*(1-coefficient_cloud_cover)*pow(vapor_partial_pressure_air/100.,0.08)); /* Edinger 1956? Apparently the vapor pressure 
													   must be in hPa (requires confirmation, I
													   haven't been able to locate the reference.
													   It is cited in Herb et al 2008.
													*/
  double surface_emissivity = 0.;
  if (surface_type == "road")
    surface_emissivity = emissivity_road_Best;
  else if (surface_type == "soil")
    surface_emissivity = emissivity_soil_Best;
  else if (surface_type == "canopy")
    surface_emissivity = emissivity_canopy_Best;
  else
    {
      std::cout << "Error: not implemented"                   << std::endl
		<< "error in get_infrared_coefficient_Best"   << std::endl;
      throw 3;
    }
  return ( steffan_boltzmann * sky_emissivity * surface_emissivity); // (W/m2K4)
}

double SurfaceCoefficients::get_infrared_outbound_coefficient_Best(const std::string surface_type)
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
      std::cout << "Error: not implemented"                          << std::endl
		<< "error in get_infrared_outbound_coefficient_Best" << std::endl;
      throw 3;
    }
  return ( steffan_boltzmann * surface_emissivity); // (W/m2K4)
}

double SurfaceCoefficients::get_convective_coefficient_Best(const std::string /*surface_type*/,
							    const double wind_speed,          // (m/s)
							    const double relative_humidity,   // (%)
							    const double air_temperature,     // (C)
							    const double temperature_surface) // (C)
{  
  double saturated_vapor_partial_pressure_surface = 611. * exp( (water_latent_heat_of_vaporization*molar_mass_water/gas_constant)*  // (Pa)
								( (1./273.15) - (1./(temperature_surface+273.15)) ) );  
  
  double saturated_vapor_partial_pressure_air     = 611. * exp( (water_latent_heat_of_vaporization*molar_mass_water/gas_constant)* // (Pa)
								( (1./273.15) - (1./(air_temperature+273.15)) ) );
  
  double mixing_ratio_surface = (molar_mass_water/molar_mass_air) * saturated_vapor_partial_pressure_surface/atmospheric_pressure; // (Pa/Pa)
  double mixing_ratio_air     = ((relative_humidity/100.) *
				 (molar_mass_water/molar_mass_air) * saturated_vapor_partial_pressure_air/atmospheric_pressure);     // (Pa/Pa)
  
  double virtual_temperature_surface    = (1. + 0.6*mixing_ratio_surface) * temperature_surface + 273.15; // (K)
  double virtual_air_temperature        = (1. + 0.6*mixing_ratio_air    ) * air_temperature     + 273.15; // (K)
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
  double virtual_temperature_difference = 0.;
  if (virtual_temperature_surface > virtual_air_temperature)
    virtual_temperature_difference = virtual_temperature_surface - virtual_air_temperature;
  else
    virtual_temperature_difference = 0.;
  /* 
     Specific humidity (q) is the number of grams of water vapor per unit mass
     of air plus water vapor. In applications, q is numerically close enough to
     the mixing ratio that one seldom needs to distinguish between them.
     [North and Erukhimova 2009, "Atmospheric Thermodynamics Elementary Physics
     and Chemistry"
  */
  return (air_density*air_specific_heat_capacity*(coefficient_forced_convection*coefficient_sheltering*wind_speed + 
						  coefficient_natural_convection*pow(virtual_temperature_difference,1./3.))); // (W/m^2K)
}

double SurfaceCoefficients::get_evaporative_flux_Best(const std::string /*surface_type*/,
						      const double air_temperature,     // (C)
						      const double wind_speed,          // (m/s)
						      const double relative_humidity,   // (%)
						      const double new_temperature_surface) // (C)
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
    capillary forces, and the vapor pressure of the air in equilibrium with the
    pore water will therefore be lower than in the air close to a free water surface"
    [Alvenas and Janssen 1997, "Model for evaporation, moisture and temperature"]
  */
  /* 
     Clausius-Clapeyron equation
     Relates equilibrium or saturation vapor pressure and temperature to the
     latent heat of the phase change, without requiring specific volume data
  */


  double saturated_vapor_partial_pressure_surface = 611. * exp( (water_latent_heat_of_vaporization*molar_mass_water/gas_constant)*  // (Pa)
  								( (1./273.15) - (1./(new_temperature_surface+273.15)) ) );  
  
  double saturated_vapor_partial_pressure_air     = 611. * exp( (water_latent_heat_of_vaporization*molar_mass_water/gas_constant)* // (Pa)
  								( (1./273.15) - (1./(air_temperature+273.15)) ) );
  
  
  // double dummy1 = (water_latent_heat_of_vaporization*molar_mass_water/gas_constant);
  // double dummy2 = 611. * exp(dummy1/273.15);

  // // double saturated_vapor_partial_pressure_surface = dummy2 * exp(-dummy1/(temperature_surface+273.15));
  // double x = -dummy1/(new_temperature_surface+273.15);
  // double y = -dummy1/(new_temperature_surface+273.15);
  // double saturated_vapor_partial_pressure_surface = dummy2 * (1. 
  // 							      + y                    /*  <-------- !! */
  // 							      + pow(x,2)/(1.*2.)
  // 							      + pow(x,3)/(1.*2.*3.)
  // 							      + pow(x,4)/(1.*2.*3.*4.)
  // 							      + pow(x,5)/(1.*2.*3.*4.*5.)
  // 							      + pow(x,6)/(1.*2.*3.*4.*5.*6.)
  // 							      + pow(x,7)/(1.*2.*3.*4.*5.*6.*7.)
  // 							      + pow(x,8)/(1.*2.*3.*4.*5.*6.*7.*8.)
  // 							      + pow(x,9)/(1.*2.*3.*4.*5.*6.*7.*8.*9.)
  // 							      + pow(x,10)/(1.*2.*3.*4.*5.*6.*7.*8.*9.*10.)
  // 							      + pow(x,11)/(1.*2.*3.*4.*5.*6.*7.*8.*9.*10.*11.));
  
  // //  double saturated_vapor_partial_pressure_air     = dummy2 * exp(-dummy1/(air_temperature+273.15));
  // x = -dummy1/(air_temperature+273.15);
  // double saturated_vapor_partial_pressure_air     = dummy2 * (1. 
  // 							      + x 
  // 							      + pow(x,2)/(1.*2.)
  // 							      + pow(x,3)/(1.*2.*3.)
  // 							      + pow(x,4)/(1.*2.*3.*4.)
  // 							      + pow(x,5)/(1.*2.*3.*4.*5.)
  // 							      + pow(x,6)/(1.*2.*3.*4.*5.*6.)
  // 							      + pow(x,7)/(1.*2.*3.*4.*5.*6.*7.)
  // 							      + pow(x,8)/(1.*2.*3.*4.*5.*6.*7.*8.)
  // 							      + pow(x,9)/(1.*2.*3.*4.*5.*6.*7.*8.*9.)
  // 							      + pow(x,10)/(1.*2.*3.*4.*5.*6.*7.*8.*9.*10.)
  // 							      + pow(x,11)/(1.*2.*3.*4.*5.*6.*7.*8.*9.*10.*11.));
  

  double mixing_ratio_surface = (molar_mass_water/molar_mass_air) * saturated_vapor_partial_pressure_surface/atmospheric_pressure; // (Pa/Pa)
  double mixing_ratio_air     = ((relative_humidity/100.) *
				 (molar_mass_water/molar_mass_air) * saturated_vapor_partial_pressure_air/atmospheric_pressure);     // (Pa/Pa)
  
  double virtual_temperature_surface = (1. + 0.6*mixing_ratio_surface) * new_temperature_surface + 273.15; // (K)
  double virtual_air_temperature     = (1. + 0.6*mixing_ratio_air    ) * air_temperature     + 273.15; // (K)
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
  double virtual_temperature_difference = 0.;
  if (virtual_temperature_surface > virtual_air_temperature)
    virtual_temperature_difference = virtual_temperature_surface - virtual_air_temperature;
  else
    virtual_temperature_difference = 0.;
  /* 
     Specific humidity (q) is the number of grams of water vapor per unit mass
     of air plus water vapor. In applications, q is numerically close enough to
     the mixing ratio that one seldom needs to distinguish between them.
     [North and Erukhimova 2009, "Atmospheric Thermodynamics Elementary Physics
     and Chemistry"
  */
  // std::cout << air_density * water_latent_heat_of_vaporization << "\t" 
  // 	    << coefficient_forced_convection*wind_speed + coefficient_natural_convection*pow(virtual_temperature_difference,1./3.) << "\t"
  // 	    << (mixing_ratio_surface - mixing_ratio_air) << std::endl
  // 	    << saturated_vapor_partial_pressure_surface << "\t" << mixing_ratio_surface << "\t" << temperature_surface + 273.15 << "\t" << virtual_temperature_surface << std::endl
  //  	    << saturated_vapor_partial_pressure_air     << "\t" << mixing_ratio_air     << "\t" << air_temperature + 273.15     << "\t" << virtual_air_temperature << std::endl
  // 	    << std::endl;
  
  return (air_density * water_latent_heat_of_vaporization * (coefficient_forced_convection  * coefficient_sheltering * wind_speed + 
							     coefficient_natural_convection * pow(virtual_temperature_difference,1./3.))
	  * (mixing_ratio_air - mixing_ratio_surface) ); // (W/m^2)
}

double SurfaceCoefficients::get_convective_coefficient_canopy_Best(const std::string /*surface_type*/,
								   const double wind_speed)        // (m/s)
{ 
  return (0.01*(wind_speed+0.3)*air_density*air_specific_heat_capacity); // (W/m^2K)
}

double SurfaceCoefficients::get_evaporative_flux_canopy_Best(const std::string /*surface_type*/,
							     const double air_temperature,     // (C)
							     const double wind_speed,          // (m/s)
							     const double relative_humidity,   // (%)
							     const double solar_radiation,     // (W/m2)
							     const double temperature_surface, // (C)
							     const bool derivative) 
{
  double maximum_solar_radiation = 1000.;
  double moisture_content = 0.24;
  double wilting_point = 0.5*moisture_content;
  double stomata_resistance = 200.*(maximum_solar_radiation/(solar_radiation+0.3*maximum_solar_radiation) 
				   + pow(wilting_point/moisture_content,2));
  
  double saturated_vapor_partial_pressure_surface = 611. * exp( (water_latent_heat_of_vaporization*molar_mass_water/gas_constant)*  // (Pa)
								( (1./273.15) - (1./(temperature_surface+273.15)) ) );  
  
  double saturated_vapor_partial_pressure_air     = 611. * exp( (water_latent_heat_of_vaporization*molar_mass_water/gas_constant)* // (Pa)
								( (1./273.15) - (1./(air_temperature+273.15)) ) );
  
  double mixing_ratio_surface = (molar_mass_water/molar_mass_air) * saturated_vapor_partial_pressure_surface/atmospheric_pressure; // (Pa/Pa)
  double mixing_ratio_air     = ((relative_humidity/100.) *
				 (molar_mass_water/molar_mass_air) * saturated_vapor_partial_pressure_air/atmospheric_pressure);     // (Pa/Pa)
  
  if (derivative == false)
    return (air_density * water_latent_heat_of_vaporization 
	    * (0.01*(wind_speed+0.3)/(1.+stomata_resistance*0.01*(wind_speed+0.3)))
	    * (mixing_ratio_air - mixing_ratio_surface) ); // (W/m^2)
  else 
    return(-1. * air_density * water_latent_heat_of_vaporization 
	   * (0.01*(wind_speed+0.3)/(1.+stomata_resistance*0.01*(wind_speed+0.3)))
	   * mixing_ratio_surface * ((water_latent_heat_of_vaporization*molar_mass_water/gas_constant)/pow(temperature_surface+273.15,2)));
  
}

#endif
