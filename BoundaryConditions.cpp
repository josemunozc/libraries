#include "BoundaryConditions.h"
#include <iostream>
#include <math.h>
BoundaryConditions::
BoundaryConditions (const bool   analytic_,
		    const double time_,
		    double experimental_air_temperature_,
		    double experimental_solar_radiation_,
		    double experimental_wind_speed_,
		    double experimental_relative_humidity_,
		    double experimental_precipitation_,
		    double experimental_surface_temperature_,
		    double experimental_surface_pressure_,
		    const bool moisture_movement_)
:
  analytic            (analytic_),
  time                (time_),
  air_temperature     (experimental_air_temperature_),
  solar_radiation     (experimental_solar_radiation_),
  wind_speed          (experimental_wind_speed_),
  relative_humidity   (experimental_relative_humidity_),
  precipitation       (experimental_precipitation_),
  surface_temperature (experimental_surface_temperature_),
  surface_pressure    (experimental_surface_pressure_),
  moisture_movement   (moisture_movement_)
{
  pi    = 3.14159265359;
  phi   = 2.*pi/(24.*3600.);
  theta = 2.*pi/(24.*3600.*365.25);
    
  if (moisture_movement==false)
    surface_pressure=-75.2025;
  
  if (analytic == true)
    {
      if (time == -1000)
	{
	  std::cout << "Error in BoundaryConditions constructor."     << std::endl 
		    << "Undefined time. Time = " << time << std::endl;
	  throw 100;
	}
      
      AnalyticSolution analytic_solution;
      air_temperature   = analytic_solution.get_analytic_air_temperature  (  time  ); // (C)
      solar_radiation   = analytic_solution.get_analytic_solar_radiation  (  time  ); // (W/m2)
      wind_speed        = analytic_solution.get_analytic_wind_speed       (/*time*/); // (m/s)
      relative_humidity = analytic_solution.get_analytic_relative_humidity(/*time*/); // (%)
      precipitation     = 0.;
    }
  
  else if ((air_temperature     == -1000) || 
	   (solar_radiation     == -1000) ||
	   (wind_speed          == -1000) ||
	   (relative_humidity   == -1000) ||
	   (surface_temperature == -1000) ||
	   (surface_pressure    == -1000) ||
	   (precipitation       == -1000))
    {
      std::cout << "Error in met data."     << std::endl 
		<< "Air temperature     = " << air_temperature     << std::endl
		<< "Solar radiation     = " << solar_radiation     << std::endl
		<< "Wind speed          = " << wind_speed          << std::endl
		<< "Relative humidity   = " << relative_humidity   << std::endl
		<< "Surface temperature = " << surface_temperature << std::endl
		<< "Surface pressure    = " << surface_pressure    << std::endl;
      throw 100;
    }
}
//**************************************************************************************************
//**************************************************************************************************
double BoundaryConditions::
get_solar_flux(const std::string surface_type,
	       const std::string author,
	       const double canopy_density)
{    
  SurfaceCoefficients surface_coefficients(surface_pressure,
					   precipitation,
					   moisture_movement);
  double absortivity = 0.0; 
  if (author=="Herb")
    {
      absortivity
	=surface_coefficients
	.get_absortivity_Herb(surface_type);
      return(absortivity*solar_radiation); // (W/m2) 
    }
  else if (author=="Jansson")
    {
      absortivity
	=surface_coefficients
	.get_absortivity_Jansson(surface_type);
      return(absortivity*solar_radiation); // (W/m2) 
    }
  else if (author=="Best")
    {
      if ( (canopy_density < 0.) ||
	   (canopy_density > 1.))
	{
	  std::cout << "Error: canopy density out of range" << std::endl
		    << "error in get_solar_flux"      << std::endl;
	  throw 100;
	}
      absortivity
	=surface_coefficients
	.get_absortivity_Best(surface_type);
      return((1.-canopy_density)*absortivity*solar_radiation); // (W/m2) 
    }
  else
    {
      std::cout << "Error: author not implemented." << std::endl
		<< "Error in get_solar_flux."       << std::endl;
      throw 3;
    }
}

double BoundaryConditions::
get_infrared_flux(const std::string surface_type,
		  const std::string author,
		  const std::string direction,
		  const double old_surface_temperature, // (C)
		  const double canopy_temperature,      // (C)
		  const double canopy_density)
{  
  if ((old_surface_temperature<-20.) ||
      (old_surface_temperature> 80.))
    {
      std::cout << "Warning: out of range.\n"
		<< "Warning in get_infrared_flux.\n"
		<< "Author: "              << author << std::endl
		<< "Surface temperature: " << old_surface_temperature << std::endl
		<< "Canopy temperature : " << canopy_temperature << std::endl
		<< "Canopy density     : " << canopy_density << std::endl;
      //throw 100;
    }

  SurfaceCoefficients surface_coefficients(surface_pressure,
					   precipitation,
					   moisture_movement);
  double outbound_coefficient=0.;
  double inbound_flux=0.;
  double outbound_flux=0.;
 
  double sky_emissivity
    =surface_coefficients
    .get_sky_emissivity(author,
			air_temperature,
			relative_humidity);
  double sky_temperature
    =pow(sky_emissivity,0.25)*(air_temperature+273.15);
  
  if (author=="Herb")
    {
      outbound_coefficient
      	=surface_coefficients
      	.get_infrared_outbound_coefficient_Herb(surface_type); // (W/m2K)
      
      // inbound_flux = 
      // 	outbound_coefficient*sky_emissivity*pow(air_temperature+273.15,4)
      // 	+ 3.*outbound_coefficient*pow(old_surface_temperature+273.15,4)
      // 	- 4.*outbound_coefficient*pow(old_surface_temperature+273.15,3)*273.15;

      // outbound_flux =
      // 	4.*outbound_coefficient*pow(old_surface_temperature+273.15,3)*surface_temperature;
      
      inbound_flux
	=4.*outbound_coefficient*pow(((sky_temperature+(old_surface_temperature+273.15))/2.),3)
	*(pow(sky_emissivity,0.25)*air_temperature
	  +273.15*(pow(sky_emissivity,0.25)-1.));
      outbound_flux
	=4.*outbound_coefficient*pow(((sky_temperature+(old_surface_temperature+273.15))/2.),3)*old_surface_temperature;
    }
  else if (author=="Jansson")
    {
      outbound_coefficient
	=surface_coefficients
	.get_infrared_outbound_coefficient_Jansson(surface_type); // (W/m2K)

      // inbound_flux = 
      // 	outbound_coefficient*sky_emissivity*pow(air_temperature+273.15,4)
      // 	+ 3.*outbound_coefficient*pow(old_surface_temperature+273.15,4)
      // 	- 4.*outbound_coefficient*pow(old_surface_temperature+273.15,3)*273.15;
      // outbound_flux =
      // 	4.*outbound_coefficient*pow(old_surface_temperature+273.15,3)*surface_temperature;
      inbound_flux
	=4.*outbound_coefficient*pow(((sky_temperature+(old_surface_temperature+273.15))/2.),3)
	*(pow(sky_emissivity,0.25)*air_temperature
	  +273.15*(pow(sky_emissivity,0.25)-1.));
      
      outbound_flux
	=4.*outbound_coefficient*pow(((sky_temperature+(old_surface_temperature+273.15))/2.),3)*old_surface_temperature;
    }
  else if (author=="Best")
    {
      if ((canopy_density<0.) ||
	  (canopy_density>1.) || 
	  (canopy_temperature<-20.) ||
	  (canopy_temperature>100.))
	{
	  std::cout << "Error: out of range.   "     << std::endl
		    << "Error in get_infrared_flux." << std::endl
		    << "Author: "             << author             << std::endl
		    << "Canopy density    : " << canopy_density     << std::endl
		    << "Canopy temperature: " << canopy_temperature << std::endl;
	  throw 100;
	}

      outbound_coefficient
	=surface_coefficients
	.get_infrared_outbound_coefficient_Best(surface_type);  // (W/m2C)
      double outbound_coefficient_canopy
	=surface_coefficients.
	get_infrared_outbound_coefficient_Best("canopy");

      inbound_flux
	=(1.-canopy_density)*outbound_coefficient*sky_emissivity*pow(air_temperature+273.15,4)
	+3.*outbound_coefficient*pow(old_surface_temperature+273.15,4)
	-4.*outbound_coefficient*pow(old_surface_temperature+273.15,3)*273.15
	+canopy_density*outbound_coefficient_canopy*pow((canopy_temperature+273.15),4);
	
      outbound_flux
	=4.*outbound_coefficient*pow(old_surface_temperature+273.15,3)*surface_temperature;
    }
  else
    {
      std::cout << "Error: not implemented."           << std::endl
		<< "Error in get_total_infrared_flux." << std::endl
		<< "Author: " << author << std::endl;;
      throw 100;
    }
    
  if (direction=="both")
    return (inbound_flux - outbound_flux);
  else if (direction=="inbound")
    return (inbound_flux);
  else if (direction=="outbound")
    return (outbound_flux);
  else
    {
      std::cout << "Error: direction not implemented."    << std::endl
		<< "Error in get_infrared_coefficient." << std::endl
		<< "Direction: " << direction << std::endl;
      throw 100;
    }
}

double BoundaryConditions::get_convective_flux(const std::string surface_type,
					       const std::string author,
					       const std::string direction,
					       const double canopy_density,
					       const double old_surface_temperature,
					       const double new_surface_temperature,
					       const bool new_estimate)
{
  SurfaceCoefficients surface_coefficients(surface_pressure,
					   precipitation,
					   moisture_movement);
  double convective_coefficient=0.;
  
  if (author=="Jansson")
    convective_coefficient
      =surface_coefficients // (W/m2K)
      .get_convective_coefficient_Jansson(surface_type,
					  wind_speed);
  else if (author=="Herb")
    convective_coefficient
      =surface_coefficients // (W/m2K)
      .get_convective_coefficient_Herb(/*surface_type,*/
				       air_temperature,
				       relative_humidity,
				       wind_speed,
				       old_surface_temperature,
				       new_surface_temperature,
				       new_estimate);
  else if (author=="Best")
    {
      if ((canopy_density<0.) ||
	  (canopy_density>1.))
	{
	  std::cout << "Error: canopy density out of range.\n"
		    << "Error in get_inbound_convective_flux.\n";
	  throw 100;
	}
      double level_of_soil_evaporation=0.0;
      convective_coefficient
	=surface_coefficients // (W/m2K)
	.get_convective_coefficient_Best(/*surface_type,*/
					 air_temperature,
					 relative_humidity,
					 wind_speed,
					 old_surface_temperature,
					 new_surface_temperature,
					 new_estimate)
	*(1.-level_of_soil_evaporation*canopy_density);
    }
  else
    {
      std::cout << "Error: author not implemented."       << std::endl
		<< "Error in get_convective_coefficient." << std::endl;
      throw 100;
    }

  if (direction=="both")
    {
      if (new_estimate)
	return(convective_coefficient
	       *(air_temperature-new_surface_temperature)); // (W/m2)
      else
	return(convective_coefficient
	       *(air_temperature-old_surface_temperature)); // (W/m2)
    }
  else if (direction=="inbound")
    return (convective_coefficient*air_temperature); // (W/m2)
  else if (direction=="outbound")
    {
      if (new_estimate)
	return (convective_coefficient*new_surface_temperature); // (W/m2)
      else
	return (convective_coefficient*old_surface_temperature); // (W/m2)
    }
  else
    {
      std::cout << "Error: direction not implemented."    << std::endl
		<< "Error in get_convective_coefficient." << std::endl;
      throw 100;
    }
}

double BoundaryConditions::
get_evaporative_flux(const std::string surface_type,
		     const std::string author,
		     const double canopy_density,
		     const double old_surface_temperature,
		     const double new_surface_temperature,
		     const bool new_estimate)
{
  SurfaceCoefficients surface_coefficients(surface_pressure,
					   precipitation,
					   moisture_movement);
  double evaporative_flux=0.;
 
  if (author=="Jansson")
    {
      evaporative_flux
	=surface_coefficients
	.get_evaporative_flux_Jansson(surface_type,
				      air_temperature,
				      relative_humidity,
				      wind_speed,
				      air_temperature,//old_surface_temperature,
				      air_temperature,//new_surface_temperature,
				      new_estimate); // (W/m2K)
      return(evaporative_flux); // (W/m2)
    }
  else if (author=="Herb")
    {
      evaporative_flux = 
	surface_coefficients
	.get_evaporative_flux_Herb(/*surface_type,*/
				   air_temperature,
				   relative_humidity,
				   wind_speed,
				   old_surface_temperature,
				   new_surface_temperature,
				   new_estimate); // (W/m2)
      return (evaporative_flux); // (W/m2)
    }
  else if (author=="Best")
    {
      if ( (canopy_density<0.) ||
	   (canopy_density>1.))
	{
	  std::cout << "Error: canopy density out of range." << std::endl
		    << "Error in get_evaporative_flux."      << std::endl;
	  throw 100;
	}
      double level_of_soil_evaporation=0.0;
      evaporative_flux
	=surface_coefficients
	.get_evaporative_flux_Best(/*surface_type,*/
				   air_temperature,
				   relative_humidity,
				   wind_speed,
				   old_surface_temperature,
				   new_surface_temperature,
				   new_estimate); // (W/m2K)
      return((1.-level_of_soil_evaporation*canopy_density)*evaporative_flux); // (W/m2)
    }
  else
    {
      std::cout << "Error: author not implemented." << std::endl
		<< "Error in get_evaporative_flux." << std::endl;
      throw 100;
    }  
}

double BoundaryConditions::
get_canopy_temperature(/*const std::string surface_type,*/
		       const std::string author,
		       const double canopy_density)
{
  std::string dummy=author;
  dummy="";
  if ((canopy_density<0.) ||
      (canopy_density>1.))
    {
      std::cout << "Error: canopy density out of range." << std::endl
		<< "Error in get_canopy_temperature."    << std::endl;
      throw 100;
    }
  if (surface_temperature==-1000.)
    {
      std::cout << "Error: undefined variable."       << std::endl
		<< "Error in get_canopy_temperature." << std::endl;
      throw 100;
    }
    
  SurfaceCoefficients surface_coefficients(surface_pressure,
					   precipitation,
					   moisture_movement);
  
  double absortivity
    =surface_coefficients
    .get_absortivity_Best ("soil");
  double coeff_infrared_in_sky
    =surface_coefficients
    .get_infrared_inbound_coefficient_Best("canopy",
					   air_temperature,
					   relative_humidity);
  double coeff_infrared_in_soil
    =surface_coefficients
    .get_infrared_outbound_coefficient_Best("soil");
  double coeff_infrared_out_canopy
    =surface_coefficients
    .get_infrared_outbound_coefficient_Best("canopy");
  double coeff_convective
    =surface_coefficients
    .get_convective_coefficient_canopy_Best(/*surface_type,*/
					    wind_speed);
  double evaporative_heat_flux
    =surface_coefficients
    .get_evaporative_flux_canopy_Best(/*surface_type,*/
				      air_temperature,
				      relative_humidity,
				      wind_speed,
				      solar_radiation,
				      surface_temperature,
				      false);
  double evaporative_heat_flux_derivative
    =surface_coefficients
    .get_evaporative_flux_canopy_Best(/*surface_type,*/
				      air_temperature,
				      relative_humidity,
				      wind_speed,
				      solar_radiation,
				      surface_temperature,
				      true);
  double old_canopy_temperature=0.;
  double new_canopy_temperature=surface_temperature;
  
  double function=0.;
  double function_derivative=0.;
  
  while (fabs(old_canopy_temperature-new_canopy_temperature)>0.0001)
    {
      old_canopy_temperature=new_canopy_temperature;
      
      function
	=(canopy_density * absortivity * solar_radiation
	 + canopy_density * coeff_infrared_in_sky * pow(air_temperature+273.15,4)
	 + canopy_density * coeff_infrared_in_soil * pow(surface_temperature+273.15,4)
	 - 2. * canopy_density * coeff_infrared_out_canopy * pow(old_canopy_temperature+273.15,4)
	 + canopy_density * coeff_convective * (air_temperature - old_canopy_temperature)
	 + canopy_density * evaporative_heat_flux);
      
      function_derivative
	=(-8. * canopy_density * coeff_infrared_out_canopy * pow(old_canopy_temperature+273.15,3)
	  - canopy_density * coeff_convective
	 + canopy_density * evaporative_heat_flux_derivative);
      
      new_canopy_temperature
	=old_canopy_temperature - function/function_derivative;
    }
  return(new_canopy_temperature);
}

double BoundaryConditions::
get_inbound_heat_flux (const std::string surface_type,
		       const std::string author,
		       const double shading_factor,
		       const double canopy_temperature,
		       const double canopy_density,
		       const double old_surface_temperature,
		       const double new_surface_temperature,
		       const bool new_estimate)
{
  print_solar_radiation=0.;
  print_inbound_convective=0.;
  print_inbound_infrared=0.;
  print_inbound_evaporative=0.;
  
  bool active_evaporation=true;
  if (precipitation<=0 && surface_type=="road")
    active_evaporation=false;

  double inbound_heat_flux=0.;
  /*
    Add up inbound solar radiation. Taking care of reducing it by the
    provided shading factor
  */
  print_solar_radiation
    =(1.-shading_factor) *
    get_solar_flux (surface_type,
		    author,
		    canopy_density);
  
  inbound_heat_flux+=print_solar_radiation;
  /*
    Add up inbound evaporative heat flux (currenlty constant)
    Two options implemented: linear and exponential (standard)
    linear option is meant to be used to estimate only the 
    evaporative flux in the new timestep, not the previous
  */
  if (active_evaporation)
    {
      print_inbound_evaporative
	=get_evaporative_flux(surface_type,
			      author,
			      canopy_density,
			      old_surface_temperature,
			      new_surface_temperature,
			      new_estimate);
      inbound_heat_flux
	+=print_inbound_evaporative;
    }
  /*
    Add up inbound convective heat flux
  */
  print_inbound_convective
    =get_convective_flux(surface_type,
			 author,
			 "inbound",
			 canopy_density,
			 old_surface_temperature,
			 new_surface_temperature,
			 new_estimate);
  inbound_heat_flux
    +=print_inbound_convective;
  /*
    Add up inbound infrared heat flux.

    The overall idea is to get an inbound heat flux to be applied
    in the construction of the right hand side of the system and
    to obtain an outbound coefficient to be used in the construction
    of the system matrix (laplace matrix).

    The original infrared heat transfer equation is:
    
    Q = surface_emissivity*steffan_boltzmann*
    (sky_temperature_k^4 - surface_temperature_k^4) ------ (1)
    
    where the subindex '_k' indicates absolute temperature.
    The sky can be considered as a blackbody at some equivalent
    sky_temperature that is related to the air temperature throug 
    the next relation:
    
    sky_temperature^4 = sky_emissivity*air_temperature_k^4
    sky_temperature   = sky_emissivity^(1/4)*air_temperature_k
    
    The particular form of the sky_emissivity term depends on
    the formulation (author) chosen. For the purposes of 
    stability, and to allow use of conventional linear algebraic
    equation solvers for solving an equation that includes (1),
    we must first liearize it. This is best done using the
    Taylor series expansion:

    Xi^4 = (Xi_old)^4 + 4*(Xi_old)^3 * (Xi-Xi_old) -------- (2)

    where Xi_old is the value of Xi corresponding to the
    previous time step. Using (2) in (1):
    
    A = surface_emissivity*steffan_boltzmann

    Q = A*sky_temperature_k^4  - A*surface_temperature_k^4
    = A*sky_temperature_k^4 - 
    A*(surface_temperature_k_old^4 
    + 4*surface_temperature_k_old^3 
    * (surface_temperature_k-surface_temperature_k_old))
    
    = A*sky_temperature_k^4
    - A*surface_temperature_k_old^4 
    - 4*A*surface_temperature_k_old^3*surface_temperature_k
    + 4*A*surface_temperature_k_old^3*surface_temperature_k_old
    
    = A*sky_temperature_k^4
    - A*surface_temperature_k_old^4 
    - 4*A*surface_temperature_k_old^3*surface_temperature_k
    + 4*A*surface_temperature_k_old^4

    = A*sky_temperature_k^4
    + 3*A*surface_temperature_k_old^4
    - 4*A*surface_temperature_k_old^3*surface_temperature_k
    
    = A*sky_temperature_k^4
    + 3*A*surface_temperature_k_old^4
    - 4*A*surface_temperature_k_old^3*(surface_temperature+273.15)

    (Note the change from K to C in the surface temperature)

    = A*sky_temperature_k^4
    + 3*A*surface_temperature_k_old^4
    - 4*A*surface_temperature_k_old^3*273.15
    - 4*A*surface_temperature_k_old^3*surface_temperature

    Now we define the infrared inbound heat flux. The idea is to 
    put everything that does not depend on the current surface
    temperature is this term. Let me define it as:

    inbound_heat_flux = A*sky_temperature_k^4 
    + 3*A*surface_temperature_k_old^4
    - 4*A*surface_temperature_k_old^3*273.15
    
    inbound_heat_flux = A*sky_emissivity*air_temperature_k^4
    + 3*A*surface_temperature_k_old^4
    - 4*A*surface_temperature_k_old^3*273.15

    inbound_heat_flux = A*sky_emissivity*air_temperature_k^4
    + (A*surface_temperature_k_old^3) * 
    (3*surface_temperature_k_old - 4*273.15)

    In the infrared outbound heat flux we put everything that
    depends on the current surface temperature. 
    
    outbound_heat_flux
    = 4*A*surface_temperature_k_old^3*surface_temperature
    
    And the outbound heat transfer coefficient is defined as:
    
    outbound_coefficient 
    = 4*A*surface_temperature_k_old^3
    = 4*surface_emissivity*steffan_boltzmann*surface_temperature_k_old^3

    Note that both, the inbound heat flux and outbound coefficient 
    depend on the previous surface temperature. In the simulation
    process we use a Crank-Nicolson time discretization scheme. On it we
    need the heat flux contributions at the current (i) and previous (i-1)
    time step. And when assembling the infrared contribution for the 
    previous time step we should, in principle, use the (i-2) time step.
    However, as the main reason to use the previous time step in the
    infrared term is to avoid instability issues and as the solution at
    (i-1) has already converged and is unlikely that will give us
    stability problems, we will use the (i-1) surface temperature to
    construct that term as well.

    Now, finally, we add up the inbound infrared heat flux (Note that the
    formulation for 'Best' requires adding an extra term from the canopy
    cover):
  */
  print_inbound_infrared
    =get_infrared_flux(surface_type,
		       author,
		       "inbound",
		       old_surface_temperature,
		       canopy_temperature,
		       canopy_density);
  inbound_heat_flux
    +=print_inbound_infrared;

  return(inbound_heat_flux);
}

void BoundaryConditions::
print_inbound_heat_fluxes(double &inbound_solar_flux,
			  double &inbound_convective_flux,
			  double &inbound_evaporative_flux,
			  double &inbound_infrared_flux,
			  double &outbound_convective_coefficient,
			  double &outbound_infrared_coefficient,
			  double &outbound_evaporative_coefficient)
{
  inbound_solar_flux       = print_solar_radiation;
  inbound_convective_flux  = print_inbound_convective;
  inbound_evaporative_flux = print_inbound_evaporative;
  inbound_infrared_flux    = print_inbound_infrared;
  outbound_convective_coefficient  = print_outbound_convection_coefficient;
  outbound_infrared_coefficient    = print_outbound_infrared_coefficient;
  outbound_evaporative_coefficient = print_outbound_evaporative_coefficient;
}


double BoundaryConditions::
get_outbound_coefficient (const std::string surface_type,
			  const std::string author,
			  const double canopy_density,
			  const double old_surface_temperature,
			  const double new_surface_temperature,
			  const bool new_estimate)
{
  print_outbound_convection_coefficient  = 0.;
  print_outbound_infrared_coefficient    = 0.;
  print_outbound_evaporative_coefficient = 0.;

  bool active_evaporation=true;
  if (precipitation<=0 && surface_type=="road")
    active_evaporation=false;

  double outbound_coefficient=0.;
  SurfaceCoefficients surface_coefficients(surface_pressure,
					   precipitation,
					   moisture_movement);

  double sky_emissivity = 
    surface_coefficients
    .get_sky_emissivity(author,
			air_temperature,
			relative_humidity);
  double sky_temperature = 
    pow(sky_emissivity,0.25)*(air_temperature+273.15);
  
  if (author=="Jansson")
    {
      print_outbound_convection_coefficient
	=surface_coefficients
	.get_convective_coefficient_Jansson(surface_type,
					    wind_speed);
      print_outbound_infrared_coefficient
	=4.*pow(((sky_temperature+(old_surface_temperature+273.15))/2.),3)
	//4.*pow(old_surface_temperature+273.15,3)
	*surface_coefficients
	.get_infrared_outbound_coefficient_Jansson(surface_type);
    }
  else if (author == "Herb")
    {
      print_outbound_convection_coefficient
	= surface_coefficients
	.get_convective_coefficient_Herb (/*surface_type,*/
					  air_temperature,
					  relative_humidity,
					  wind_speed,
					  old_surface_temperature,
					  new_surface_temperature,
					  new_estimate);/* (W/m2C) */;
      
      print_outbound_infrared_coefficient
	= 4.*pow(((sky_temperature+(old_surface_temperature+273.15))/2.),3)
	//4.*pow(old_surface_temperature+273.15,3)
	*surface_coefficients
	.get_infrared_outbound_coefficient_Herb (surface_type);
    }
  else if (author=="Best")
    {
      double level_of_soil_evaporation=0.0;
      print_outbound_convection_coefficient
	=(1.-level_of_soil_evaporation*canopy_density) *
	surface_coefficients
	.get_convective_coefficient_Best (/*surface_type,*/
					  air_temperature,
					  relative_humidity,
					  wind_speed,
					  old_surface_temperature,
					  new_surface_temperature,
					  new_estimate);/* (W/m2C) */
      print_outbound_infrared_coefficient
	=/*(1.-canopy_density) */
	4.*pow(old_surface_temperature+273.15,3)
	*surface_coefficients
	.get_infrared_outbound_coefficient_Best (surface_type);
    }
  else
    {
      std::cout << "Error: Author not implemented."  << std::endl
		<< "Error in get_inbound_heat_flux." << std::endl;
      throw 100;
    }
  
  if (active_evaporation)
    {
      
      if (author=="Jansson")
	print_outbound_evaporative_coefficient
	  =surface_coefficients
	  .get_evaporative_coefficient_Jansson(surface_type,
					       wind_speed);
      else if (author == "Herb")
	print_outbound_evaporative_coefficient
	  =surface_coefficients
	  .get_evaporative_coefficient_Herb(/*surface_type,*/air_temperature,
					    relative_humidity,wind_speed,
					    old_surface_temperature,
					    new_surface_temperature,
					    new_estimate);
      else if (author == "Best")
	{
	  double level_of_soil_evaporation=0.0;
	  print_outbound_evaporative_coefficient
	    = (1.-level_of_soil_evaporation*canopy_density) *
	    surface_coefficients
	    .get_evaporative_coefficient_Best(/*surface_type,*/air_temperature,
					      relative_humidity,wind_speed,
					      old_surface_temperature,
					      new_surface_temperature,
					      new_estimate);
	}
      else
	{
	  std::cout << "Error: Author not implemented."  << std::endl
		    << "Error in get_inbound_heat_flux." << std::endl;
	  throw 100;
	}
    }
  
  outbound_coefficient += 
    print_outbound_convection_coefficient;
  outbound_coefficient += 
    print_outbound_infrared_coefficient;
  // outbound_coefficient +=
  //   print_outbound_evaporative_coefficient;
  return (outbound_coefficient);
}


double BoundaryConditions::
get_evaporative_mass_flux (const std::string surface_type,
			   const std::string author,
			   const double canopy_density,
			   const double old_surface_temperature,
			   const double new_surface_temperature,
			   const bool new_estimate)
{
  SurfaceCoefficients surface_coefficients(surface_pressure,
					   precipitation,
					   moisture_movement);
  double evaporative_mass_flux=0.;
  
  if (author=="Jansson")
    {
      evaporative_mass_flux
	=surface_coefficients
	.get_evaporative_mass_flux_Jansson(surface_type,
					   air_temperature,
					   relative_humidity,
					   wind_speed,
					   old_surface_temperature,
					   new_surface_temperature,
					   new_estimate); // (W/m2K)
    }
  else if (author == "Herb")
    {
      evaporative_mass_flux =
	surface_coefficients
	.get_evaporative_mass_flux_Herb(/*surface_type,*/
					air_temperature,
					relative_humidity,
					wind_speed,
					old_surface_temperature,
					new_surface_temperature,
					new_estimate); // (W/m2)
    }
  else if (author == "Best")
    {
      if ( (canopy_density < 0.) ||
	   (canopy_density > 1.))
	{
	  std::cout << "Error: canopy density out of range." << std::endl
		    << "Error in get_evaporative_flux."      << std::endl;
	  throw 100;
	}
      evaporative_mass_flux =
	(1.-canopy_density)*
	surface_coefficients
	.get_evaporative_mass_flux_Best (/*surface_type,*/
					 air_temperature,
					 relative_humidity,
					 wind_speed,
					 old_surface_temperature,
					 new_surface_temperature,
					 new_estimate); // (W/m2K)
    }
  else
    {
      std::cout << "Error: author not implemented." << std::endl
		<< "Error in get_evaporative_flux." << std::endl;
      throw 100;
    }
  return (evaporative_mass_flux); // (W/m2)
}
