#include "/home/zerpiko/Dropbox/PhD/libraries/AnalyticSolution.h"
#include "/home/zerpiko/Dropbox/PhD/libraries/SurfaceCoefficients.h"
#include <math.h>
#include <stdio.h>
//#include <vector>

AnalyticSolution::AnalyticSolution (const double thermal_conductivity_,
				    const double thermal_capacity_,
				    const double density_,
				    const std::string surface_type_,
				    const bool validation_
				    // const double thermal_conductivity_ =    1.2,
				    // const double thermal_capacity_     =  840.0,
				    // const double density_              = 1960.0,
				    // const std::string surface_type_    = "soil",
				    // const bool validation_             = false
				    )
{
  // surface_type (surface_type_),
  // validation   (validation_)

  surface_type = surface_type_;
  validation   = validation_;

  pi    = 3.14159265359;            // (dimensionless)
  phi   = 2.*pi/(24.*3600.);        // (1/s)
  theta = 2.*pi/(24.*3600.*365.25); // (1/s)

  if (!validation)
    {
      thermal_conductivity = thermal_conductivity_; // (W/mK)
      thermal_capacity     = thermal_capacity_;     // (J/kgK)
      density              = density_;              // (kg/m3)
      thermal_diffusivity  = (thermal_conductivity/
			      (density*thermal_capacity)); // (m2/s)
      
      Rsummer_daily_average = 204.23; // average of june     and july    monthly-daily averages 1985-2004 (W/m2)
      Rwinter_daily_average =  21.31; // average of december and january monthly-daily averages 1985-2004 (W/m2)
      Tave_ann    =  9.52;    // 20 years (1month average) // (C)
      Tamp_summer = 5.54/2.;  // monthly-daily averages 1985-2004 (C)
      Tamp_winter = 8.37/2.;  // monthly-daily averages 1985-2004 (C)
      Tave_summer = 15.44;    // monthly-daily averages 1985-2004 (C)
      Tave_winter =  3.60;    // monthly-daily averages 1985-2004 (C)
    }
  else
    {
      thermal_conductivity = 1;
      thermal_capacity     = 800;
      density              = 2000;
      thermal_diffusivity  = (thermal_conductivity/
			      (density*thermal_capacity));

      Rsummer_daily_average = 250;
      Rwinter_daily_average =  20;
      Tave_ann    =  9.8;
      Tamp_summer =  2.5;
      Tamp_winter =  5.;
      Tave_summer = 16.;
      Tave_winter =  3.60;
    }
  RadiationConstantA = (Rsummer_daily_average-Rwinter_daily_average)/2; // (W/m2)
  RadiationConstantB = (Rsummer_daily_average+Rwinter_daily_average)/2; // (W/m2)
  RadiationConstantC = (2.*pi)/(3.*pi+4.);                              // (dimensionless)
  RadiationConstantD = (4.-pi)/(2.*pi);                                 // (dimensionless)
  Rave_factor        = (3*pi+4)/4.;                                     // (dimensionless)
  Tamp_amp    = 0.5*(Tamp_summer - Tamp_winter); // (C)
  Tave_amp    = 0.5*(Tamp_summer + Tamp_winter); // (C)
  Tamp_ave    = 0.5*(Tave_summer - Tave_winter); // (C)
  Tave_ave    = 0.5*(Tave_summer + Tave_winter); // (C)
  
  air_temperature_correction_factor = 1;
}

double AnalyticSolution::get_value (const double depth,
				    const double time,
				    int iterations = 2500)
{
  double domain_lenght;        // (m)
  double initial_temperature;  // (C)
  double average_soil_surface_temperature; // (C)
  
  if (!validation)
    {
      domain_lenght       = 20.;
      initial_temperature = 10.9;
      //average_soil_surface_temperature = 8.234946;
      average_soil_surface_temperature = 9.604456;  // <<-- using emissivity with cloud factor
    }
  else
    {
      domain_lenght       = 20.;
      initial_temperature = 14.;
      //average_soil_surface_temperature = 8.903705;
      average_soil_surface_temperature = 10.12734; // <<-- using emissivity with cloud factor
    }
  
  double relative_humidity = get_analytic_relative_humidity (/*const double time*/);
  double wind_speed        = get_analytic_wind_speed        (/*const double time*/);
   
  //Jansson coefficients
  SurfaceCoefficients surface_coefficients;
  double absortivity                     =  surface_coefficients // (dimensionless)
    .get_absortivity_Jansson(surface_type);
  
  double outbound_convective_coefficient = surface_coefficients  // (W/m2K)
    .get_convective_coefficient_Jansson(surface_type,
					wind_speed)  
    + surface_coefficients
    .get_infrared_coefficient_Jansson(surface_type,
				      Tave_ann,
				      relative_humidity,
				      Tave_ann,
				      average_soil_surface_temperature);
  double inbound_convective_coefficient  = surface_coefficients  //(W/m2K)
    .get_convective_coefficient_Jansson(surface_type,
					wind_speed)  
    + (surface_coefficients
       .get_infrared_coefficient_Jansson(surface_type,
					 Tave_ann,
					 relative_humidity,
					 Tave_ann,
					 average_soil_surface_temperature)
       * pow( (surface_coefficients.get_infrared_inbound_coefficient_Jansson( surface_type,
									      relative_humidity,
									      Tave_ann)
	       / surface_coefficients.get_infrared_outbound_coefficient_Jansson (surface_type)),0.25));
  
  
  double evaporative_heat_flux  =  surface_coefficients          // (W/m2)
    .get_evaporative_flux_Jansson(surface_type,
				  Tave_ann,
				  wind_speed,
				  relative_humidity,
				  average_soil_surface_temperature); // (W/m2)

  double infrared_error = surface_coefficients
    .get_infrared_coefficient_Jansson(surface_type,
				      Tave_ann,
				      relative_humidity,
				      Tave_ann,
				      average_soil_surface_temperature)
    * 273.15 * (1-pow( (surface_coefficients.get_infrared_inbound_coefficient_Jansson( surface_type,
										       relative_humidity,
										       Tave_ann)
			/ surface_coefficients.get_infrared_outbound_coefficient_Jansson (surface_type)),0.25));
  
  double H1 = outbound_convective_coefficient/thermal_conductivity; // (1/m)
  double hi = inbound_convective_coefficient;                       // (W/m2K)  
  //--------------------------------------------------------------------------
  double Soil_Temperature = 0; // initialization of the temperature variable
  double beta_m = 0;           // eigenvalues used in the analytical solution
  
  double kernel_const = 0;
  double kernel_X     = 0;
  double kernel_0     = 0;
  double kernel_L_int = 0;
  
  double C0 = 0;
  double C1 = 0;
  
  if (iterations<=0)
    {
      std::cout << "Error. Wrong number of iterations" << std::endl;
      throw 100;
    }
  
  for (unsigned int i=1; i<=(unsigned)iterations; i++)
    {
      beta_m = (roots_BtanB(beta_m*domain_lenght,H1*domain_lenght)/
		domain_lenght);                                     // (1/m)
      
      kernel_const = sqrt(2 * (pow(H1,2) + pow(beta_m,2))/(H1 + domain_lenght*(pow(H1,2) + pow(beta_m,2)))); // (1/m^0.5)
      kernel_X     = kernel_const * cos(beta_m*(domain_lenght-depth));                                    // (1/m^0.5)
      kernel_0     = kernel_const * cos(beta_m*domain_lenght);                                            // (1/m^0.5)
      kernel_L_int = kernel_const * sin(beta_m*domain_lenght) * (1/beta_m);                               // (m^0.5)
      
      C0 = thermal_diffusivity*pow(beta_m,2); // (1/s)
      C1 = exp(-C0*time);       // (dimensionless)
      
      Soil_Temperature = Soil_Temperature
        + kernel_X * kernel_L_int * C1 * initial_temperature                            // this line corresponds to the initial condition (C)
        + kernel_X * kernel_0 * (thermal_diffusivity/thermal_conductivity) *            // the following lines correspond to the b.c. at x=0 (m2C/Ws)
	(air_temperature_correction_factor * 
	 (int_exp_cos_0_t  (-0.50*hi*Tamp_amp,C0,  phi - theta,time)//this line corresponds to the ambient temperature (Ws/m2)
	  + int_exp_cos_0_t(-0.50*hi*Tamp_amp,C0,  phi + theta,time)//this line corresponds to the ambient temperature (Ws/m2)
	  + int_exp_sin_0_t( 0.25*hi*Tamp_amp,C0,  phi - theta,time)//this line corresponds to the ambient temperature (Ws/m2)
	  + int_exp_sin_0_t(-0.25*hi*Tamp_amp,C0,  phi + theta,time)//this line corresponds to the ambient temperature (Ws/m2)
	  + int_exp_cos_0_t(-1.00*hi*Tave_amp,C0,  phi        ,time)//this line corresponds to the ambient temperature (Ws/m2)
	  + int_exp_cos_0_t( 1.00*hi*Tamp_ave,C0,        theta,time)//this line corresponds to the ambient temperature (Ws/m2)
	  + int_exp_sin_0_t( 0.50*hi*Tamp_ave,C0,        theta,time)//this line corresponds to the ambient temperature (Ws/m2)
	  + Tave_ave*hi*(1/C0)*(1-C1))                               //this line corresponds to the ambient temperature (Ws/m2)
	 + absortivity*Rave_factor*RadiationConstantC*
	 (int_exp_cos_0_t  ( 0.25*RadiationConstantA,C0,2*phi + theta,time)   //this line corresponds to radiation (Ws/m2)
	  + int_exp_cos_0_t( 0.25*RadiationConstantA,C0,2*phi - theta,time)   //this line corresponds to radiation (Ws/m2)
	  + int_exp_cos_0_t(-0.50*RadiationConstantA,C0,  phi + theta,time)   //this line corresponds to radiation (Ws/m2)
	  + int_exp_cos_0_t(-0.50*RadiationConstantA,C0,  phi - theta,time)   //this line corresponds to radiation (Ws/m2)
	  + int_exp_cos_0_t( 0.50*RadiationConstantA*(1+2*RadiationConstantD),C0,        theta,time)//this line corresponds to radiation (Ws/m2)
	  + int_exp_cos_0_t( 0.50*RadiationConstantB,C0,2*phi        ,time)   //this line corresponds to radiation (Ws/m2)
	  + int_exp_cos_0_t(     -RadiationConstantB,C0,  phi        ,time)   //this line corresponds to radiation (Ws/m2)
	  + 0.50*RadiationConstantB*(1 + 2*RadiationConstantD)*(1/C0)*(1-C1)) //this line corresponds to radiation  (Ws/m2)
	 + (evaporative_heat_flux-infrared_error)*(1/C0)*(1-C1)                                //this line corresponds to evaporation (Ws/m2)
	 );
    }
  return (Soil_Temperature);
}

double AnalyticSolution::get_integral_value(const double time,
					    const int iterations = 2500)
{
  double domain_lenght;       // (m)
  double initial_temperature; // (C)
  double average_soil_surface_temperature; // (C)
  
  if (!validation)
    {
      domain_lenght       = 12.875;
      initial_temperature = 10.9;
      //average_soil_surface_temperature = 8.234946;
      average_soil_surface_temperature = 9.604456; // <<-- using emissivity with cloud factor
    }
  else
    {
      domain_lenght       = 20.;
      initial_temperature = 14.;
      //average_soil_surface_temperature = 8.903705;
      average_soil_surface_temperature = 10.12734; // <<-- using emissivity with cloud factor
    }

  double relative_humidity = get_analytic_relative_humidity (/*const double time*/);
  double wind_speed        = get_analytic_wind_speed        (/*const double time*/);
  
  //Jansson coefficients
  SurfaceCoefficients surface_coefficients;
  double absortivity                     = surface_coefficients  // (dimensionless)
    .get_absortivity_Jansson(surface_type);
  double outbound_convective_coefficient = surface_coefficients  // (W/m2K)
    .get_convective_coefficient_Jansson(surface_type,
					wind_speed)  
    + surface_coefficients
    .get_infrared_coefficient_Jansson(surface_type,
				      Tave_ann,
				      relative_humidity,
				      Tave_ann,
				      average_soil_surface_temperature);
  double inbound_convective_coefficient  = surface_coefficients  // (W/m2K)
    .get_convective_coefficient_Jansson(surface_type,
					wind_speed)  
    + (surface_coefficients
       .get_infrared_coefficient_Jansson(surface_type,
					 Tave_ann,
					 relative_humidity,
					 Tave_ann,
					 average_soil_surface_temperature)
       * pow( (surface_coefficients.get_infrared_inbound_coefficient_Jansson( surface_type,
									      relative_humidity,
									      Tave_ann)
	       / surface_coefficients.get_infrared_outbound_coefficient_Jansson (surface_type)),0.25));
  double evaporative_heat_flux = surface_coefficients // (W/m2)
    .get_evaporative_flux_Jansson (surface_type,
				   Tave_ann,
				   wind_speed,
				   relative_humidity,
				   average_soil_surface_temperature); // (W/m2)
  double infrared_error = surface_coefficients
    .get_infrared_coefficient_Jansson(surface_type,
				      Tave_ann,
				      relative_humidity,
				      Tave_ann,
				      average_soil_surface_temperature)
    * 273.15 * (1-pow( (surface_coefficients.get_infrared_inbound_coefficient_Jansson( surface_type,
										       relative_humidity,
										       Tave_ann)
			/ surface_coefficients.get_infrared_outbound_coefficient_Jansson (surface_type)),0.25)); // (W/m2)

  double H1 = outbound_convective_coefficient/thermal_conductivity; // (1/m)
  double hi = inbound_convective_coefficient;                       // (W/mK)  
  //--------------------------------------------------------------------------
  double Soil_Temperature = 0; // initialization of the temperature variable
  double beta_m = 0;           // eigenvalues used in the analytical solution

  if (iterations<=0)
    {
      std::cout << "Error. Wrong number of iterations" << std::endl;
      std::cout << "error in get_integral_value" << std::endl;
      throw 100;
    } 
 
  for (unsigned int i=1; i<=(unsigned)iterations; i++)
    {
      
      beta_m = (roots_BtanB(beta_m*domain_lenght,H1*domain_lenght)/ 
		domain_lenght);                                     // (1/m)
      
      double kernel_const = sqrt(2 * (pow(H1,2) + pow(beta_m,2))/(H1 + domain_lenght*(pow(H1,2) + pow(beta_m,2)))); // (1/m^0.5)
      double int_kernel_X = (kernel_const/beta_m) * sin(beta_m*domain_lenght);                                    // (1/m^0.5)
      double kernel_0     = kernel_const * cos(beta_m*domain_lenght);                                            // (1/m^0.5)
      double kernel_L_int = kernel_const * sin(beta_m*domain_lenght) * (1/beta_m);                               // (m^0.5)
      
      double C0 = thermal_diffusivity*pow(beta_m,2); // (1/s)
      double C1 = exp(-C0*time);       // (dimensionless)
      
      Soil_Temperature = Soil_Temperature
        + int_kernel_X * kernel_L_int * C1 * initial_temperature                            // this line corresponds to the initial condition (C)
        + int_kernel_X * kernel_0 * (thermal_diffusivity/thermal_conductivity) *            // the following lines correspond to the b.c. at x=0 (m2C/Ws)
	(air_temperature_correction_factor * 
	 (int_exp_cos_0_t  (-0.50*hi*Tamp_amp,C0,  phi - theta,time)//this line corresponds to the ambient temperature (Ws/m2)
	  + int_exp_cos_0_t(-0.50*hi*Tamp_amp,C0,  phi + theta,time)//this line corresponds to the ambient temperature (Ws/m2)
	  + int_exp_sin_0_t( 0.25*hi*Tamp_amp,C0,  phi - theta,time)//this line corresponds to the ambient temperature (Ws/m2)
	  + int_exp_sin_0_t(-0.25*hi*Tamp_amp,C0,  phi + theta,time)//this line corresponds to the ambient temperature (Ws/m2)
	  + int_exp_cos_0_t(-1.00*hi*Tave_amp,C0,  phi        ,time)//this line corresponds to the ambient temperature (Ws/m2)
	  + int_exp_cos_0_t( 1.00*hi*Tamp_ave,C0,        theta,time)//this line corresponds to the ambient temperature (Ws/m2)
	  + int_exp_sin_0_t( 0.50*hi*Tamp_ave,C0,        theta,time)//this line corresponds to the ambient temperature (Ws/m2)
	  + Tave_ave*hi*(1/C0)*(1-C1))                               //this line corresponds to the ambient temperature (Ws/m2)
	 + absortivity*Rave_factor*RadiationConstantC*
	 (int_exp_cos_0_t  ( 0.25*RadiationConstantA,C0,2*phi + theta,time)   //this line corresponds to radiation (Ws/m2)
	  + int_exp_cos_0_t( 0.25*RadiationConstantA,C0,2*phi - theta,time)   //this line corresponds to radiation (Ws/m2)
	  + int_exp_cos_0_t(-0.50*RadiationConstantA,C0,  phi + theta,time)   //this line corresponds to radiation (Ws/m2)
	  + int_exp_cos_0_t(-0.50*RadiationConstantA,C0,  phi - theta,time)   //this line corresponds to radiation (Ws/m2)
	  + int_exp_cos_0_t( 0.50*RadiationConstantA*(1+2*RadiationConstantD),C0,        theta,time)//this line corresponds to radiation (Ws/m2)
	  + int_exp_cos_0_t( 0.50*RadiationConstantB,C0,2*phi        ,time)   //this line corresponds to radiation (Ws/m2)
	  + int_exp_cos_0_t(     -RadiationConstantB,C0,  phi        ,time)   //this line corresponds to radiation (Ws/m2)
	  + 0.50*RadiationConstantB*(1 + 2*RadiationConstantD)*(1/C0)*(1-C1)) //this line corresponds to radiation  (Ws/m2)
	 + (evaporative_heat_flux-infrared_error)*(1/C0)*(1-C1)                                //this line corresponds to evaporation (Ws/m2)
	 );
    }
  return (Soil_Temperature);
}

double AnalyticSolution::roots_BcotB (const double initial_value,
				      const double C)
{
  double error = 1000;
  
  double a = initial_value;
  double b = 0;
  double c = initial_value + 3.14159265359;
  
  double fb = 0;
  double fc = c*cos(c) + C*sin(c);
  
  while (fabs(c-a)>0.1)
    {
      b = (a+c)/2;
      fb = b*cos(b) + C*sin(b);
      
      if (fc*fb<0)
	a=b;
      else
	{
	  c=b;
	  fc=fb;
	}
    }
  
  double new_value = b;
  double f         = 0;
  double fx        = 0;
  double old_value = 0;

  while(error>0.000001)
    {
      old_value = new_value;    
      f  =  old_value*cos(old_value) + C*sin(old_value);
      fx = -old_value*sin(old_value) + (1+C)*cos(old_value);
      
      new_value = old_value - f/fx;
      
      error = new_value - old_value;
    }
  
  return (new_value);
}

double AnalyticSolution::roots_BtanB (const double initial_value,
				      const double C)
{
  double error = 1000;
  
  double a = initial_value;
  double b = 0;
  double c = initial_value + 3.14159265359;
  
  double fb = 0;
  double fc = c*sin(c) - C*cos(c);
  
  while (fabs(c-a)>0.1)
    {
      b = (a+c)/2;
      fb = b*sin(b) - C*cos(b);
      
      if (fb*fc<0)
	a=b;
      else
	{
	  c=b;
	  fc=fb;
	}
    }
  
  double new_value = b;
  double f         = 0;
  double fx        = 0;
  double old_value = 0;

  while(error>0.000001)
    {
      old_value = new_value;    
      f  = old_value*sin(old_value) - C*cos(old_value);
      fx = old_value*cos(old_value) + (1+C)*sin(old_value);
      
      new_value = old_value - f/fx;
      
      error = new_value - old_value;
    }
  
  return (new_value);
}

double AnalyticSolution::int_exp_cos_0_t(const double constant_1 /*(W/m2)*/,
					 const double constant_2 /*(1/s)*/,
					 const double period     /*(1/s)*/,
					 const double time       /*(s)*/)
{
  // Integral of C*exp(constant_2*time)*cos(period*time)
  // This function provides the result for the integral of the function
  //
  //  constat_1*int{ exp(constat_2*time')*cos(period*time') }dtime'
  //
  // where 'constant_1', 'constant_2' and 'period' are coefficients. 
  // The result is calculated from the exact solution:
  //
  //   exp(constant_2*t')*( constant_2*cos(period*time') + period*sin(period*time') ) / (constant_2^2 + time^2)
  //
  // evaluated between the limits time' = 0 and time' = time
  
  //NOTE I'm multiplying everithing for exp(-constant_2*time)
  
  return (constant_1*( (( constant_2*cos(period*time) + period*sin(period*time) )/(pow(constant_2,2) + pow(period,2))) -
		       (exp(-constant_2*time)*(constant_2/(pow(constant_2,2) + pow(period,2))))));
  
}

double AnalyticSolution::int_exp_sin_0_t(const double constant_1 /*(W/m2)*/,
					 const double constant_2 /*(1/s)*/,
					 const double period     /*(1/s)*/,
					 const double time       /*(s)*/)
{
  // Integral of constant_1*exp(constant_2*time)*sin(period*time)
  // This function provides the result for the integral of the function
  //
  //  constant_1*int{ exp(constant_2*time')*sin(period*time') }dtime'
  //
  // where 'constant_1', 'constant_2' and 'period' are coefficients. The result is calculated from
  // the exact solution:
  //
  //   exp(constant_2*time')*( constant_2*sin(period*time') - period*cos(period*time') ) / (constant_2^2 + period^2)
  //
  // evaluated between the limits time' = 0 and time' = time
  
  // NOTE I'm multiplying everything for exp(-constant_2*time) after performing the
  // integration
  
  return ( constant_1*( (( constant_2*sin(period*time) - period*cos(period*time) )/(pow(constant_2,2) + pow(period,2))) + 
			(exp(-constant_2*time)*(period/(pow(constant_2,2) + pow(period,2)))))); 
}


void AnalyticSolution::AnalyticSolution::output (const std::string &filename)
{
  FILE       *file_stream = fopen(filename.c_str(),"a");
  if (!file_stream)
    throw 2;
  std::cout << std::endl
	    << " \t name of output file: "
	    << filename << std::endl
	    << std::endl;  
  
  double time  = 3600;
  double L     = 15;
  double dx    = 0.01;
  double value = 0;
  
  for (unsigned int t=0; t<=7; t++)
    {
      for (unsigned int i=0; i<=(L*100); i++)
	{
	  value = get_value (i*dx,time + t*24*3600);
	  
	  fprintf (file_stream,"%f\t%f\n",
		   i*dx,value);
	}
      std::cout << "\t time =" << t  << std::endl;
    }
  if (fclose(file_stream))
    throw 3;
}

// return value from analytical equation for solar radiation for
// a given time in seconds
double AnalyticSolution::get_analytic_solar_radiation(const double time)
{    
  return ((RadiationConstantA*cos(theta*time)+RadiationConstantB) *
	  Rave_factor * RadiationConstantC * 
	  (pow(cos(phi*time),2) - cos(phi*time) + RadiationConstantD)); // (W/m2)
}
// return value from analytical equation for air temperature for
// a given time in seconds
double AnalyticSolution::get_analytic_air_temperature  (const double time)
{  
  return (air_temperature_correction_factor * ((Tamp_ave*(cos(theta*time) + 0.5*sin(theta*time)) + Tave_ave) // (C)
					       - (Tamp_amp*(cos(theta*time) + 0.5*sin(theta*time)) + Tave_amp) * cos(phi*time)));
}
// return value for average wind speed based on TRL measurements. Currently not analytical
// expresion for wind speed implemented
double AnalyticSolution::get_analytic_wind_speed (/*const double time*/)
{
  return (1.1413); // (m/s)
}
// return value for average relative humidity based on TRL measurements. Currently not analytical
// expresion for relative humidity implemented
double AnalyticSolution::get_analytic_relative_humidity (/*const double time*/)
{
  return (80.66); // (%)
}
// Precipitation is a variable particularlly difficult to provide an analytical expression for.
// So, I don't bother with trying I just return a 0 value. It is interesting to compare
// results with and without evaporation
double AnalyticSolution::get_analytic_precipitation (/*const double time*/)
{
  return (0); // (mm/time_step)
}

void AnalyticSolution::composite_region (const double t0,
					 const unsigned int layers,
					 const std::vector<double> X,
					 const std::vector<double> thermal_conductivity,
					 const std::vector<double> thermal_diffusivity,
					 std::vector< std::vector<double> > &system_matrix,
					 std::vector<double> &rhs)
{
  const unsigned int nodes  = layers+1;
  const double PI = 3.14159265359;
  const double T0 = 20.;
  const double TL = 0.;

  // std::cout << "Layers: " << layers << std::endl;
  // std::cout << "Nodes : " << nodes  << std::endl;
  // std::cout << "Time  : " << t0 << std::endl;

  std::vector<double> T;
  std::vector< std::vector<double> >
    dVi_1(layers-1,std::vector<double>(2,0)),
    dVi  (layers-1,std::vector<double>(2,0));

  for (unsigned int i=0; i<nodes; i++)
    {
      if (i==0)
  	T.push_back(T0);
      else if (i==nodes-1)
  	T.push_back(TL);
      else
  	T.push_back(1.0);
    }

  if ((X.size() != nodes) ||
      (T.size() != nodes))
    {
      std::cout << "X size : " << X.size() << std::endl
  		<< "T size : " << T.size() << std::endl;
      throw 100;
    }
  /*
    The following are the integrals of Vi's with respect of time
  */
  for (unsigned int i=1; i<layers; i++)
    {
      double SUM_dVi_1_0 = 0.;
      double SUM_dVi_1_1 = 0.;
      double SUM_dVi_0 = 0.;
      double SUM_dVi_1 = 0.;
      
      for (unsigned int k=1; k<4000; k++)
  	{
  	  double n = (double) k;
	  
  	  SUM_dVi_1_0 += pow(-1.,n)*
  	    (1.-exp(-thermal_diffusivity[i-1]*pow(n*PI/(X[i]-X[i-1]),2)*t0))/pow(n,2);
  	  SUM_dVi_1_1 += 
  	    (1.-exp(-thermal_diffusivity[i-1]*pow(n*PI/(X[i]-X[i-1]),2)*t0))/pow(n,2);
	  
  	  SUM_dVi_0 += 
  	    (1.-exp(-thermal_diffusivity[i]*pow(n*PI/(X[i+1]-X[i]),2)*t0))/pow(n,2);
  	  SUM_dVi_1 += pow(-1.,n)*
  	    (1.-exp(-thermal_diffusivity[i]*pow(n*PI/(X[i+1]-X[i]),2)*t0))/pow(n,2);
  	}
      
      dVi_1[i-1][0] = ( 1./(X[i]-X[i-1]) + (1./t0)*2.*(X[i]-X[i-1])/(thermal_diffusivity[i-1]*pow(PI,2)) * SUM_dVi_1_0);
      dVi_1[i-1][1] = ( 1./(X[i]-X[i-1]) + (1./t0)*2.*(X[i]-X[i-1])/(thermal_diffusivity[i-1]*pow(PI,2)) * SUM_dVi_1_1);
      
      dVi[i-1][0]   = (-1./(X[i+1]-X[i]) - (1./t0)*2.*(X[i+1]-X[i])/(thermal_diffusivity[i]*pow(PI,2)) * SUM_dVi_0);
      dVi[i-1][1]   = (-1./(X[i+1]-X[i]) - (1./t0)*2.*(X[i+1]-X[i])/(thermal_diffusivity[i]*pow(PI,2)) * SUM_dVi_1);
    }

  // std::cout << std::endl;
  // for (unsigned int i=0; i<dVi_1.size(); i++)
  //   std::cout << dVi_1[i][0] << "\t" << dVi_1[i][1] << std::endl
  // 	      << dVi  [i][0] << "\t" << dVi  [i][1] << std::endl
  // 	      << std::endl;
  // std::cout << std::endl;

  /*
    Construct the system of equations
    we have layers-1 equations since we have 
  */
  // std::vector< std::vector<double> > system_matrix;
  // std::vector< double > rhs;
  for (unsigned int i=1; i<(nodes-1); i++)
    {
      std::vector<double> system_row;
      if (i==1)
  	{
  	  rhs.push_back (-thermal_conductivity[i-1]*dVi_1[i-1][0] * T[i-1]);
	  
  	  system_row.push_back((thermal_conductivity[i  ]*dVi  [i-1][0] - 
  				thermal_conductivity[i-1]*dVi_1[i-1][1]) * T[i  ]);
  	  system_row.push_back(-thermal_conductivity[i  ]*dVi  [i-1][1] * T[i+1]);
	  
  	  for (unsigned int j=1; j<(nodes-3); j++)
  	    system_row.push_back(0.);
  	}
      else if (i==layers-1)
  	{
  	  rhs.push_back ( thermal_conductivity[i  ]*dVi  [i-1][1] * T[i+1]);
  	  for (unsigned int j=1; j<(nodes-3); j++)
  	    system_row.push_back(0.);
  	  system_row.push_back( thermal_conductivity[i-1]*dVi_1[i-1][0] * T[i-1]);
  	  system_row.push_back((thermal_conductivity[i  ]*dVi  [i-1][0] - 
  				thermal_conductivity[i-1]*dVi_1[i-1][1]) * T[i  ]);
  	}
      else
  	{
  	  rhs.push_back (0);
	  
  	  for (unsigned int j=1; j<(i-1); j++)
  	    system_row.push_back(0.);
  	  system_row.push_back( thermal_conductivity[i-1]*dVi_1[i-1][0] * T[i-1]);
  	  system_row.push_back((thermal_conductivity[i  ]*dVi  [i-1][0] - 
  				thermal_conductivity[i-1]*dVi_1[i-1][1]) * T[i  ]);
  	  system_row.push_back(-thermal_conductivity[i  ]*dVi  [i-1][1] * T[i+1]);
  	  for (unsigned int j=(i+1); j<(nodes-1); j++)
  	    system_row.push_back(0.);
  	}
      system_matrix.push_back(system_row);
    }


  /*
    Solve system matrix
  */
  // for (unsigned int i=0; i<system_matrix.size(); i++)
  //   {
  //     for (unsigned int j=0; j<system_matrix[i].size(); j++)
  // 	std::cout << system_matrix[i][j] << "\t";
  //     std::cout << std::endl;
  //   }

  // for (unsigned int i=0; i<rhs.size(); i++)
  //   std::cout << rhs[i] << "\t";
  // std::cout << std::endl;
}
