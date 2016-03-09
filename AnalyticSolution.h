#ifndef __ANALYTIC_SOLUTION_INCLUDED__
#define __ANALYTIC_SOLUTION_INCLUDED__

#include <string>
#include <vector>
#include "SurfaceCoefficients.h"

class AnalyticSolution
{
 public:
  AnalyticSolution (const double thermal_conductivity_ =    1.2,
		    const double thermal_capacity_     =  840.0,
		    const double density_              = 1960.0,
		    const std::string surface_type_    = "soil",
		    const bool heat_source_            = false,
		    const bool validation_             = false,
		    const int type_of_weather_         = 1);
  double get_value (const double depth,
		    const double time,
		    int iterations = 5000);
  
  double get_integral_value (const double time,
			     int iterations);
  
  double get_analytic_solar_radiation   (const double time);
  double get_analytic_air_temperature   (const double time);
  double get_analytic_wind_speed        (/*const double time*/);
  double get_analytic_relative_humidity (/*const double time*/);
  double get_analytic_precipitation     (/*const double time*/);
  
  void output (const std::string &filename);  
  
  void composite_region (const double time,
			 const unsigned int layers,
			 const std::vector<double> X,
			 const std::vector<double> thermal_conductivity,
			 const std::vector<double> thermal_diffusivity,
			 std::vector< std::vector<double> > &matrix,
			 std::vector<double> &right_hand_side);
  
 private:  
  double roots_BcotB     (const double initial_value,
			  const double constant);
  
  double roots_BtanB     (const double initial_value,
			  const double constant);
  
  double int_exp_cos_0_t (const double constat_1,
			  const double constant_2,
			  const double period,
			  const double time);
  
  double int_exp_sin_0_t (const double constat_1,
			  const double constant_2,
			  const double period,
			  const double time);

  double thermal_conductivity; // soil thermal conductivity (W/mK)
  double thermal_capacity;     // soil heat capacity (J/kg)
  double density;              // soil density (kg/m3)
  double thermal_diffusivity;  // soil thermal diffusivity (m2/s)
  std::string surface_type;    // type of surface required
  bool heat_source;
  bool validation;
  int type_of_weather;
  
  double pi;
  double phi;   // period for radiation and ambient temperature functions (1/s)
  double theta; // period for annual variation (1/s)
  
  double Rsummer_ave;
  double Rwinter_ave;
  double Rannual_ave;
  
  double Rsummer_daily_average;
  double Rwinter_daily_average;
  
  double RadiationConstantA;
  double RadiationConstantB;
  double RadiationConstantC;
  double RadiationConstantD;
  double Rave_factor;

  double Tave_ann;
  double Tmax_summer;
  double Tmin_summer;
  double Tmax_winter;
  double Tmin_winter;

  double Tamp_summer;
  double Tamp_winter;
  double Tave_summer;
  double Tave_winter;
  double Tamp_amp;
  double Tave_amp;
  double Tamp_ave;
  double Tave_ave;
  
  double air_temperature_correction_factor;
};

#endif
