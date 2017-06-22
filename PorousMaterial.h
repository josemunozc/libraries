#ifndef __POROUS_MATERIAL_INCLUDED__
#define __POROUS_MATERIAL_INCLUDED__
/*
* This class defines a material that in principle can be composed of three 
* constituents: solid, liquid and gas.
*/
#include "Material.h"

class PorousMaterial:public Material
{
public:
  PorousMaterial(std::string material_name_,
		 double porosity_=0.,
		 double degree_of_saturation_=0.);

  PorousMaterial(double solids_thermal_conductivity_,
  		 double solids_density_,
  		 double solids_specific_heat_capacity_,
		 double porosity_=0.,
		 double degree_of_saturation_=0.);

  virtual double thermal_conductivity();
  virtual double thermal_conductivity(const std::string relationship);
  virtual double volumetric_heat_capacity(double temperature);
  
  void set_liquid_properties(double liquid_thermal_conductivity_,
			     double liquid_density_,
			     double liquid_specific_heat_capacity_);
  void set_gas_properties(double gas_thermal_conductivity_,
			  double gas_density_,
			  double gas_specific_heat_capacity_);
  void set_ice_properties(double ice_thermal_conductivity_,
			  double ice_density_,
			  double ice_specific_heat_capacity_);
  void set_freezing_properties(double freezing_point_,
			       double coefficient_alpha_,
			       double reference_temperature_,
			       double latent_heat_of_fusion_);
  double degree_of_saturation_ice(const double temperature);
  double degree_of_saturation_ice_derivative(const double temperature);
  double thermal_energy(const double temperature);
  
private:
  void init_liquid_gas_and_ice_parameters();
  
  void init_freezing_parameters();
  
  double porosity;
  double degree_of_saturation;
  
  double liquid_thermal_conductivity;
  double liquid_density;
  double liquid_specific_heat_capacity;

  double gas_thermal_conductivity;
  double gas_density;
  double gas_specific_heat_capacity;

  double ice_thermal_conductivity;
  double ice_density;
  double ice_specific_heat_capacity;
  
  double freezing_point;
  double coefficient_alpha;
  double reference_temperature;
  double latent_heat_of_fusion;
};

#endif
