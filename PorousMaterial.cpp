#include <string>
#include <math.h>
#include <iostream>
#include "PorousMaterial.h"

PorousMaterial::PorousMaterial(std::string material_name_,
                               double porosity_,
                               double degree_of_saturation_)
  :Material(material_name_),
   porosity(porosity_),
   degree_of_saturation(degree_of_saturation_){
  init_liquid_gas_and_ice_parameters();
  init_freezing_parameters();
}

PorousMaterial::PorousMaterial(double solids_thermal_conductivity_,
                               double solids_density_,
                               double solids_specific_heat_capacity_,
                               double porosity_,
                               double degree_of_saturation_)
  :Material(solids_thermal_conductivity_,
   solids_density_,
   solids_specific_heat_capacity_),
   porosity(porosity_),
   degree_of_saturation(degree_of_saturation_){
  init_liquid_gas_and_ice_parameters();
  init_freezing_parameters();
}

void PorousMaterial::init_liquid_gas_and_ice_parameters(){
  Material material("water");
  liquid_thermal_conductivity  =material.thermal_conductivity();// [W/mK]
  liquid_density               =material.density();// [kg/m3]
  liquid_specific_heat_capacity=material.specific_heat_capacity();// [J/kgK]

  material=Material("air");
  gas_thermal_conductivity  =material.thermal_conductivity();// [W/mK]
  gas_density               =material.density();// [kg/m3]
  gas_specific_heat_capacity=material.specific_heat_capacity();// [J/kgK]

  material=Material("ice");
  ice_thermal_conductivity  =material.thermal_conductivity();//(@  0C) [W/mK]
  ice_density               =material.density();//(@-30C)[kg/m3]
  ice_specific_heat_capacity=material.specific_heat_capacity();//(@-30C)[J/kgK]
}

void PorousMaterial::init_freezing_parameters(){
  freezing_point       =   273.15; // [K]
  coefficient_alpha    =   -10.0; // [nondimensional?]
  reference_temperature=   293.0; // [K] 123.4(C);//
  latent_heat_of_fusion=334000.0; // [J/K]
}

void PorousMaterial::set_liquid_properties(double liquid_thermal_conductivity_,
                                           double liquid_density_,
                                           double liquid_specific_heat_capacity_){
  liquid_thermal_conductivity=liquid_thermal_conductivity_;
  liquid_density=liquid_density_;
  liquid_specific_heat_capacity=liquid_specific_heat_capacity_;
}

void PorousMaterial::set_gas_properties(double gas_thermal_conductivity_,
                                        double gas_density_,
                                        double gas_specific_heat_capacity_){
  gas_thermal_conductivity=gas_thermal_conductivity_;
  gas_density=gas_density_;
  gas_specific_heat_capacity=gas_specific_heat_capacity_;
}

void PorousMaterial::set_ice_properties(double ice_thermal_conductivity_,
                                        double ice_density_,
                                        double ice_specific_heat_capacity_){
  ice_thermal_conductivity=ice_thermal_conductivity_;
  ice_density=ice_density_;
  ice_specific_heat_capacity=ice_specific_heat_capacity_;
}

void PorousMaterial::set_freezing_properties(double freezing_point_,
                                             double coefficient_alpha_,
                                             double reference_temperature_,
                                             double latent_heat_of_fusion_){
  freezing_point=freezing_point_;
  coefficient_alpha=coefficient_alpha_;
  reference_temperature=reference_temperature_;
  latent_heat_of_fusion=latent_heat_of_fusion_;
}

double PorousMaterial::thermal_conductivity(){
  return thermal_conductivity("donazzi");
}

double PorousMaterial::thermal_conductivity(const std::string relationship){
  double thermal_conductivity=0.;
  if (relationship.compare("donazzi")==0){
      /*
       * Estimation of thermal conductivity from Donazzi (1979)
       * Donazzi neglects the contribution of air. But the formulation
       * is applicable to unsaturated soils because it includes degree
       * of saturation.
       * */
      thermal_conductivity=1./(pow(1./liquid_thermal_conductivity,porosity)*
                           pow(1./solids_thermal_conductivity,1.-porosity)*
                           exp(3.08*(1.-degree_of_saturation)*porosity));
    }
  else if (relationship.compare("haigh")==0){
      /*
       * Estimation of thermal conductivity from Haigh (2012)
       *
       */
      double void_ratio=porosity/(1.-porosity);
      double E=(2.*void_ratio-1.)/3.;//xi
      double B=(1./3.)*acos((2.*(1.+3.*E)*(1.-degree_of_saturation)
                  -pow(1.+E,3.))/pow(1.+E,3.));
      double X=0.5*(1.+E)*(1.+cos(B)-pow(3.,0.5)*sin(B));
      double a_w=liquid_thermal_conductivity/solids_thermal_conductivity;
      double a_a=gas_thermal_conductivity/solids_thermal_conductivity;

      thermal_conductivity=1.58*solids_thermal_conductivity*(2.*pow(1.+E,2.)*
                           ((a_w/pow(1.-a_w,2.))*log(((1+E)+(a_w-1.)*X)/(E+a_w))+
                            (a_a/pow(1.-a_a,2.))*log((1.+E)/((1+E)+(a_a-1.)*X)))
                           +(2.*(1.+E)/((1.-a_w)*(1.-a_a)))*((a_w-a_a)*X-(1.-a_a)*a_w));
  }
  else if (relationship.compare("bulk")==0){
      /*
       * This only returns the same value that is given as input in the input
       * file. Is meant to quickly let you test values for this parameter or
       * for materials that are not porous (like plastics).
       */
      thermal_conductivity=solids_thermal_conductivity;
    }
  else{
      std::cout << "Error. Wrong thermal conductivity relationship "
                   "requested not currently implemented.\n";
      throw -1;
    }
  return thermal_conductivity;
}

double PorousMaterial::degree_of_saturation_ice(const double temperature){
  if (temperature+273.15<=freezing_point)
    return 1.-pow(1.-(temperature+273.15-freezing_point),coefficient_alpha);
  else
    return 0.;
}

double PorousMaterial::degree_of_saturation_ice_derivative(const double temperature){
  if (temperature+273.15<=freezing_point)
    return coefficient_alpha*pow(1.-(temperature+273.15-freezing_point),
                                 coefficient_alpha-1.);
  else
    return 0.;
}

double PorousMaterial::volumetric_heat_capacity(double temperature){
  double Hc=
    (1.-degree_of_saturation_ice(temperature))*
    porosity*degree_of_saturation*liquid_specific_heat_capacity*liquid_density
    +porosity*gas_specific_heat_capacity*gas_density*(1.-degree_of_saturation)
    +solids_specific_heat_capacity*solids_density*(1.-porosity)
    +porosity*degree_of_saturation*degree_of_saturation_ice(temperature)*
    ice_specific_heat_capacity*ice_density;
  double a=
    (temperature+273.15-reference_temperature)*
    (degree_of_saturation*ice_density*ice_specific_heat_capacity
     -degree_of_saturation*liquid_density*liquid_specific_heat_capacity);
  double b=
    degree_of_saturation*ice_density*latent_heat_of_fusion;
  return Hc+porosity*degree_of_saturation_ice_derivative(temperature)*(a-b);;
}

double PorousMaterial::thermal_energy(const double temperature){
  double Hc=
    (1.-degree_of_saturation_ice(temperature))*
    porosity*degree_of_saturation*liquid_specific_heat_capacity*liquid_density
    +porosity*gas_specific_heat_capacity*gas_density*(1.-degree_of_saturation)
    +solids_specific_heat_capacity*solids_density*(1.-porosity)
    +porosity*degree_of_saturation*degree_of_saturation_ice(temperature)*
    ice_specific_heat_capacity*ice_density;
  return (Hc*(temperature+273.15-reference_temperature)
          -latent_heat_of_fusion*porosity*degree_of_saturation*
          degree_of_saturation_ice(temperature)*ice_density);
}
