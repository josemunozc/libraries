#include <string>
#include <map>
#include <iostream>
#include "Material.h"


Material::Material(std::string material_name_)
{
  solids_thermal_conductivity=0.;
  solids_density=0.;
  solids_specific_heat_capacity=0.;
  
  std::map<std::string,std::array<double,3> >::iterator it;
  it = material_data.find(material_name_);
  if(it!=material_data.end())
    {
      solids_thermal_conductivity=it->second[0];
      solids_density=it->second[1];
      solids_specific_heat_capacity=it->second[2];
    }
  else
    {
      std::cout << "Error. Material "
		<< material_name_ << " not found\n";
    }
}

Material::Material(double solids_thermal_conductivity_,
		   double solids_density_,
		   double solids_specific_heat_capacity_)
  :
  solids_thermal_conductivity(solids_thermal_conductivity_),
  solids_density(solids_density_),
  solids_specific_heat_capacity(solids_specific_heat_capacity_)
{
  
}

double Material::thermal_conductivity()
{
  return solids_thermal_conductivity;
}

double Material::density()
{
  return solids_density;
}

double Material::specific_heat_capacity()
{
  return solids_specific_heat_capacity;
}

double Material::volumetric_heat_capacity()
{
  return solids_specific_heat_capacity*solids_density;
}

double Material::thermal_diffusivity()
{
  return solids_thermal_conductivity/(solids_specific_heat_capacity*solids_density);
}
