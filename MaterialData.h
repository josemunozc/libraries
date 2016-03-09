#ifndef __MATERIAL_DATA_INCLUDED__
#define __MATERIAL_DATA_INCLUDED__

/*----------------------------------------- material data -----------------------------------------*/
/*
  Currently these are assumed values. Eventually they're going to be updated according with 
  experimental results.
  
  A linear dependence for the clay thermal properties have been implemented using data from 
  Pielke 1984, table 11-3
  
  Clay soil   thermal conductivity     specific heat capacity     density     thermal diffusivity
  (40% pore space) (W/mK)                     (J/kgK)              (kg/m3)         (m2/s)
  dry               0.25                        890                 1600            0.18E-6
  mc=0.1            0.63                       1005                 1700            0.37E-6
  mc=0.2            1.12                       1172                 1800            0.53E-6
  mc=0.3            1.33                       1340                 1900            0.52E-6
  mc=0.4            1.58                       1550                 2000            0.51E-6

  The saturation moisture content for soil is 0.482 (Garrat 1994, Table 9). So, the previous
  table covers almost the whole range of moisture contents for clay. It can be seen that the
  values are almost linear, specially for the specific heat capacity and density. So, I'll
  use linear interpolation to calculate the properties for moisture content in the two ranges
  defined in table 7.
*/

#include <vector>
class MaterialData
{
public:
  MaterialData (int dim_, bool system_, double moisture_content_,
		bool moisture_active_=false);
  double get_soil_thermal_diffusivity  (const unsigned int material_id) const;
  double get_soil_heat_capacity        (const unsigned int material_id) const;
  double get_soil_density              (const unsigned int material_id) const;
  double get_soil_thermal_conductivity (const unsigned int material_id) const;
  
  double get_pipe_thermal_conductivity () const;
 private:
  void check_index (const unsigned int index) const;

  const int dim;
  const bool system;
  const double moisture_content;
  const bool moisture_active;
  unsigned int number_of_boundaries;
  unsigned int number_of_materials;
  double pipe_thermal_conductivity;
  std::vector<double> soil_thermal_conductivity;
  std::vector<double> soil_specific_heat_capacity;
  std::vector<double> soil_density;

  double clay_thermal_conductivity;
  double clay_specific_heat_capacity;
  double clay_density;
  double clay_thermal_diffusivity;
};

#endif
