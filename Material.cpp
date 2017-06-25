#include <string>
#include <map>
#include <vector>
#include <iostream>
#include "Material.h"

Material::Material(std::string material_name_)
{
  material_data["dummy_1"].push_back(1.);
  material_data["dummy_1"].push_back(2.);
  material_data["dummy_1"].push_back(3.);

  material_data["dummy_2"].push_back(4.);
  material_data["dummy_2"].push_back(5.);
  material_data["dummy_2"].push_back(6.);

  material_data["quartz_1"].push_back(8.79);
  material_data["quartz_1"].push_back(2660.);
  material_data["quartz_1"].push_back(2010.);

  material_data["pvc_1"].push_back(0.22);
  material_data["pvc_1"].push_back(1200.);
  material_data["pvc_1"].push_back(1200.);

  material_data["glass_beads"].push_back(0.8);
  material_data["glass_beads"].push_back(2500.);
  material_data["glass_beads"].push_back(1175.);

  material_data["pvc_2"].push_back(0.16);
  material_data["pvc_2"].push_back(1440.);
  material_data["pvc_2"].push_back(900.);

  material_data["clay_trl"           ].push_back(   1.2 );
  material_data["clay_trl"           ].push_back(1960.  );
  material_data["clay_trl"           ].push_back( 840.  );

  material_data["wearing_course"     ].push_back(   1.2  );
  material_data["wearing_course"     ].push_back(2400.   );
  material_data["wearing_course"     ].push_back(1200.   );

  material_data["binder_course"      ].push_back(   1.2  );
  material_data["binder_course"      ].push_back(2400.   );
  material_data["binder_course"      ].push_back(1200.   );

  material_data["concrete"           ].push_back(   1.4  );
  material_data["concrete"           ].push_back(2100.   );
  material_data["concrete"           ].push_back( 840.   );

  material_data["granular_material_1"].push_back(   1.4  );
  material_data["granular_material_1"].push_back(2100.   );
  material_data["granular_material_1"].push_back( 840.   );

  material_data["insulation_trl"     ].push_back(   0.034);
  material_data["insulation_trl"     ].push_back(  30.   );
  material_data["insulation_trl"     ].push_back(1130.   );

  material_data["sand_trl"           ].push_back(   0.33 );
  material_data["sand_trl"           ].push_back(2240.   );
  material_data["sand_trl"           ].push_back( 840.   );

  material_data["snow_bulk"          ].push_back(   0.15 );//snow (dry) 0.050-0.250 W/mK (from Wikipedia)
  material_data["snow_bulk"          ].push_back( 450.   );//100-800 kg/m3 (from Wikipedia)
  material_data["snow_bulk"          ].push_back(2090.   );

  material_data["ice"                ].push_back(   2.22 );//(@  0C) [W/mK]
  material_data["ice"                ].push_back( 920.   );//(@-30C)[kg/m3]
  material_data["ice"                ].push_back(1844.   );//(@-30C)[J/kgK]

  material_data["water"              ].push_back(   0.57 );
  material_data["water"              ].push_back(1000.   );
  material_data["water"              ].push_back(4186.   );

  material_data["air"                ].push_back(   0.025);
  material_data["air"                ].push_back(   1.25 );
  material_data["air"                ].push_back(   1.256);
  
  solids_thermal_conductivity=0.;
  solids_density=0.;
  solids_specific_heat_capacity=0.;
  
  std::map<std::string,std::vector<double> >::iterator it;
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

double Material::volumetric_heat_capacity(double temperature)
{
  temperature=temperature+0.; // Not used in this fuction but necessary in derived classes.
  return solids_specific_heat_capacity*solids_density;
}

double Material::thermal_diffusivity()
{
  return solids_thermal_conductivity/(solids_specific_heat_capacity*solids_density);
}
