#ifndef __MATERIAL_INCLUDED__
#define __MATERIAL_INCLUDED__

#include <map>

class Material
{
public:
  Material(std::string material_name);
  Material(double solids_thermal_conductivity_,
	   double solids_density_,
	   double solids_specific_heat_capacity_);
  virtual double thermal_conductivity();
  virtual double density();
  virtual double specific_heat_capacity();
  virtual double volumetric_heat_capacity();
  virtual double thermal_diffusivity();
protected:
  double solids_thermal_conductivity;
  double solids_density;
  double solids_specific_heat_capacity;
  /*
   * material name - thermal conductivity [W/mK] - density [kg/m3] - specific heat capacity [J/kgK]
   */
  std::map<std::string,std::array<double,3> > material_data
  ={
    {"dummy_1"            ,{1.   , 2.     ,    3.  }},
    {"dummy_2"            ,{4.   , 5.     ,    6.  }},
    {"quartz_1"           ,{8.79 , 2660.  , 2010.  }},
    {"pvc_1"              ,{0.22 , 1200.  , 1200.  }},
    {"glass_beads"        ,{0.8  , 2500.  , 1175.  }},
    {"pvc_2"              ,{0.16 , 1440.  ,  900.  }},
    {"clay_trl"           ,{1.21 , 1960.  ,  840.  }},
    {"wearing_course"     ,{1.2  , 2400.  , 1200.  }},
    {"binder_course"      ,{1.2  , 2400.  , 1200.  }},
    {"concrete"           ,{1.4  , 2100.  ,  840.  }},
    {"granular_material_1",{1.4  , 2100.  ,  840.  }},
    {"insulation_trl"     ,{0.034,   30.  , 1130.  }},
    {"sand_trl"           ,{0.33 , 2240.  ,  840.  }}
  };
};

#endif
