#include "MaterialData.h"
#include "iostream"
#include "vector"
#include "math.h"

MaterialData::
MaterialData (int dim_,bool system_,double moisture_content_, bool moisture_active_):
  dim (dim_),
  system (system_),
  moisture_content (moisture_content_),
  moisture_active (moisture_active_),
  pipe_thermal_conductivity (0.4)
{
  /*
    Define ranges for the different materials (volumes in 3d, surfaces in 2d). This
    depends on the imported mesh. Different meshes will have different number of
    materials and different indexes. In the current case we also have a certain 
    set of indexes that belong to boundaries (surfaces in 2d, lines in 2d). The
    boundaries come before the materials in the ordering. The tables we are using
    to store the values for thermal properties have indexes that run from 0 to N
    where N is the number of materials. For this reason, when we receive the
    material number from the main program we need to substract the number of
    boundaries.
    Next we define the number of boundaries and the number of materials for a 2d
    and 3d cases.
  */
  if (dim==1)
    {
      number_of_boundaries = 2;
      number_of_materials  = 1;
    }
  else if (dim==2)
    {
      number_of_boundaries = 7; // Boundary indexes from 1 to 7
      number_of_materials  = 7; // Material indexes from 8 to 14
    }
  else if (dim==3)
    {
      number_of_boundaries = 11; // Boundary indexes from  1 to 11
      number_of_materials  =  9; // Material indexes from 12 to 20
    }
  else
    {
      std::cout << "Error in MaterialData: Not implemented. dim= "
		<< dim << std::endl;
      throw 100;
    }
  /*
    Define the thermal properties of clay based on the moisture
    content and the values presented at the beginning of 
    Material_Data.h
  */
  if(moisture_active==true)
    {
      double MC[5] = {    0.00   ,    0.10   ,    0.20    ,    0.30    ,     0.40};
      double K[5]  = {    0.25   ,    0.63   ,    1.12    ,    1.33    ,     1.58};
      double C[5]  = {  890.00   , 1005.00   , 1172.00    , 1340.00    ,  1550.00};
      double D[5]  = { 1600.00   , 1700.00   , 1800.00    , 1900.00    ,  2000.00};
      double A[5]  = {    0.18E-6,    0.37E-6,    0.53E-6 ,    0.52E-6 ,     0.51E-6};
   
      if ((moisture_content<MC[0])||
	  ((moisture_content-MC[4])>0.0001))
	std::cout << "Error in MaterialData: Moisture content out of range\n"
		  << "Moisture content: " << moisture_content 
		  << " is out of range: " << MC[0] << "---" << MC[4] << std::endl;
    
      for (unsigned int i=1; i<5; i++)
	{
	  if (moisture_content<MC[i] && moisture_content>=MC[i-1])
	    {
	      clay_thermal_conductivity
		=((K[i]-K[i-1])/(MC[i]-MC[i-1]))*(moisture_content-MC[i-1])+K[i-1];
	      clay_specific_heat_capacity
		=((C[i]-C[i-1])/(MC[i]-MC[i-1]))*(moisture_content-MC[i-1])+C[i-1];
	      clay_density
		=((D[i]-D[i-1])/(MC[i]-MC[i-1]))*(moisture_content-MC[i-1])+D[i-1];
	      clay_thermal_diffusivity
		=((A[i]-A[i-1])/(MC[i]-MC[i-1]))*(moisture_content-MC[i-1])+A[i-1];
	      break;
	    }
	}
    }
  else
    {
      clay_thermal_conductivity
	=1.2;
      clay_specific_heat_capacity
	=840.;
      clay_density
	=1960.;
      clay_thermal_diffusivity
	=clay_thermal_conductivity/(clay_specific_heat_capacity*clay_density);
    }

  soil_thermal_conductivity.resize   (number_of_materials); // Indexes from 0 to number_of_materials-1
  soil_specific_heat_capacity.resize (number_of_materials); // Indexes from 0 to number_of_materials-1
  soil_density.resize                (number_of_materials); // Indexes from 0 to number_of_materials-1

  
  if (dim==1)
    {
      soil_thermal_conductivity[0]//soil (1.21 silty clay, assuming same value for backfill and unperturbed soil)
	=clay_thermal_conductivity;
      soil_specific_heat_capacity[0]//soil (840  silty clay, assuming same value for backfill and unperturbed soil)
	=clay_specific_heat_capacity;
      soil_density[0]// soil (1960 silty clay, assuming same value for backfill and unperturbed soil)
	=clay_density;
    }
  else if (dim==2)
    {
      /*********************************
       Thermal Conductivity in 2D (W/mK)
      **********************************/
      soil_thermal_conductivity[0]=1.200; // wearing course 0.85
      soil_thermal_conductivity[1]=1.200; // binder course 0.85
      soil_thermal_conductivity[6]=clay_thermal_conductivity;// soil (1.21 silty clay, assuming same value for backfill and unperturbed soil)
      if (system)
	{
	  /*If the pipe arrays and insulation layer are present*/
	  soil_thermal_conductivity[2]=1.400;//collector pipes (concrete) 1.4
	  soil_thermal_conductivity[3]=1.400;//granular material type 1   1.4
	  soil_thermal_conductivity[4]=0.034;//insulation                 0.034
	  soil_thermal_conductivity[5]=0.330;//storage pipes (sand)       0.330
	}
      else
	{
	  /*If they are not present*/
	  soil_thermal_conductivity[2]=1.400;//concrete 1.4
	  soil_thermal_conductivity[3]=1.400;//concrete 1.4
	  soil_thermal_conductivity[4]=clay_thermal_conductivity; // soil (1.210 silty clay, assuming same value for backfill and unperturbed soil)
	  soil_thermal_conductivity[5]=clay_thermal_conductivity; // soil (1.210 silty clay, assuming same value for backfill and unperturbed soil)
	}
      /***********************************
       Specific Heat Capacity in 2D (J/kgK)
      ************************************/
      soil_specific_heat_capacity[0]=1200;//wearing course 850
      soil_specific_heat_capacity[1]=1200;//binder course 850
      soil_specific_heat_capacity[6]=clay_specific_heat_capacity; // soil (840 silty clay, assuming same value for backfill and unperturbed soil)
      if (system)
	{
	  /*If the pipe arrays and insulation layer are present*/
	  soil_specific_heat_capacity[2]= 840;// collector pipes (concrete)
	  soil_specific_heat_capacity[3]= 840;// granular material type 1
	  soil_specific_heat_capacity[4]=1130;// insulation 1130
	  soil_specific_heat_capacity[5]= 840;// storage pipes (sand);
	}
      else
	{
	  /*If they are not present*/
	  soil_specific_heat_capacity[2]=840;// concrete
	  soil_specific_heat_capacity[3]=840;// concrete
	  soil_specific_heat_capacity[4]=clay_specific_heat_capacity; // soil (840 silty clay, assuming same value for backfill and unperturbed soil)
	  soil_specific_heat_capacity[5]=clay_specific_heat_capacity; // soil (840 silty clay, assuming same value for backfill and unperturbed soil)
	}
      /**********************
       Density in 2D (kg/m3)
      ***********************/
      soil_density[0]=2400;//wearing course 2400
      soil_density[1]=2400;//binder course  2400
      soil_density[6]=clay_density;//soil (1960 silty clay, assuming same value for backfill and unperturbed soil)
      if (system)
	{
	  /*If the pipe arrays and insulation layer are present*/
	  soil_density[2]=2100;//collector pipes (concrete) 2100
	  soil_density[3]=2100;//granular material type 1   2100
	  soil_density[4]=  30;//insulation 30
	  soil_density[5]=2240;//storage pipes (sand)       2240
	}
      else
	{
	  /*If they are not present*/
	  soil_density[2]=2100;//concrete 2100
	  soil_density[3]=2100;//concrete 2100
	  soil_density[4]=clay_density;//soil (1960 silty clay, assuming same value for backfill and unperturbed soil)
	  soil_density[5]=clay_density;//soil (1960 silty clay, assuming same value for backfill and unperturbed soil)
	}
    }
  else if (dim==3)
    {
      /*********************************
       Thermal Conductivity in 3D (W/mK)
      **********************************/
      soil_thermal_conductivity[0] = clay_thermal_conductivity; // soil (1.210 silty clay, assuming same value for backfill and unperturbed soil)
      soil_thermal_conductivity[7] = 0.850; // binder course 0.85
      soil_thermal_conductivity[8] = 0.850; // wearing course 0.85
      if (system)
	{
	  /*If the pipe arrays and insulation layer are present*/
	  soil_thermal_conductivity[1] = 0.330; // storage pipes (sand above pipes)
	  soil_thermal_conductivity[2] = 0.330; // storage pipes (sand)
	  soil_thermal_conductivity[3] = 0.034; // insulation
	  soil_thermal_conductivity[4] = 1.400; // granular material type 1
	  soil_thermal_conductivity[5] = 1.400; // collector pipes (concrete above pipes)
	  soil_thermal_conductivity[6] = 1.400; // collector pipes (concrete)
	}
      else
	{
	  /*If they are not present*/
	  soil_thermal_conductivity[1] = clay_thermal_conductivity; // soil (1.210 silty clay, assuming same value for backfill and unperturbed soil)
	  soil_thermal_conductivity[2] = clay_thermal_conductivity; // soil (1.210 silty clay, assuming same value for backfill and unperturbed soil)
	  soil_thermal_conductivity[3] = clay_thermal_conductivity; // soil (1.210 silty clay, assuming same value for backfill and unperturbed soil)
	  soil_thermal_conductivity[4] = 1.400; // concrete
	  soil_thermal_conductivity[5] = 1.400; // concrete
	  soil_thermal_conductivity[6] = 1.400; // concrete
	}
      /***********************************
       Specific Heat Capacity in 3D (J/kgK)
      ************************************/
      soil_specific_heat_capacity[0] =  clay_specific_heat_capacity; // soil (840 silty clay, assuming same value for backfill and unperturbed soil)
      soil_specific_heat_capacity[7] =  850; // binder course
      soil_specific_heat_capacity[8] =  850; // wearing course
      if (system)
	{
	  /*If the pipe arrays and insulation layer are present*/
	  soil_specific_heat_capacity[1] =  840; // storage pipes (sand above pipes);
	  soil_specific_heat_capacity[2] =  840; // storage pipes (sand);
	  soil_specific_heat_capacity[3] = 1130; // insulation
	  soil_specific_heat_capacity[4] =  840; // granular material type 1
	  soil_specific_heat_capacity[5] =  840; // collector pipes (concrete above pipes)
	  soil_specific_heat_capacity[6] =  840; // collector pipes (concrete)

	}
      else
	{
	  /*If they are not present*/
	  soil_specific_heat_capacity[1] =  clay_specific_heat_capacity; // soil (840 silty clay, assuming same value for backfill and unperturbed soil)
	  soil_specific_heat_capacity[2] =  clay_specific_heat_capacity; // soil (840 silty clay, assuming same value for backfill and unperturbed soil)
	  soil_specific_heat_capacity[3] =  clay_specific_heat_capacity; // soil (840 silty clay, assuming same value for backfill and unperturbed soil)
	  soil_specific_heat_capacity[4] =  840; // concrete
	  soil_specific_heat_capacity[5] =  840; // concrete
	  soil_specific_heat_capacity[6] =  840; // concrete
	}
      /**********************
       Density in 3D (kg/m3)
      ***********************/
      soil_density[0] = clay_density; // soil (1960 silty clay, assuming same value for backfill and unperturbed soil)
      soil_density[7] = 2400; // binder course
      soil_density[8] = 2400; // wearing course
      if (system)
	{
	  /*If the pipe arrays and insulation layer are present*/
	  soil_density[1] = 2240; // storage pipes (sand above pipes)
	  soil_density[2] = 2240; // storage pipes (sand)
	  soil_density[3] = 30;   // insulation
	  soil_density[4] = 2100; // granular material type 1
	  soil_density[5] = 2100; // collector pipes (concrete above pipes)
	  soil_density[6] = 2100; // collector pipes (concrete)
	}
      else
	{
	  /*If they are not present*/
	  soil_density[1] = clay_density; // soil (1960 silty clay, assuming same value for backfill and unperturbed soil)
	  soil_density[2] = clay_density; // soil (1960 silty clay, assuming same value for backfill and unperturbed soil)
	  soil_density[3] = clay_density; // soil (1960 silty clay, assuming same value for backfill and unperturbed soil)
	  soil_density[4] = 2100; // concrete
	  soil_density[5] = 2100; // concrete
	  soil_density[6] = 2100; // concrete
	}
    }
  else
    {
      std::cout << "Error in MaterialData: Not implemented" << std::endl;
      throw 100;
    }
}

void MaterialData::
check_index (const unsigned int index) const {
  if (dim==1)
    {
      if (index!=0)
	std::cout << "Error in MaterialData.\n"
		  << "Wrong index: "
		  << index << ".\n" 
		  << "Only index = 0 available"
		  << " for dim=" << dim << std::endl;
    }
  else
    {
      if ((index<0) ||
	  (index>number_of_materials-1))
	std::cout << "Error in MaterialData.\n"
		  << "Index: "
		  << index << " out of range: " 
		  << "0 --- " << (index>number_of_materials-1)
		  << " for dim=" << dim << std::endl;
    }
}

double MaterialData::
get_soil_thermal_diffusivity (const unsigned int material_id) const {
  unsigned int index = material_id-number_of_boundaries-1;
  check_index (index);
  return ((soil_thermal_conductivity[index])/(soil_density[index]*soil_specific_heat_capacity[index]));
}

double MaterialData::
get_soil_heat_capacity (const unsigned int material_id) const {
  unsigned int index = material_id-number_of_boundaries-1;
  check_index (index);
  return (soil_specific_heat_capacity[index]);
}

double MaterialData::
get_soil_density (const unsigned int material_id) const {
  unsigned int index = material_id-number_of_boundaries-1;
  check_index (index);
  return (soil_density[index]);
}

double MaterialData::
get_soil_thermal_conductivity (const unsigned int material_id) const {
  unsigned int index = material_id-number_of_boundaries-1;
  check_index (index);
  return (soil_thermal_conductivity[index]);
}

double MaterialData::
get_pipe_thermal_conductivity () const {
  return (pipe_thermal_conductivity);
}
