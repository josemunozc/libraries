#include <Names.h>

void Names::get_depths(std::vector<double> &output_vector,
		       std::string depths_location="road",
		       const int number_of_values)
{
  std::vector<double> data_vector;
  if (depths_location=="road")
    data_vector=road_depths;
  else if (depths_location=="soil")
    data_vector=soil_depths;
  else
    {
      std::cout << "Error in Names::get_depths. depths_location"
		<< " not implemented\n";
      throw -1;
    }

  if (number_of_values<0)
    output_vector=data_vector;
  else
    {
      if (number_of_values>(signed)data_vector.size())
	{
	  std::cout << "Error in Names::get_depths. Too many"
		    << " values requested\n";
	  throw -1;
	}
      for (int i=0; i<number_of_values; i++)
	output_vector.push_back(data_vector[i]);
    }
}

void Names::get_bh_temperatures_filenames(std::vector<std::string> &output_vector,
					  std::string borehole_location="A",
					  const int number_of_filenames)
{
  std::vector<std::string> data_vector;
  if (borehole_location=="A")
    data_vector=bh_A_temperatures;
  else
    {
      std::cout << "Error in Names::get_bh_temperatures_filenames."
		<< " bh_location not implemented\n";
    }
  
  if (number_of_filenames<0)
    output_vector=data_vector;
  else
    {
      if (number_of_filenames>(signed)data_vector.size())
	{
	  std::cout << "Error in Names::get_depths. Too many"
		    << " values requested\n";
	  throw -1;
	}
      for (int i=0; i<number_of_filenames; i++)
	output_vector.push_back(data_vector[i]);
    }
}

void Names::get_met_data_filenames(std::vector<std::string> &output_vector,
				   const int preheating_step,
				   std::string met_data_type)
{
  std::vector<std::string> data_vector;
  if (met_data_type=="met_office_data")
    data_vector=met_data_met_office;
  else if (met_data_type=="trl_met_data")
    data_vector=met_data_trl;
  else if (met_data_type=="badc_daily_met_data")
    data_vector=met_data_badc_daily_averages;
  else
    {
      std::cout << "Error in Names::get_met_data_filenames."
		<< " met_data_type not implemented\n";
      throw -1;
    }

  if ((preheating_step<0) ||
      (preheating_step>8))
    {
      std::cout << "Error in Names::get_met_data_filenames."
		<< " preheating_step not implemented.\n";
      throw -1;
    }

  if ((met_data_type=="met_office_data") &&
      (preheating_step>4))
    {
      std::cout << "Error in Names::get_met_data_filenames."
		<< " preheating_step not implemented for met_office_data.\n";
      throw -1;
    }

  if (met_data_type=="badc_daily_met_data")
    {
      std::cout << "Error in Names::get_met_data_filenames."
		<< " badc_daily_met_data not yet implemented.\n";
      throw -1;
    }

  if (preheating_step==1)
    for (unsigned int i=0; i<12; i++)
      output_vector.push_back(data_vector[i]);
  else if (preheating_step==2)
    for (unsigned int i=0; i<8; i++)
      output_vector.push_back(data_vector[i]);
  else if (preheating_step==3)
    {
      output_vector.push_back(data_vector[8]);
      output_vector.push_back(data_vector[9]);
      output_vector.push_back(data_vector[10]);
      output_vector.push_back(data_vector[12]);
    }
  else if (preheating_step==4)
    output_vector.push_back(data_vector[13]);
  else if (preheating_step==5)
    output_vector.push_back(data_vector[14]);
  else if (preheating_step==6)
    output_vector.push_back(data_vector[15]);
  else if (preheating_step==7)
    output_vector.push_back(data_vector[16]);
  else if (preheating_step==8)
    output_vector.push_back(data_vector[17]);
}

Names::Names (std::string input_path_):
  input_path(input_path_)
{
  soil_depths.push_back(  0.025 );
  soil_depths.push_back(  0.125 );
  soil_depths.push_back(  0.825 );
  soil_depths.push_back(  0.875 );
  soil_depths.push_back(  1.025 );
  soil_depths.push_back(  1.175 );
  soil_depths.push_back(  1.375 );
  soil_depths.push_back(  1.875 );
  soil_depths.push_back(  3.875 );
  soil_depths.push_back(  7.875 );
  soil_depths.push_back( 12.875 );

  road_depths.push_back(  0.000); // 0
  road_depths.push_back(  0.010); // 1
  road_depths.push_back(  0.025); // 2
  road_depths.push_back(  0.050); // 3
  road_depths.push_back(  0.075); // 4
  road_depths.push_back(  0.100); // 5

  road_depths.push_back(  0.120); // 6
  road_depths.push_back(  0.1325); // 7
  road_depths.push_back(  0.160); // 8
  road_depths.push_back(  0.180); // 9

  road_depths.push_back(  0.200); // 10
  road_depths.push_back(  0.250); // 11
  road_depths.push_back(  0.300); // 12
  road_depths.push_back(  0.350); // 13
  road_depths.push_back(  0.400); // 14
  road_depths.push_back(  0.450); // 15
  road_depths.push_back(  0.500); // 16
  road_depths.push_back(  0.550); // 17
  road_depths.push_back(  0.600); // 18
  road_depths.push_back(  0.650); // 19
  road_depths.push_back(  0.700); // 20
  road_depths.push_back(  0.750); // 21
  road_depths.push_back(  0.800); // 22
  road_depths.push_back(  0.8475); // 23

  road_depths.push_back(  0.875); // 24
  road_depths.push_back(  0.925); // 25
  road_depths.push_back(  0.975); // 26
  road_depths.push_back(  1.025); // 27
  road_depths.push_back(  1.175); // 28
  road_depths.push_back(  1.375); // 29
  road_depths.push_back(  1.875); // 30
  road_depths.push_back(  2.875); // 31
  road_depths.push_back(  3.875); // 32
  road_depths.push_back(  7.875); // 33
  road_depths.push_back( 12.875); // 34
  /*
    Borehole temperatures
  */  
  bh_A_temperatures
    .push_back (input_path+"/trl_borehole_data/borehole_A_average_data/average_borehole_A_05_09.txt");
  bh_A_temperatures
    .push_back (input_path+"/trl_borehole_data/borehole_A_average_data/average_borehole_A_05_10.txt");
  bh_A_temperatures
    .push_back (input_path+"/trl_borehole_data/borehole_A_average_data/average_borehole_A_05_11.txt");
  bh_A_temperatures
    .push_back (input_path+"/trl_borehole_data/borehole_A_average_data/average_borehole_A_05_12.txt");
  bh_A_temperatures
    .push_back (input_path+"/trl_borehole_data/borehole_A_average_data/average_borehole_A_06_01.txt");
  bh_A_temperatures
    .push_back (input_path+"/trl_borehole_data/borehole_A_average_data/average_borehole_A_06_02.txt");
  bh_A_temperatures
    .push_back (input_path+"/trl_borehole_data/borehole_A_average_data/average_borehole_A_06_03.txt");
  bh_A_temperatures
    .push_back (input_path+"/trl_borehole_data/borehole_A_average_data/average_borehole_A_06_04.txt");
  bh_A_temperatures
    .push_back (input_path+"/trl_borehole_data/borehole_A_average_data/average_borehole_A_06_05.txt");
  bh_A_temperatures
    .push_back (input_path+"/trl_borehole_data/borehole_A_average_data/average_borehole_A_06_06.txt");
  bh_A_temperatures
    .push_back (input_path+"/trl_borehole_data/borehole_A_average_data/average_borehole_A_06_07.txt");
  bh_A_temperatures
    .push_back (input_path+"/trl_borehole_data/borehole_A_average_data/average_borehole_A_06_08.txt");
  /*
    Metereological data from Met Office
  */
  met_data_met_office.push_back (input_path+"/met_data/met_office/met_office_2005_09.txt"); // 0
  met_data_met_office.push_back (input_path+"/met_data/met_office/met_office_2005_10_using_2006_data.txt"); // 1
  met_data_met_office.push_back (input_path+"/met_data/met_office/met_office_2005_11.txt"); // 2
  met_data_met_office.push_back (input_path+"/met_data/met_office/met_office_2005_12.txt"); // 3
  met_data_met_office.push_back (input_path+"/met_data/met_office/met_office_2006_01.txt"); // 4
  met_data_met_office.push_back (input_path+"/met_data/met_office/met_office_2006_02.txt"); // 5
  met_data_met_office.push_back (input_path+"/met_data/met_office/met_office_2006_03.txt"); // 6
  met_data_met_office.push_back (input_path+"/met_data/met_office/met_office_2006_04.txt"); // 7
  met_data_met_office.push_back (input_path+"/met_data/met_office/met_office_2006_05.txt"); // 8
  met_data_met_office.push_back (input_path+"/met_data/met_office/met_office_2006_06.txt"); // 9
  met_data_met_office.push_back (input_path+"/met_data/met_office/met_office_2006_07.txt"); // 10
  met_data_met_office.push_back (input_path+"/met_data/met_office/met_office_2006_08.txt"); // 11
  met_data_met_office.push_back (input_path+"/met_data/met_office/met_office_2006_08_up_to_day_22.txt"); // 12
  met_data_met_office.push_back (input_path+"/met_data/met_office/hourly_average_met_office_data_from_23_08_2005_to_14_11_2005.txt"); // 13
  /*
    Metereological data from TRL
  */
  met_data_trl.push_back (input_path+"/met_data/trl/hourly_average_met_data_05_09.txt"); // 0
  met_data_trl.push_back (input_path+"/met_data/trl/hourly_average_met_data_05_10.txt"); // 1
  met_data_trl.push_back (input_path+"/met_data/trl/hourly_average_met_data_05_11.txt"); // 2
  met_data_trl.push_back (input_path+"/met_data/trl/hourly_average_met_data_05_12.txt"); // 3
  met_data_trl.push_back (input_path+"/met_data/trl/hourly_average_met_data_06_01.txt"); // 4
  met_data_trl.push_back (input_path+"/met_data/trl/hourly_average_met_data_06_02.txt"); // 5
  met_data_trl.push_back (input_path+"/met_data/met_office/met_office_2006_03.txt"); // 6
  //<<-- replace this file since the correspoding TRL's is from 2007 and it seems that that year was 
  // particularly hot
  //met_data_trl.push_back (input_path+"/met_data/trl/hourly_average_met_data_06_03.txt"); // 6
  met_data_trl.push_back (input_path+"/met_data/trl/hourly_average_met_data_06_04.txt"); // 7
  met_data_trl.push_back (input_path+"/met_data/trl/hourly_average_met_data_06_05.txt"); // 8
  met_data_trl.push_back (input_path+"/met_data/trl/hourly_average_met_data_06_06.txt"); // 9
  met_data_trl.push_back (input_path+"/met_data/trl/hourly_average_met_data_06_07.txt"); // 10
  met_data_trl.push_back (input_path+"/met_data/met_office/met_office_2006_08.txt");     // 11
  met_data_trl.push_back (input_path+"/met_data/met_office/met_office_2006_08_up_to_day_22.txt"); // 12
  met_data_trl.push_back (input_path+"/met_data/trl/hourly_average_met_data_from_23_08_2005_to_14_11_2005.txt"); // 13
  met_data_trl.push_back (input_path+"/met_data/trl/hourly_average_met_data_from_15_11_2005_to_20_02_2006.txt"); // 14
  met_data_trl.push_back (input_path+"/met_data/trl/hourly_average_met_data_from_21_02_2006_to_26_04_2006.txt"); // 15
  met_data_trl.push_back (input_path+"/met_data/trl/hourly_average_met_data_from_27_04_2006_to_31_10_2006.txt"); // 16
  met_data_trl.push_back (input_path+"/met_data/trl/hourly_average_met_data_from_01_11_2006_to_28_02_2007.txt"); // 17
  /*
    DAILY averages from the BADC
  */
  met_data_badc_daily_averages.push_back(input_path+"/met_data/met_office/daily_averages_solar_air_temperature_1985.dat");
  met_data_badc_daily_averages.push_back(input_path+"/met_data/met_office/daily_averages_solar_air_temperature_1986.dat");
  met_data_badc_daily_averages.push_back(input_path+"/met_data/met_office/daily_averages_solar_air_temperature_1987.dat");
  met_data_badc_daily_averages.push_back(input_path+"/met_data/met_office/daily_averages_solar_air_temperature_1988.dat");
  met_data_badc_daily_averages.push_back(input_path+"/met_data/met_office/daily_averages_solar_air_temperature_1989.dat");
  met_data_badc_daily_averages.push_back(input_path+"/met_data/met_office/daily_averages_solar_air_temperature_1990.dat");
  met_data_badc_daily_averages.push_back(input_path+"/met_data/met_office/daily_averages_solar_air_temperature_1991.dat");
  met_data_badc_daily_averages.push_back(input_path+"/met_data/met_office/daily_averages_solar_air_temperature_1992.dat");
  met_data_badc_daily_averages.push_back(input_path+"/met_data/met_office/daily_averages_solar_air_temperature_1993.dat");
  met_data_badc_daily_averages.push_back(input_path+"/met_data/met_office/daily_averages_solar_air_temperature_1994.dat");
  met_data_badc_daily_averages.push_back(input_path+"/met_data/met_office/daily_averages_solar_air_temperature_1995.dat");
  met_data_badc_daily_averages.push_back(input_path+"/met_data/met_office/daily_averages_solar_air_temperature_1996.dat");
  met_data_badc_daily_averages.push_back(input_path+"/met_data/met_office/daily_averages_solar_air_temperature_1997.dat");
  met_data_badc_daily_averages.push_back(input_path+"/met_data/met_office/daily_averages_solar_air_temperature_1998.dat");
  met_data_badc_daily_averages.push_back(input_path+"/met_data/met_office/daily_averages_solar_air_temperature_1999.dat");
  met_data_badc_daily_averages.push_back(input_path+"/met_data/met_office/daily_averages_solar_air_temperature_2000.dat");
  met_data_badc_daily_averages.push_back(input_path+"/met_data/met_office/daily_averages_solar_air_temperature_2001.dat");
  met_data_badc_daily_averages.push_back(input_path+"/met_data/met_office/daily_averages_solar_air_temperature_2002.dat");
  met_data_badc_daily_averages.push_back(input_path+"/met_data/met_office/daily_averages_solar_air_temperature_2003.dat");
  met_data_badc_daily_averages.push_back(input_path+"/met_data/met_office/daily_averages_solar_air_temperature_2004.dat");
  met_data_badc_daily_averages.push_back(input_path+"/met_data/met_office/daily_averages_solar_air_temperature_2005.dat");
  met_data_badc_daily_averages.push_back(input_path+"/met_data/met_office/daily_averages_solar_air_temperature_2006.dat");
}
