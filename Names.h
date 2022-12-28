#ifndef __NAMES_INCLUDED__
#define __NAMES_INCLUDED__

#include<vector>
#include<string>
#include <iostream> // cout

class Names{
  /*
   * This class defines soil and road depths available in each data file.
   * It also defines the file location of meteorological files and soil
   * temperature files.
   * In the constructor all this information is stored in vectors containing
   * strings with the location of the file.
   * Then these vectors are return by specific functions within the class.
   * THIS CLASS DOES NOT TEST THE IF THE FILE EXIST
   *
   * Update: I included a function that checks if the file exists.
   * Hopefully this will be useful. At the moment the problem is that
   * the locations are hardcoded and they vary depending on the user's
   * directory structure (e.g. office desktop vs cluster)
   */
public:
  Names (std::string input_path_);
  void get_depths(std::vector<double> &output_vector,
  std::string depths_location,
  const int number_of_entries=-1);
  void get_bh_temperatures_filenames(std::vector<std::string> &output_vector,
                                     std::string borehole_location,
                                     const int number_of_filenames=-1);
  void get_met_data_filenames(std::vector<std::string> &output_vector,
                              const int preheating_step,
                              std::string met_data_type="trl_met_data");

private:
  std::string input_path;
  std::vector<double> soil_depths;
  std::vector<double> road_depths;
  std::vector<std::string> bh_A_temperatures;
  std::vector<std::string> met_data_trl;
  std::vector<std::string> met_data_met_office;
  std::vector<std::string> met_data_badc_daily_averages;
};
#endif
