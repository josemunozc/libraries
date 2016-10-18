#include "DataTools.h"
#include <string>
#include <fstream>
#include <vector>
#include <iostream>

DataTools::DataTools()
{

}

void DataTools::open_file(std::ofstream &file,
			  const std::string filename)
{
  file.open(filename.c_str());
  if (!file.is_open())
    throw 2;
}

void DataTools::close_file(std::ofstream &file)
{
  file.close();
  if (file.is_open())
    throw 3;
  file.clear();
}

double DataTools::interpolate_data(std::vector< std::pair<double,double> > table,
				   const double x)
{
  const double INF = 1.e100;
  // Assumes that "table" is sorted by .first
  // Check if x is out of bound
  if (x > table.back().first)
    {
      std::cout << "Warning. In data_tools.h -> interpolate_data. value out of bound." 
		<< std::endl;
      return (INF);
    }
  if (x < table[0].first)
    {
      std::cout << "Warning. In data_tools.h -> interpolate_data. value out of bound." 
		<< std::endl;
      return (-INF);
    }
  std::vector<std::pair<double, double> >::iterator it, it2;
  // INFINITY is defined in math.h in the glibc implementation
  it = std::lower_bound(table.begin(), table.end(), std::make_pair(x, -INF));
  // Corner case
  if (it == table.begin()) return (it->second);
  it2 = it;
  --it2;
  return (it2->second + (it->second - it2->second)*(x - it2->first)/(it->first - it2->first));
}
