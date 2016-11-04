#include "DataTools.h"
#include <string>
#include <sstream>
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


void DataTools::read_data(const std::vector< std::string >   &filenames,
			  std::vector< std::vector<int> >    &date_and_time,
			  std::vector< std::vector<double> > &data,
			  bool day_number_column = true)// by default we assume that the file has a column for the day number
{
	/*
    We clear the containers that will store met
    data and date data from each set of files
	 */
	data.clear();
	date_and_time.clear();
	/*
    In this loop, we will go through each set of 
    files, reading them, storing the entries in 
    temporal vectors and storing those vectors in
    the above containers
	 */
	for (unsigned int j=0; j<filenames.size(); j++)
	{
		// Open files and check that were correctly opened
		std::ifstream file (filenames[j].c_str());
		if (!file.is_open())
		{
			std::cout << "Cannot open file: " << filenames[j] << std::endl;
			throw 2;
		}
		// Define some variables that will help us with the
		// reading process
		std::string line_file;
		std::string file_token;
		// We read each line in each file. If EOF is reach
		// in any of the files, the loop is finished
		while (std::getline(file,line_file))
		{
			/*
	    there are two possible options (so far) for the file 
	    we are passing to this function. 1) the file is a met
	    file with date and time in the first and second column:
	    dd/mm/yyyy\tHH/MM/SS\tD1\tD2\t...
	    2) the file is a general file with only numbers in
	    columns:
	    D1\tD2\tD3\t.....
	    In the first case we look for the characters '\' or
	    ':' in the first line (we are assuming that all lines
	    follow the same pattern). If we find them then is the 
	    first case, if not, then is the second.
			 */
			std::stringstream file_iss;
			file_iss << line_file;
			std::vector<int>    row_date_and_time;
			std::vector<double> row_data;
			if ((line_file.find('/') != std::string::npos) ||
					(line_file.find(':') != std::string::npos))
			{
				std::stringstream file_iss1;
				int c = 0;

				while (std::getline(file_iss,file_token,'\t'))
				{
					file_iss1 << file_token;
					if (c==0)
						while (std::getline(file_iss1,file_token,'/')) row_date_and_time.push_back(atoi(file_token.c_str()));
					if (c==1)
						while (std::getline(file_iss1,file_token,':')) row_date_and_time.push_back(atoi(file_token.c_str()));
					if (c>1)
						while (std::getline(file_iss1,file_token,':'))
						{
							if (day_number_column && c>2)
								row_data.push_back(atof(file_token.c_str()));
							else if (!day_number_column)
								row_data.push_back(atof(file_token.c_str()));
						}
					file_iss1.clear();
					c++;
				}
				file_iss.clear();
				// Store temporal vectors to corresponding containers
				date_and_time.push_back (row_date_and_time);
				data.push_back          (row_data);
			}
			else
			{
				while (std::getline(file_iss,file_token,'\t'))
					row_data.push_back(atof(file_token.c_str()));
				data.push_back(row_data);
			}
		}

		// close files and check that were succesfully closed
		file.close();
		if (file.is_open())
			throw 3;
	}
	// One more check, this time check if we store the same number of lines
	// (just in case)
	// if (data.size() != date_and_time.size())
	//   throw 4;
}
