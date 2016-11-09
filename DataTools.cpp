#include "DataTools.h"
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <iostream>
#include <iomanip>
#include <math.h>

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

void DataTools::print_data(std::ostream &outFile,
			   const std::vector< std::vector<double> > &data,
			   std::vector< std::vector<int> >*date_and_time,
			   std::string type_of_data,
			   int lines_to_print)
{
  /*
   * This function can be used to print two types of data.
   * One that comes as regular numbers arranged in columns
   * and another that is correlated with a date and time
   * like (meteorological data). By default the function
   * assumes the former and only one container is required.
   * In the second case, two containers are required, one
   * with the data and one with the date and time, and they
   * need to be check if have the same size.
   **/
  std::vector< std::vector<int> >& date_and_time_ref=*date_and_time;
  if (type_of_data.compare("date_and_time")!=0 &&
      type_of_data.compare("data_only")!=0)
    {
      std::cout << "Error, wrong type of data specified in"
		<< "print_data function in DataTools.cpp" << std::endl;
      throw -1;
    }
  
  if (type_of_data.compare("date_and_time")==0 &&
      data.size()!=date_and_time_ref.size())
    {
      for (unsigned int i=0; i<date_and_time_ref[0].size(); i++)
	std::cout << date_and_time_ref[0][i] << "\t";
      for (unsigned int i=0; i<data[0].size(); i++)
	std::cout << data[0][i] << "\t";
      std::cout << std::endl;
      
      std::cout << "Error in print_data function. Vector mismatch." << std::endl;
      std::cout << "Vector 1 size : " << date_and_time_ref.size() << std::endl;
      std::cout << "Vector 2 size : " << data.size() << std::endl;      
      throw -1;
    }
  /*
   * If we call the function with a negative number of
   * lines to print, it will print them all.
   **/
  if (lines_to_print<0)
    lines_to_print=data.size();
  /*
   * If we are calling the function for date_and_time_ref data, we assume
   * 6 entries in each line of date_and_time_ref container, 3 for date and
   * 3 for time. Let's check for this
   **/
  if (type_of_data.compare("date_and_time")==0 &&
      date_and_time_ref[0].size()!=6)
    throw 4;
  /*
   * print formated output to selected ostream
   **/
  for (unsigned int i = 0; i<(unsigned int)lines_to_print; i++)
    {
      if (type_of_data.compare("date_and_time")==0)
	{
	  outFile << std::setw(2) << std::setfill('0') << date_and_time_ref[i][0] << "/" 
		  << std::setw(2) << std::setfill('0') << date_and_time_ref[i][1] << "/" 
		  << std::setw(2) << std::setfill('0') << date_and_time_ref[i][2] << "\t"
		  << std::setw(2) << std::setfill('0') << date_and_time_ref[i][3] << ":" 
		  << std::setw(2) << std::setfill('0') << date_and_time_ref[i][4] << ":" 
		  << std::setw(2) << std::setfill('0') << date_and_time_ref[i][5] << "\t";
	}
      outFile.setf( std::ios::fixed, std::ios::floatfield );
      for (unsigned int j=0; j<data[0].size(); j++)
	{
	  if(fabs(data[i][j])<0.1)
	    outFile << std::setw(7) << std::setfill(' ')
		    << std::scientific << data[i][j] << "\t";
	  else
	    outFile << std::setw(7) << std::setfill(' ')
		    << std::setprecision(2) << std::fixed << data[i][j] << "\t";
	}
      outFile << "\n";
    }
}

void DataTools::read_data(const std::vector< std::string >   &filenames,
			  std::vector< std::vector<int> >    &date_and_time,
			  std::vector< std::vector<double> > &data,
			  bool day_number_column=true)
{
  /*
   * We clear the containers that will store met
   * data and date data from each set of files
   **/
  data.clear();
  date_and_time.clear();
  /*
   * In this loop, we will go through each set of 
   * files, reading them, storing the entries in 
   * temporal vectors and storing those vectors in
   * the above containers
   **/
  for (unsigned int j=0; j<filenames.size(); j++)
    {
      /*
       * Open files and check that were correctly opened
       **/
      std::ifstream file (filenames[j].c_str());
      if (!file.is_open())
	{
	  std::cout << "Cannot open file: " << filenames[j] << std::endl;
	  throw 2;
	}
      /*
       * Define some variables that will help us with the
       * reading process
       **/
      std::string line_file;
      std::string file_token;
      /*
       * We read each line in each file. If EOF is reach
       * in any of the files, the loop is finished
       **/
      while (std::getline(file,line_file))
	{
	  /*
	   * there are two possible options (so far) for the file 
	   * we are passing to this function. 1) the file is a met
	   * file with date and time in the first and second column:
	   * dd/mm/yyyy\tHH/MM/SS\tD1\tD2\t...
	   * 2) the file is a general file with only numbers in
	   * columns:
	   * D1\tD2\tD3\t.....
	   * In the first case we look for the characters '\' or
	   * ':' in the first line (we are assuming that all lines
	   * follow the same pattern). If we find them then is the 
	   * first case, if not, then is the second.
	   **/
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
		    while (std::getline(file_iss1,file_token,'/'))
		      row_date_and_time.push_back(atoi(file_token.c_str()));
		  if (c==1)
		    while (std::getline(file_iss1,file_token,':'))
		      row_date_and_time.push_back(atoi(file_token.c_str()));
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
	      /*
	       * Store temporal vectors to corresponding containers
	       **/
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
      /*
       * close files and check that were succesfully closed
       **/
      file.close();
      if (file.is_open())
	throw 3;
    }
  /*
   * One more check, this time check if we store the same number of lines
   * (just in case)
   * if (data.size() != date_and_time.size())
   *   throw 4;
   **/
}

double DataTools::interpolate_data(std::vector< std::pair<double,double> > table,
				   const double x)
{
  const double INF = 1.e100;
  /*
   * Assumes that "table" is sorted by .first
   * Check if x is out of bound
   **/
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
  /*
   * INFINITY is defined in math.h in the glibc implementation
   **/
  it = std::lower_bound(table.begin(), table.end(), std::make_pair(x, -INF));
  // Corner case
  if (it == table.begin()) return (it->second);
  it2 = it;
  --it2;
  return (it2->second + (it->second - it2->second)*(x - it2->first)/(it->first - it2->first));
}


