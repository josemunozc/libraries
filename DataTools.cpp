#include "DataTools.h"

DataTools::DataTools()
{

}

void DataTools::open_file(std::ofstream &file,const std::string filename)
{
	file.open(filename.c_str());
	if (!file.is_open())
    throw 2;
}

void DataTools::open_file(std::ifstream &file,const std::string filename)
{
	file.open(filename.c_str());
	if (!file.is_open())
	{
		std::cout << "Error opening file: " << filename << "\n";
		throw 2;
	}
}

void DataTools::close_file(std::ofstream &file)
{
  file.close();
  if (file.is_open())
    throw 3;
  file.clear();
}

void DataTools::close_file(std::ifstream &file)
{
  file.close();
  if (file.is_open())
    throw 3;
  file.clear();
}

void DataTools::print_data(std::ostream &outFile,
			   const std::vector< std::vector<double> > &data,
			   int lines_to_print)
{
	if (lines_to_print<0)
		lines_to_print=data.size();

	if (data.size()!=0)
		for (int i=0; i<lines_to_print; i++)
		{
			outFile.setf(std::ios::fixed, std::ios::floatfield);
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
	else
	{
		std::cout << "In DataTools::print_data. Matrix to print is empty\n";
	}
}

void DataTools::print_data(std::ostream &outFile,
			   const std::vector< std::vector<double> > &data,
			   std::vector< std::vector<int> > &date_and_time,
			   int lines_to_print)
{
  if (data.size()!=date_and_time.size())
    {
	  for (unsigned int i=0; i<date_and_time[0].size(); i++)
		  std::cout << date_and_time[0][i] << "\t";
	  for (unsigned int i=0; i<data[0].size(); i++)
		  std::cout << data[0][i] << "\t";
	  std::cout << std::endl;
      
      std::cout << "Error in print_data function. Vector mismatch." << std::endl;
      std::cout << "Vector 1 size : " << date_and_time.size() << std::endl;
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
  if (date_and_time[0].size()!=6)
    throw 4;
  /*
   * print formated output to selected ostream
   **/
  for (unsigned int i = 0; i<(unsigned int)lines_to_print; i++)
    {
	  outFile << std::setw(2) << std::setfill('0') << date_and_time[i][0] << "/"
			  << std::setw(2) << std::setfill('0') << date_and_time[i][1] << "/"
			  << std::setw(2) << std::setfill('0') << date_and_time[i][2] << "\t"
			  << std::setw(2) << std::setfill('0') << date_and_time[i][3] << ":"
			  << std::setw(2) << std::setfill('0') << date_and_time[i][4] << ":"
			  << std::setw(2) << std::setfill('0') << date_and_time[i][5] << "\t";
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
		std::vector< std::vector<double> > &data,
		std::vector< std::vector<int> >    &date_and_time,
		bool day_number_column)
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

void DataTools::read_data(const std::vector< std::string >   &filenames,
		std::vector< std::vector<double> > &data)
{
	data.clear();
	for (unsigned int j=0; j<filenames.size(); j++)
	{
		std::ifstream file;
		open_file(file,filenames[j].c_str());

		std::string line_file;
		std::string file_token;
		while (std::getline(file,line_file))
		{
			std::stringstream file_iss;
			file_iss << line_file;
			std::vector<double> row_data;
			while (std::getline(file_iss,file_token,'\t'))
				row_data.push_back(atof(file_token.c_str()));
			data.push_back(row_data);
		}
		close_file(file);
	}
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
    		  << " x= " << x << std::endl;
      return (INF);
    }
  if (x < table[0].first)
    {
      std::cout << "Warning. In data_tools.h -> interpolate_data. value out of bound." 
    		  << " x= " << x << std::endl;
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

bool DataTools::file_exists(const std::string& name)
{
	/*
	 * Code copied from:
	 * http://stackoverflow.com/questions/12774207/
	 * fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
	 * */
	struct stat buffer;
	return (stat (name.c_str(), &buffer)==0);
}

bool DataTools::dir_exists(const std::string& name)
{
  /*
   * Based on:
   * http://stackoverflow.com/questions/4980815/
   * c-determining-if-directory-not-a-file-exists-in-linux
   */
  struct stat buffer;
  bool directory_exists=false;
  
  if (stat (name.c_str(), &buffer)==0 && S_ISDIR(buffer.st_mode)==true)
    directory_exists=true;
  
  return (directory_exists);
}

time_t DataTools::mktimeUTC(struct tm* timeinfo)
{
  /*
   * Problem: The time seems to be shifted by 1 hour 
   *
   * Answer: saw  in StackOverflow: 
   * http://stackoverflow.com/questions/3660983/c-time-t-problem/3661129
   *
   * "mktime takes a struct tm giving a local time and returns the number
   * of seconds since January 1, 1970, 0:00 UTC. Therefore, your 
   * GetDate(1970,1,1,0,0,0); call will return 0 if your local time zone 
   * is UTC, but may return other values for other time zones.
   *
   * Edit: For a UTC version of mktime or your GetDate, try the following (untested):
   *    Call getenv to save the current value of the TZ environment variable (if any).
   *    Call putenv to change the TZ environment variable to "UTC".
   *    Call tzset to make your changes active.
   *    Call mktime.
   *    Restore the old value of TZ, then call _tzset again."
   */
  // *** enter in UTC mode
  char* oldTZ = getenv("TZ");
  char tzUTC[] = "TZ=UTC";
  putenv(tzUTC);
  tzset();
  // ***
  
  time_t ret = mktime ( timeinfo );

  // *** Restore previous TZ
  if(oldTZ == NULL)
    {
      char tz[] = "TZ=";
      putenv(tz);
    }
  else
    {
      char buff[255];
      sprintf(buff,"TZ=%s",oldTZ);
      putenv(buff);
    }
  tzset();
  // ***
  return ret;
}

void DataTools::date_to_seconds(const std::vector< std::vector<int> > &date_and_time,
				std::vector< std::vector<int> > &date_in_seconds)
{
  time_t rawtime;// = time(0);
  time(&rawtime);
  struct tm *ref_date; // 01/07/2005   00:00:00
  
  // time_t rawtime = 0;
  // struct tm * ref_date; // 01/07/2005   00:00:00
  
  ref_date = localtime(&rawtime);  
  ref_date->tm_hour = 0;
  ref_date->tm_min  = 0;
  ref_date->tm_sec  = 0;
  ref_date->tm_mon  = 7-1;
  ref_date->tm_mday = 1;
  ref_date->tm_year = 2005-1970 + 70;
  time_t ref_seconds = mktimeUTC(ref_date);
  
  // std::cout << ref_date->tm_isdst << " " << ref_date->tm_yday << " " << ref_date->tm_wday << "\t"
  // << ref_date->tm_mday  << "/" << ref_date->tm_mon  << "/" << ref_date->tm_year << "\t"
  // << ref_date->tm_hour  << ":" << ref_date->tm_min  << ":" << ref_date->tm_sec  << "\t" 
  // << ref_seconds  << std::endl;
  
  // Translate date and time entries to seconds since 01/07/2005
  for (unsigned int i=0; i<date_and_time.size(); i++)
    {
      time_t dummy_time = 0;
      // time(&dummy_time);
      struct tm *date;
      
      date = localtime(&dummy_time);
	
      date->tm_mday = date_and_time[i][0];
      date->tm_mon  = date_and_time[i][1]-1;
      date->tm_year = date_and_time[i][2]-1970+70;

      date->tm_hour = date_and_time[i][3];
      date->tm_min  = date_and_time[i][4];
      date->tm_sec  = date_and_time[i][5];
      date->tm_isdst = 1;
            
      std::vector<int> temp1;
      temp1.clear();
      temp1.push_back((int)difftime(mktimeUTC(date),ref_seconds));
      date_in_seconds.push_back(temp1);      
      // if (i<15)
      // std::cout << mktimeUTC(date) << "\t" << ref_seconds << "\t" << temp1[0] << std::endl;
      // if (i<15)
      // std::cout 
      // << date->tm_mday << "/" << date->tm_mon << "/" << date->tm_year << "\t" 
      // << date->tm_hour << ":" << date->tm_min << ":" << date->tm_sec << "\t" 
      // << date_and_time[i][0] << "/" << date_and_time[i][1]-1 << "/" 
      // << date_and_time[i][2]-1900 << "\t"
      // << date_and_time[i][3] << ":" << date_and_time[i][4] << ":" << date_and_time[i][5] << "\t"
      // << mktimeUTC(date) << "\t" << ref_seconds << "\t" << temp1[0] << std::endl;
    }
}

void DataTools::seconds_to_date(std::vector< std::vector<int> > &date_and_time,
				const std::vector< std::vector<int> > &date_in_seconds)
{
  time_t rawtime = 0;
  struct tm * ref_date; // 01/07/2005   00:00:00
  ref_date = localtime(&rawtime);
  ref_date->tm_hour = 0;
  ref_date->tm_min  = 0;
  ref_date->tm_sec  = 0;
  ref_date->tm_mon  = 7-1;
  ref_date->tm_mday = 1;
  ref_date->tm_year = 2005-1970+70;

  time_t ref_seconds = mktimeUTC(ref_date);

  for (unsigned int i=0; i<date_in_seconds.size(); i++)
    {
      time_t time_entry = (time_t)date_in_seconds[i][0] + ref_seconds;
      struct tm * ptm;
      
      ptm = gmtime (&time_entry);
      std::vector<int> temp1;
      temp1.clear();
      temp1.push_back(ptm->tm_mday);
      temp1.push_back(ptm->tm_mon+1);
      temp1.push_back(ptm->tm_year+1900);
      temp1.push_back(ptm->tm_hour);
      temp1.push_back(ptm->tm_min);
      temp1.push_back(ptm->tm_sec);
      
      date_and_time.push_back(temp1);

      // if (i<15)
      // std::cout << ptm->tm_mday << "/" << ptm->tm_mon+1 << "/" << ptm->tm_year+1900 << "\t"
      // << ptm->tm_hour << ":" << ptm->tm_min << ":" << ptm->tm_sec << "\t"
      // << date_in_seconds[i][0] << std::endl;
    }
}

void DataTools::read_met_data(std::vector< std::vector<int> >    &date_and_time,
			      std::vector< std::vector<double> > &met_data,
			      int time_step,
			      const int preheating_step,
			      const std::string met_data_type)
{
	date_and_time.clear();
	met_data.clear();

	bool day_number = false;

	std::string met_data_filename = "met_data_ph_step_";
	std::stringstream ph_step;
	ph_step << preheating_step;
	met_data_filename += ph_step.str() + ".txt";
	/*
	 * Check if we have a file with the previous name. If so
	 * read from this file the met data. We are assuming that
	 * the file contains the corresponding met data for the
	 * preheating step of interest. There is not way to check
	 * if this actually true. I'm just guessing that there
	 * shouldn't be any other file there that happens to have
	 * the same name and doesn't contain the desired met data.
	 */
	std::ifstream file (met_data_filename.c_str());
	if (file.good())
	{
		std::vector<std::string> met_data_file;
		met_data_file.push_back(met_data_filename);

		read_data(met_data_file,
			  met_data,
			  date_and_time,
			  day_number);
	}
	else
	{
		Names names;
		std::vector<std::string> met_data_filenames;

		names.get_met_data_filenames(met_data_filenames,
				preheating_step,
				met_data_type);

		read_data(met_data_filenames,
			  met_data,
			  date_and_time,
			  day_number);

		if (time_step<3600)
		{
			std::vector< std::vector<int> > date_and_time_in_seconds;

			std::vector< std::vector<int> >    temp_date_and_time;
			std::vector< std::vector<int> >    temp_date_and_time_in_seconds;
			std::vector< std::vector<double> > temp_met_data;

			date_to_seconds(date_and_time,
					date_and_time_in_seconds);
			if (date_and_time_in_seconds.size() != date_and_time.size())
				throw 3;

			// generate the output times we want
			int seconds = date_and_time_in_seconds[0][0];
			while (seconds <= date_and_time_in_seconds.back()[0])
			{
				std::vector<int> temp_date_and_time_in_seconds_row;
				temp_date_and_time_in_seconds_row.push_back (seconds);
				seconds += time_step;
				temp_date_and_time_in_seconds.push_back(temp_date_and_time_in_seconds_row);
			}
			// generate the corresponding dates
			seconds_to_date( temp_date_and_time,
					temp_date_and_time_in_seconds);

			std::vector<std::vector<std::pair<double,double> > >tables;
			for (unsigned int j=0; j<met_data[0].size(); j++)
			{
				std::vector<std::pair<double,double> > table;
				for (unsigned int i=0; i<met_data.size(); i++)
					table.push_back(std::make_pair((double)date_and_time_in_seconds[i][0],met_data[i][j]));
				tables.push_back(table);
			}
			// interpolate data according to the previous
			// generated time entries
			for (unsigned int i=0; i<temp_date_and_time_in_seconds.size(); i++)
			{
				std::vector<double> temp1;
				temp1.clear();
				for (unsigned int j=0; j<met_data[0].size(); j++)
					temp1.push_back( interpolate_data(tables[j],(double)temp_date_and_time_in_seconds[i][0]));
				temp_met_data.push_back(temp1);
			}
			met_data = temp_met_data;
			date_and_time = temp_date_and_time;
		}
		/*
		 * If we need to consider several years we do it like this
		 */
		if (preheating_step==1)
		{
			unsigned int one_year_size = date_and_time.size();
			for  (unsigned int i=0; i<7; i++)
				for (unsigned int j=0; j<one_year_size; j++)
				{
					std::vector<int> date;
					date=date_and_time[j];
					date[2]= date[2]+i+1;

					date_and_time.push_back(date);
					met_data.push_back(met_data[j]);
				}
		}
	}
}
