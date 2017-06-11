#ifndef __DATA_TOOLS_INCLUDED__
#define __DATA_TOOLS_INCLUDED__

#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <stdlib.h> // atoi
#include <algorithm> // lower_bound
#include <iostream> // cout
#include <iomanip> // setw
#include <math.h> // fabs
#include <sys/stat.h> // stat
#include <ctime> // time_t, tm
#include <Names.h>

class DataTools
{
public:
  DataTools();
  void open_file(std::ofstream &file,const std::string filename);
  void open_file(std::ifstream &file,const std::string filename);
  void close_file(std::ofstream &file);
  void close_file(std::ifstream &file);

  void print_data (std::ostream &outFile,
		   const std::vector< std::vector<double> > &data,
		   int lines_to_print=-1);

  void print_data (std::ostream &outFile,
		   const std::vector< std::vector<double> > &data,
		   std::vector< std::vector<int> > &date_and_time,
		   int lines_to_print=-1);

  void read_data(const std::vector< std::string >   &filenames,
		 std::vector< std::vector<double> > &data,
		 std::vector< std::vector<int> >    &date_and_time,
		 bool day_number_column=false);

  void read_data(const std::vector< std::string >   &filenames,
		 std::vector< std::vector<double> > &data);

  void read_data(const std::string &filename,
		  std::vector< std::vector<double> > &data);

  double interpolate_data(std::vector< std::pair<double,double> > table,
			  const double x);

  bool file_exists (const std::string& name);

  bool dir_exists (const std::string& name);
  
  time_t mktimeUTC(struct tm* timeinfo);

  void date_to_seconds(const std::vector< std::vector<int> > &date_and_time,
		       std::vector< std::vector<int> > &date_in_seconds);

  void seconds_to_date(std::vector< std::vector<int> > &date_and_time,
		       const std::vector< std::vector<int> > &date_in_seconds);

  void read_met_data(std::vector< std::vector<int> >    &date_and_time,
		     std::vector< std::vector<double> > &met_data,
		     int time_step,
		     const int preheating_step,
		     const std::string met_data_type,
		     const std::string input_path);
private:  
};

// void make_scatter_plots(const std::string script_filename,
// 			const std::string data_output_filename,
// 			const std::vector< std::string > graph_filenames,
// 			const std::vector< double > &xdata,
// 			const std::vector< double > &ydata,
// 			const std::vector< std::vector<int> > &date_and_time,
// 			std::string main_variable_name,
// 			std::string main_variable_units,
// 			std::string optional_variable_name = "",
// 			std::string optional_variable_units = "",
// 			double optional_variable_value = 0.)
// {
//   if(xdata.size()!=ydata.size())
//     throw 4;
//   /*
//     Generate GNUplot scripts to plot the actual
//     scatter graphs.
//     First we calculate the coefficients for the 
//     least squares line adjusment. The idea is to
//     include in the scatter graph a line of the form
//     Y = Ax + B
//     The following code is to calculate A and B.
//   */
//   double Xmax  = 0;
//   double Xmin  = 10000;
//   double Ymax  = 0;
//   double Ymin  = 10000;
//   double Xmean = 0;
//   double Ymean = 0;
//   double XYmean = 0;
//   double X2mean = 0;
//   double Y2mean = 0;
    
//   double MSE   = 0;
//   double RMSD  = 0;
//   double NRMSD = 0;
  
//   for (unsigned int i=0; i<xdata.size(); i++)
//     {
//       double x = xdata[i];
//       double y = ydata[i];
      
//       if (x > Xmax)
//   	Xmax = x;
//       if (x < Xmin)
//   	Xmin = x;
//       if (y > Ymax)
//   	Ymax = y;
//       if (y < Ymin)
//   	Ymin = y;
      
//       Xmean  += x;
//       Ymean  += y;
//       XYmean += x*y;
//       X2mean += x*x;
//       Y2mean += y*y;
	
//       MSE += pow(x-y,2); 
//     }
  
//   int n = xdata.size();

//   Xmean = Xmean/n;
//   Ymean = Ymean/n;
//   XYmean = XYmean/n;
//   X2mean = X2mean/n;
//   Y2mean = Y2mean/n;
  
//   MSE   = MSE/n;
//   RMSD  = sqrt(MSE);
//   NRMSD = RMSD/(Xmax-Xmin);

//   double alpha  = 0;
//   double beta   = 0;
//   double gamma  = 0;
//   double beta_1 = 0;
//   double beta_2 = 0;
//   double SS_res = 0;
//   double SS_tot = 0;
  
//   // Fitting the regression line  Y = alpha + beta*X  AND  Y = gamma*X
//   for (unsigned int i=0; i<xdata.size(); i++)
//     {
//       double x = xdata[i];
//       double y = ydata[i];
	
//       beta_1 += (x-Xmean)*(y-Ymean);
//       beta_2 += pow(x-Xmean,2);
//     }
//   beta  = beta_1/beta_2;
//   alpha = Ymean - beta*Xmean;
//   gamma = XYmean/X2mean;

//   double Fmax = 0;
//   double Fmin = 0;
//   for (unsigned int i=0; i<xdata.size(); i++)
//     {
//       double x = xdata[i];
//       double y = ydata[i];
//       double f = alpha + beta*x;

//       if (f > Fmax)
//   	Fmax = f;
//       if (f < Fmin)
//   	Fmin = f;

//       SS_tot += pow(y-Ymean,2); //the total sum of squares (proportional to the sample variance)      
//       SS_res += pow(y-f,2);     //the sum of squares of residuals, also called the residual sum of squares.
//     }
//   double R2 = 0;
//   R2 = 1-SS_res/SS_tot;
//   /*
//     Generate files with data to draw scatter plots    
//     date     time     file1_entry     file2_entry     difference     normalized difference     day     month     year
//   */
//   std::vector< std::vector<double> > processed_data;
//   processed_data.clear();
//   for (unsigned int i=0; i<xdata.size(); i++)
//     {
//       std::vector<double> temp;
//       temp.clear();
//       temp.push_back(xdata[i]);
//       temp.push_back(ydata[i]);
//       temp.push_back(ydata[i] - (alpha + beta*xdata[i]));
//       temp.push_back(pow(ydata[i] - (alpha + beta*xdata[i]),2));
//       temp.push_back(pow(ydata[i]-Ymean,2));
//       temp.push_back((double) date_and_time[i][0]);
//       temp.push_back((double) date_and_time[i][1]);
//       temp.push_back((double) date_and_time[i][2]);

//       processed_data.push_back(temp);
//     }
  
//   std::ofstream data_output_file (data_output_filename.c_str());
//   if (!data_output_file.is_open())
//     throw 2;
//   print_data (0,
// 	      date_and_time,
// 	      processed_data,
// 	      /*std::cout*/data_output_file);
//   data_output_file.close();
//   if (data_output_file.is_open())
//     throw 3;

//   // Now, let us write the actual script for gnuplot.
//   // We make use of AWK program. The idea is to plot
//   // with a different color the data for corresponding
//   // to different months

//   std::ofstream outFile (script_filename.c_str());
//   if (!outFile.is_open())
//     throw 2;
  
//   outFile << "#*********************************************************"
//   	  << std::endl
//   	  << " \t #name of output file: "  << data_output_filename << std::endl
//   	  << std::endl
//   	  << "\t #Max X value  : " << Xmax  << std::endl
//   	  << "\t #Min X value  : " << Xmin  << std::endl
//   	  << "\t #Mean X value : " << Xmean << std::endl
//   	  << "\t #Max Y value    : " << Ymax  << std::endl
//   	  << "\t #Min Y value    : " << Ymin  << std::endl
//   	  << "\t #Mean Y value   : " << Ymean << std::endl
//   	  << "\t #MSE   : " << MSE   << std::endl
//   	  << "\t #RMSD  : " << RMSD  << std::endl
//   	  << "\t #NRMSD : " << NRMSD << std::endl
//   	  << "\t #R2    : " << R2    << std::endl
//   	  << std::endl
//   	  << "\t #Equation: \t" << "Xvalue = (" << alpha << ") + (" << beta << ")*Yvalue" 
//   	  << std::endl
//   	  << "#*********************************************************"
//   	  << std::endl;
   
//   outFile << "a =" << alpha << std::endl
//   	  << "b =" << beta  << std::endl
//   	  << "c =" << gamma << std::endl
//   	  << "set title " << "\"" << main_variable_name << " -  Experimental vs Predicted";
  
//   if (optional_variable_value != 0.)
//     outFile << "\\n" << optional_variable_name << " " << optional_variable_value << " " << optional_variable_units;
//   outFile << "\"" << std::endl
//   	  << "set xlabel \"Experimental data " << main_variable_units << "\"" << std::endl
//   	  << "set ylabel \"Predicted data "   << main_variable_units << "\"" << std::endl
//   	  << "set size square"   << std::endl;
  
//   double range_max = 0;
//   double range_min = 0;
//   double tics = 0;
//   double Max = 0;
//   double Min = 0;

//   if (Xmax>Ymax)
//     Max=Xmax;
//   else
//     Max=Ymax;

//   if (Xmin<Ymin)
//     Min=Xmin;
//   else
//     Min=Ymin;

//   if ((Max-Min)<10)
//     {
//       range_max = ceil  (Max);
//       range_min = floor (Min);
//       tics = 0.5;
//     }
//   else if ((Max-Min)<100)
//     {
//       range_max = 10.*ceil  (Max/10);
//       range_min = 10.*floor (Min/10);
//       tics = 2;
//     }
//   else if ((Max-Min)<500)
//     {
//       range_max = 10.*ceil  (Max/10);
//       range_min = 10.*floor (Min/10);
//       tics = 50;
//     }
//   else 
//     {
//       range_max = 100.*ceil  (Max/100);
//       range_min = 100.*floor (Min/100);
//       tics = 100;
//     }
//   /*
//     Scatter plot
//   */  
//   outFile << "set terminal pngcairo dashed size 1700,1000 font \"Arial, 18\"" << std::endl
//   	  << "set output \"" << graph_filenames[0] << "\"" << std::endl;
//   outFile << "set xrange [" << range_min << ":" << range_max << "]" << std::endl
//   	  << "set yrange [" << range_min << ":" << range_max << "]" << std::endl
//   	  << "set xtics " << tics << std::endl
//   	  << "set ytics " << tics << std::endl
//   	  << "set grid"           << std::endl
//   	  << "set key vert at " << range_max + tics/5 << "," << range_max << " left Left" << std::endl
//   	  << "set label 1 \"********************\\n RMSD = " << RMSD << main_variable_units << "\\n NRMSD = " << NRMSD << "\\n R2    = " << R2 
// 	  << " \\n********************\" at " << range_max + tics/5 << "," << 0.5*(range_max+range_min) - tics << std::endl
// 	  << "set lmargin 0" << std::endl
//   	  << std::endl;
//   outFile << "plot\\" << std::endl
//   	  << "\"< awk '{if(($9 ==  \\\"9.000\\\")) print}' " << data_output_filename << "\" u 3:4 with points lt  9 pt  9 t \"September\",\\" << std::endl
//   	  << "\"< awk '{if(($9 == \\\"10.000\\\")) print}' " << data_output_filename << "\" u 3:4 with points lt 10 pt 10 t \"October  \",\\" << std::endl
//   	  << "\"< awk '{if(($9 == \\\"11.000\\\")) print}' " << data_output_filename << "\" u 3:4 with points lt 11 pt 11 t \"November \",\\" << std::endl
//   	  << "\"< awk '{if(($9 == \\\"12.000\\\")) print}' " << data_output_filename << "\" u 3:4 with points lt 12 pt 12 t \"December \",\\" << std::endl
//   	  << "\"< awk '{if(($9 ==  \\\"1.000\\\")) print}' " << data_output_filename << "\" u 3:4 with points lt 13 pt 13 t \"January  \",\\" << std::endl
//   	  << "\"< awk '{if(($9 ==  \\\"2.000\\\")) print}' " << data_output_filename << "\" u 3:4 with points lt 14 pt 14 t \"February \",\\" << std::endl
//   	  << "\"< awk '{if(($9 ==  \\\"3.000\\\")) print}' " << data_output_filename << "\" u 3:4 with points lt 15 pt 15 t \"March    \",\\" << std::endl
//   	  << "\"< awk '{if(($9 ==  \\\"4.000\\\")) print}' " << data_output_filename << "\" u 3:4 with points lt 16 pt 16 t \"April    \",\\" << std::endl
//   	  << "\"< awk '{if(($9 ==  \\\"5.000\\\")) print}' " << data_output_filename << "\" u 3:4 with points lt 17 pt 17 t \"May      \",\\" << std::endl
//   	  << "\"< awk '{if(($9 ==  \\\"6.000\\\")) print}' " << data_output_filename << "\" u 3:4 with points lt 18 pt 18 t \"June     \",\\" << std::endl
//   	  << "\"< awk '{if(($9 ==  \\\"7.000\\\")) print}' " << data_output_filename << "\" u 3:4 with points lt 19 pt 19 t \"July     \",\\" << std::endl
// 	  << "\"< awk '{if(($9 ==  \\\"8.000\\\")) print}' " << data_output_filename << "\" u 3:4 with points lt 20 pt 20 t \"August   \",\\" << std::endl
//       	  << "x"       << " lt 2 lc 4 lw 1 t \"Ideal\",\\" << std::endl
//   	  << "a + b*x" << " lt 1 lc 7 lw 2 t \"Y = " << alpha << "\\n    + " << beta << "*X\"" 
//   	  << std::endl
//   	  << std::endl;
  
//   /*
//     Residual Plot
//   */
//   if (Fmax>Ymax)
//     Max=Fmax;
//   else
//     Max=Ymax;

//   if (Fmin<Ymin)
//     Min=Fmin;
//   else
//     Min=Ymin;
  
//   if ((Max-Min)<10)
//     {
//       range_max = ceil  (Max);
//       range_min = floor (Min);
//       tics = 0.5;
//     }
//   else if ((Max-Min)<100)
//     {
//       range_max = 10.*ceil  (Max/10);
//       range_min = 10.*floor (Min/10);
//       tics = 2;
//     }
//   else if ((Max-Min)<500)
//     {
//       range_max = 100.*ceil  (Max/100);
//       range_min = 100.*floor (Min/100);
//       tics = 50;
//     }
//   else 
//     {
//       range_max = 100.*ceil  (Max/100);
//       range_min = 100.*floor (Min/100);
//       tics = 100;
//     }

//   outFile << "set terminal pngcairo size 1300,1000 font \"Arial, 18\"" << std::endl
//   	  << "set output \"" << graph_filenames[1] << "\"" << std::endl;

//   outFile << "set title " << "\"" << main_variable_name << " - Residuals";
//   if (optional_variable_value != 0.)
//     outFile << "\\n" << optional_variable_name << " " << optional_variable_value << " " << optional_variable_units;
//   outFile << "\"" << std::endl;
  
//   outFile << "set xrange [" << range_min << ":" << range_max << "]" << std::endl
//   	  << "set yrange [" << range_min << ":" << range_max << "]" << std::endl
// 	  << "unset label 1" << std::endl
//   	  << "set key vert at " << range_max + tics/5 << "," << range_max << " left Left" << std::endl
//   	  << "set grid"           << std::endl
// 	  << "set xlabel \"Predicted data " << main_variable_units << "\"" << std::endl
//   	  << "set ylabel \"Residuals "       << main_variable_units << "\"" << std::endl
	  
// 	  << "set lmargin 0" << std::endl
// 	  << std::endl; 
  
//   outFile << "plot\\" << std::endl
//   	  << "\"< awk '{if(($9 ==  \\\"9.000\\\")) print}' " << data_output_filename << "\" u 4:5 with points lt  9 pt  9 t \"September\",\\" << std::endl
//   	  << "\"< awk '{if(($9 == \\\"10.000\\\")) print}' " << data_output_filename << "\" u 4:5 with points lt 10 pt 10 t \"October  \",\\" << std::endl
//   	  << "\"< awk '{if(($9 == \\\"11.000\\\")) print}' " << data_output_filename << "\" u 4:5 with points lt 11 pt 11 t \"November \",\\" << std::endl
//   	  << "\"< awk '{if(($9 == \\\"12.000\\\")) print}' " << data_output_filename << "\" u 4:5 with points lt 12 pt 12 t \"December \",\\" << std::endl
//   	  << "\"< awk '{if(($9 ==  \\\"1.000\\\")) print}' " << data_output_filename << "\" u 4:5 with points lt 13 pt 13 t \"January  \",\\" << std::endl
//   	  << "\"< awk '{if(($9 ==  \\\"2.000\\\")) print}' " << data_output_filename << "\" u 4:5 with points lt 14 pt 14 t \"February \",\\" << std::endl
//   	  << "\"< awk '{if(($9 ==  \\\"3.000\\\")) print}' " << data_output_filename << "\" u 4:5 with points lt 15 pt 15 t \"March    \",\\" << std::endl
//   	  << "\"< awk '{if(($9 ==  \\\"4.000\\\")) print}' " << data_output_filename << "\" u 4:5 with points lt 16 pt 16 t \"April    \",\\" << std::endl
//   	  << "\"< awk '{if(($9 ==  \\\"5.000\\\")) print}' " << data_output_filename << "\" u 4:5 with points lt 17 pt 17 t \"May      \",\\" << std::endl
//   	  << "\"< awk '{if(($9 ==  \\\"6.000\\\")) print}' " << data_output_filename << "\" u 4:5 with points lt 18 pt 18 t \"June     \",\\" << std::endl
//   	  << "\"< awk '{if(($9 ==  \\\"7.000\\\")) print}' " << data_output_filename << "\" u 4:5 with points lt 19 pt 19 t \"July     \",\\" << std::endl
// 	  << "\"< awk '{if(($9 ==  \\\"8.000\\\")) print}' " << data_output_filename << "\" u 4:5 with points lt 20 pt 20 t \"August   \"\\" << std::endl	
//   	  << std::endl
//   	  << std::endl;
//   /*
//     Plot time series
//   */
//   outFile << "set terminal pngcairo size 1700,1000 font \"Arial, 18\"" << std::endl
//   	  << "set output \"" << graph_filenames[2] << "\"" << std::endl;
  
//   outFile << "set title " << "\"" << main_variable_name << " - Time series";
//   if (optional_variable_value != 0.)
//     outFile << "\\n" << optional_variable_name << " " << optional_variable_value << " " << optional_variable_units;
//   outFile << "\"" << std::endl;

//   outFile << "set xdata time" << std::endl
// 	  << "set timefmt \"%d/%m/%Y\\t%H:%M:%S\"" << std::endl
// 	  << "set format x \"%d/%m/%y\"" << std::endl
// 	  << "set xtics auto" << std::endl
// 	  << "set size nosquare"   << std::endl;
  
//   outFile << "set autoscale" << std::endl
// 	  << "set yrange [" << range_min << ":" << range_max << "]" << std::endl
// 	  << "set key default" << std::endl
//   	  << "set grid"           << std::endl
// 	  << "set xlabel \"Time\"" << std::endl
//   	  << "set ylabel \"" << main_variable_units << "\"" << std::endl
// 	  << "set lmargin" << std::endl
// 	  << std::endl;
  
//   outFile << "plot\\" << std::endl
// 	  << "\"" << data_output_filename << "\" u 1:4 w lp lt 25 t \"Predicted   \",\\" << std::endl
// 	  << "\"" << data_output_filename << "\" u 1:3 w lp lt 1  t \"Experimental\"\\" << std::endl;

//   outFile.close();
//   if (outFile.is_open())
//     throw 3;


// }
#endif
