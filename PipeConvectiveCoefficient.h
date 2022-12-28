#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679
class ConvectiveCoefficient
{
 public:
  ConvectiveCoefficient(double radius, double flow, double EG/*, double T*/);
  /*
   * The convective heat transfer coefficient for the aqueous solution inside the pipes 
   * depends on several things. For starters, the percentage of ethylene glycol in the
   * solution, the solution mean temperature. This defines some physical properties
   * such as: thermal conductivity, dynamic viscosity, specific heat capacity, specific
   * gravity and density. Which in turn are used together with the diameter of the
   * pipes and the mean velocity of the solution to calculate non-dimensional numbers
   * such as Reynolds, Prandtl, and Nusselt as well as the friction factor that define
   * the heat transfer coefficient.
   * The temperature is the more variable parameter in the problem and all the physical
   * properties depend on it. Variations on other parameters such as mean velocity and
   * pipe radius can help to optimize the system performance.
   * 
   * This is the function that give us the convective coefficient, it calls all the other
   * private fuctions that do the actual calculations. This function, the constructor and
   * print_data, that prints all the information to the screen, are the only available
   * from outside this class.
   */					 
  void   print_data();
  void   EG_temperature(double T);
  double get_convective_coefficient();
  
 private:
  /*
   * First we calculate the physical properties of the aqueous ethylene-glycol solution.
   * The following functions do the job. They are accessed from
   * get_convective_coefficient() and modify directly private variables on the class,
   * so there is no need for a return argument.
   */
  void ethylene_glycol_thermal_conductivity  (); 
  void ethylene_glycol_specific_heat_capacity(); 
  void ethylene_glycol_dynamic_viscosity     (); 
  void ethylene_glycol_density               ();
  /*
   * With the physical properties known we can calculate the nondimensional numbers:
   * Reynolds, Prandtl, Nusselt and the friction factor. The function receives a bool 
   * parameter to print to screen the results and check if they are correct.
   */
  void nondimensional_numbers                ();
  /*
   * The first set of private variables are initialized in the constructor with the
   * arguments passed from the call to this class
   */
  double pipe_external_diameter; // meters
  double pipe_internal_diameter; // meters
  double volumetric_flow;        // liters/second
  double ethylene_glycol;        // percentage
  double mean_temperature;       // celsius
  /*
   * The next set of variables are those which are going to be calculated in the object,
   * based on the arguments passed to the constructor. They are initialized to zero in
   * the constructor. This is in case print is called before get_convective_coefficient.
   */
  double h_eg;   // watts    / (square metre - kelvin)
  double k_eg;   // watts    / (metre - kelvin)
  double Cp_eg;  // joules   / (kilogram - kelvin)
  double mu_eg;  // kilogram / (metre - second)
  double rho_eg; // kilogram / (cubic metre)

  double reynolds; // nondimensional
  double prandtl;  // nondimensional
  double nusselt;  // nondimensional
  double f;        // nondimensional
};
/*
 * This is the constructor of the class. It receives four arguments: the external radius
 * of the pipe, the volumetric flow going through the pipes, the percentage of
 * ethylene-glycol present on the solution and its mean temperature.
 */
ConvectiveCoefficient::ConvectiveCoefficient (double radius,
                                              double flow,
                                              double EG/*, double T*/):
  pipe_external_diameter(2.*radius),
  pipe_internal_diameter(2.*radius - 2.*0.0023),
  volumetric_flow(flow),
  ethylene_glycol(EG)
  //mean_temperature(T)
{
  h_eg   = 0.0;
  k_eg   = 0.0;
  Cp_eg  = 0.0;
  mu_eg  = 0.0;
  rho_eg = 0.0;

  reynolds = 0.0;
  prandtl  = 0.0;
  nusselt  = 0.0;
  f        = 0.0;
  
  mean_temperature = 0.0;
}
/* 
 * The calculus of the thermal conductivity for an aqueous ethylene-glycol solution is
 * based on a product guide document from MEGlobal. It is given by:
 * k_eg = A + B*T + C*T*T								 
 * where the coefficients A,B and C depend on the percentage of ethylene glycol present
 * in the solution and T is its temperature. The result is in Btu(F)/hr(ft2)F and needs
 * to be converted to SI units:
 * 1 Btu(F)/hr(ft2)F = 1.7307 W/mK
 */
void ConvectiveCoefficient::ethylene_glycol_thermal_conductivity(){
  double a=0.;
  double b=0.;
  double c=0.;
  double X[3] = { 10.00       ,  20.00        , 30.00       };
  double A[3] = {  0.30433     ,  0.28697     ,  0.27038    };
  double B[3] = {  0.00089727  ,  0.0006635   ,  0.00045096 };
  double C[3] = { -0.0000036114, -0.0000029292, -0.000002316};

  for (unsigned int i = 1; i<3; i++){
      if (ethylene_glycol <= X[i] && ethylene_glycol >= X[i-1]){
	      a = ((ethylene_glycol - X[i-1])/(X[i]-X[i-1]))*(A[i]-A[i-1]) + A[i-1];
	      b = ((ethylene_glycol - X[i-1])/(X[i]-X[i-1]))*(B[i]-B[i-1]) + B[i-1];
	      c = ((ethylene_glycol - X[i-1])/(X[i]-X[i-1]))*(C[i]-C[i-1]) + C[i-1];
	      break;
	    }
    }

  k_eg = 1.7307*(a+b*mean_temperature+c*mean_temperature*mean_temperature);
}
/* 
 * The calculus of the thermal specific heat capacity for an aqueous ethylene-glycol
 * solution is based on a product guide document from MEGlobal. It is given by:
 * Cp_eg = a + b*T + c*T*T								 
 * where the coefficients A,B and C depend on the percentage of ethylene glycol present
 * in the solution and T is its temperature. The result is in Btu/lbF and needs to be
 * converted to SI units:
 * 1 Btu/lbF = 4186.8 J/kgK
*/
void ConvectiveCoefficient::ethylene_glycol_specific_heat_capacity( ){
  double a=0.;
  double b=0.;
  double c=0.;
  double X[3] = { 10.00        , 20.00     , 30.00     };
  double A[3] = { 0.97236      , 0.93576   , 0.89889   };
  double B[3] = { 0.00018001   , 0.00039963, 0.00051554};
  double C[3] = { 0.00000057049, 0.0       , 0.0       };

  for (unsigned int i = 1; i<3; i++){
      if (ethylene_glycol <= X[i] && ethylene_glycol >= X[i-1]){
	      a = ((ethylene_glycol - X[i-1])/(X[i]-X[i-1]))*(A[i]-A[i-1]) + A[i-1];
	      b = ((ethylene_glycol - X[i-1])/(X[i]-X[i-1]))*(B[i]-B[i-1]) + B[i-1];
	      c = ((ethylene_glycol - X[i-1])/(X[i]-X[i-1]))*(C[i]-C[i-1]) + C[i-1];
	      break;
	    }
    }

  Cp_eg = 4186.8*(a+b*mean_temperature+c*mean_temperature*mean_temperature);
}
/* The calculus of the dynamic viscosity for an aqueous ethylene-glycol solution is
 * based on a product guide document from MEGlobal. It is given by:
 *  \mu_eg = 10^(a-b/(x+c))
 * where the coefficients A,B and C depend on the temperature of the solution and x is
 * the percentage of ethylene-glycol present in the solution. The result is in
 * centipoise and needs to be converted to SI units:
 * 1cP = 0.001 kg/ms
 */
void ConvectiveCoefficient::ethylene_glycol_dynamic_viscosity( ){
  double a=0.;
  double b=0.;
  double c=0.;
  double T[3] = {-1.111   , 10.00   , 37.7777 };
  double A[3] = {-3.770236,-4.489869,-3.96839 };
  double B[3] = { 1495.186, 1941.309, 1596.092};
  double C[3] = {-368.93  ,-422.768 ,-420.283 };
  /*
   * Throw and error if the mean temperature is above expected
   */
  if (mean_temperature>T[2])
    {
      /* std::cout << "Mean temperature: " << mean_temperature  */
      /* 		<< " is higher than expected: " << T[2] */
      /* 		<< std::endl; */
      //throw 1000;
      mean_temperature=T[2];
    }  
  if (mean_temperature<T[0])
    mean_temperature=T[0];
  
  for (unsigned int i=1; i<3; i++){
      if (mean_temperature<=T[i] && mean_temperature>=T[i-1]){
	      a = ((mean_temperature - T[i-1])/(T[i]-T[i-1]))*(A[i]-A[i-1]) + A[i-1];
	      b = ((mean_temperature - T[i-1])/(T[i]-T[i-1]))*(B[i]-B[i-1]) + B[i-1];
	      c = ((mean_temperature - T[i-1])/(T[i]-T[i-1]))*(C[i]-C[i-1]) + C[i-1];
	      break;
	    }
    }
  
  mu_eg = 0.001*pow(10,a-b/(ethylene_glycol+c));
}
/* The calculus of the specific gravity for an aqueous ethylene-glycol solution is based
 * on a product guide document from MEGlobal. It is given by:
 * \rho_eg = a+b*x+c*T*T								 
 * where the coefficients A,B and C depend on the temperature of the solution and x is
 * the percentage of ethylene-glycol present in the solution.
 * This give us the specific gravity that is the ratio of the density of the solution
 * to the density of water at 60F.
*/
void ConvectiveCoefficient::ethylene_glycol_density( ){
  double a=0.;
  double b=0.;
  double c=0.;
  double T[3] = {-1.111       , 10.00        , 37.7777      };
  double A[3] = { 0.98147     ,  0.99873     ,  0.99284     };
  double B[3] = { 0.002498    ,  0.0016424   ,  0.0014017   };
  double C[3] = {-0.0000091168, -0.0000040019, -0.0000029868};
  /*
   * Throw and error if the mean temperature is above expected
   */
  if (mean_temperature>T[2]){
      // std::cout << "Mean temperature: " << mean_temperature 
      // 		<< " is higher than expected: " << T[2]
      // 		<< std::endl;
      //throw 1000;
      mean_temperature=T[2];
    }
  if (mean_temperature<T[0])
    mean_temperature=T[0];

  for (unsigned int i = 1; i<3; i++){
      if (mean_temperature <= T[i] && mean_temperature >= T[i-1]){
	      a = ((mean_temperature - T[i-1])/(T[i]-T[i-1]))*(A[i]-A[i-1]) + A[i-1];
	      b = ((mean_temperature - T[i-1])/(T[i]-T[i-1]))*(B[i]-B[i-1]) + B[i-1];
	      c = ((mean_temperature - T[i-1])/(T[i]-T[i-1]))*(C[i]-C[i-1]) + C[i-1];
	      break;
	    }
    }
  /*
   * The specific gravity is multiplied by the density of water at 60F (15.55C) in kg/m3	*/
  rho_eg = 998.8025*(a+b*ethylene_glycol+c*pow(ethylene_glycol,2)); }
/*
 * This function does the calculation of nondimensional numbers: Reynolds
 * Prandtl, Nusselt and resistance factor. These numbers may require 
 * different equations based on the physical properties of the solution
 * and geometrical variables of the pipe that defines the regime of the
 * fluid (laminar, turbulent or transitional).
 * Depending on the Reynolds number, the fluid can be laminar, transitional
 * or turbulent. For the laminar and turbulent regimes there are defined
 * relations for the resistance factor but not for the transitional regime.
 * According to Nikuradse (1933) [saw in Engineering Fluid Mechanics by 
 * John A. Roberson] the behaviour in this region is linear and so what I
 * do here is to calculate a linear expresion for the region based on 
 * the values of the resistance coefficient on the extremes.
 * It is important to mention the criterion used to define where the fluid
 * is laminar, transitional or turbulent. For example, Cengel (Heat Transfer:
 * A practical approach) says that the transitional region is 2300<Re<10000
 * and that lower values assure laminar behaviours and greater values assure
 * turbulent ones. However, Roberson mentions that altough it is possible
 * to maintain laminar regimes in high Reynolds values, these regimes are
 * unstable and subjected to become turbulent when there are vibrations
 * present. In most engineering applications there are vibrations involved
 * (as is the case of a highway) and it is expected that the transitional
 * region to be narrower. For this reason I choose the limits offered by
 * Roberson, he defines the transitional region as 2000<Re<4000 for smooth
 * pipes.
 */
void ConvectiveCoefficient::nondimensional_numbers( ){
  //reynolds = rho_eg*volumetric_flow*(pipe_internal_diameter)/
  //(1000*mu_eg*numbers::PI*pow(pipe_internal_diameter/2,2));
  reynolds = rho_eg*volumetric_flow*(pipe_internal_diameter)/
      (1000*mu_eg*PI*pow(pipe_internal_diameter/2,2));
  prandtl  = mu_eg * Cp_eg / k_eg;

  if (reynolds < 2000.) // laminar region
    f = 64./reynolds;
  if (reynolds >= 2000. && reynolds <= 4000.) // transitional region
    f = (((64./2000)-(pow(0.79*log(4000)-1.64,-2)))/(2000.-4000.))*
        (reynolds-2000.) + (64./2000);
  if (reynolds >  4000. && reynolds <= 100000.) // turbulent region (up to Re=100000)
    f = pow(0.79*log(reynolds)-1.64,-2);
  
  nusselt = (f/8.)*(reynolds-1000.)*prandtl/
      (1.+12.7*pow(f/8.,0.5)*(pow(prandtl,2./3.)-1.)); 
}

void ConvectiveCoefficient::print_data(){
  std::cout << std::endl
        << "Data for convective heat transfer coefficient in pipes:" << std::endl
        << "\t internal diameter: " << pipe_internal_diameter << " m"   << std::endl
        << "\t volumetric flow: "   << volumetric_flow        << " l/s" << std::endl
        << "\t mean temperature: "  << mean_temperature       << " C"   << std::endl
        << "\t ethyle-glycol: "     << ethylene_glycol        << " %"   << std::endl
        << "Physical parameters of the aqueous ehtylene-glycol solution" << std::endl
        << "\t k_eg  : " << k_eg   << " W/mK"  << std::endl
        << "\t Cp_eg : " << Cp_eg  << " J/kgK" << std::endl
        << "\t mu_eg : " << mu_eg  << " kg/ms" << std::endl
        << "\t rho_eg: " << rho_eg << " kg/m3" << std::endl
        << "Nondimensional numbers" << std::endl
        << "\t Re: " << reynolds << std::endl
        << "\t Pr: " << prandtl  << std::endl
        << "\t f : " << f        << std::endl
        << "\t Nu: " << nusselt  << std::endl
        << "Convective heat transfer coefficient" << std::endl
        << "\t h_eg: " << h_eg << " W/m2K" << std::endl
        << std::endl;
}

double ConvectiveCoefficient::get_convective_coefficient(){
  ethylene_glycol_thermal_conductivity  (); 
  ethylene_glycol_specific_heat_capacity(); 
  ethylene_glycol_dynamic_viscosity     (); 
  ethylene_glycol_density               (); 

  nondimensional_numbers                ();
  
  h_eg = nusselt * k_eg / pipe_internal_diameter;
  
  return (h_eg);
}

void ConvectiveCoefficient::EG_temperature(double T){
  mean_temperature = T;
}
