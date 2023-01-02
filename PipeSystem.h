class PipeSystem
{
public:
	PipeSystem(const double time_step_);

	void pipe_heat_flux(const std::vector<double> &new_avg_pipe_temperature,
						const std::vector<double> &old_avg_pipe_temperature,
						std::vector<double> &current_new_inlet__temperatures_pipes,
						std::vector<double> &current_new_outlet_temperatures_pipes,
						std::vector<double> &current_new_pipe_heat_flux,
						bool print_info = false);

private:
	const unsigned int n_pipes;
	const double time_step;
	const double lenght_of_pipe_section;	   //=30.+0.5*2.5*3.14159265359;  //(m)
	const double pipe_external_radius;		   //=0.0125; //(m)
	const double pipe_internal_radius;		   //=0.0102; //(m)
	const double pipe_flow_rate;			   //=0.155;//0.155;        //(kg/s)
	const double water_specific_heat_capacity; //=4181.3;//(J/kgK)
	const double pi;						   //=3.14159265359;

	double pipe_internal_area;
	double pipe_cross_section_area;
	double flow_velocity;
	double time_per_loop;
	double number_of_loops;
	double pipe_wall_heat_coefficient;
};

PipeSystem::PipeSystem(const double time_step_) : n_pipes(40),
												  time_step(time_step_),
												  lenght_of_pipe_section(30. + 0.5 * 2.5 * 3.14159265359),
												  pipe_external_radius(0.0125),
												  pipe_internal_radius(0.0102),
												  pipe_flow_rate(0.155),
												  water_specific_heat_capacity(4181.3),
												  pi(3.14159265359)
{
	/**
	 * Units and approximate values
	 * pipe_internal_area         - (m2)  - 2.17 m2
	 * pipe_cross_section_area    - (m2)  - 3.26E-4 m2
	 * flow_velocity              - (m/s) - 0.4375 m/s      0.4294m/s (0.14 l/s)
	 * time_per_loop              - (s)   - 68.57 s
	 * number_of_loops            - ()    - 900/68.57 ~ 13.125
	 * pipe_wall_heat_coefficient - (W/m2K) - ~192.85 (W/m2K) with k=0.4
	 */
	pipe_internal_area = 2. * pi * pipe_internal_radius * lenght_of_pipe_section;										 //(m2) 2.17 m2
	pipe_cross_section_area = pi * pow(pipe_internal_radius, 2);														 //(m2) 3.26E-4 m2
	flow_velocity = (pipe_flow_rate / 1000.) / pipe_cross_section_area;													 //(m/s) 0.4375 m/s      0.4294m/s (0.14 l/s)
	time_per_loop = lenght_of_pipe_section / flow_velocity;																 //(s) 68.57 s
	number_of_loops = (time_step / time_per_loop);																		 // 900/68.57= 13.125
	pipe_wall_heat_coefficient = 1. / ((pipe_internal_radius * log(pipe_external_radius / pipe_internal_radius)) / 0.4); // with k=0.4  ~192.85 (W/m2K)
}

void PipeSystem::pipe_heat_flux(
	const std::vector<double> &new_avg_pipe_temperature,
	const std::vector<double> &old_avg_pipe_temperature,
	std::vector<double> &current_new_inlet__temperatures_pipes,
	std::vector<double> &current_new_outlet_temperatures_pipes,
	std::vector<double> &current_new_pipe_heat_flux,
	bool print_info)
{
	/*
	 * Sanity check that all vectors match the number of pipes
	 */
	if ((n_pipes != new_avg_pipe_temperature.size()) ||
		(n_pipes != old_avg_pipe_temperature.size()) ||
		(n_pipes != current_new_inlet__temperatures_pipes.size()) ||
		(n_pipes != current_new_outlet_temperatures_pipes.size()) ||
		(n_pipes != current_new_pipe_heat_flux.size()))
	{
		std::cout << "Error, the number of pipes should be "
				  << n_pipes << std::endl
				  << "Error in PipeSystem::pipe_heat_flux"
				  << std::endl;
		throw -1;
	}
	/*
	 * Get average temperatures at collector and storage pipes for current and previous time steps
	 */
	double new_avg_collector_temperature = 0;
	double new_avg_storage___temperature = 0;
	double old_avg_collector_temperature = 0;
	double old_avg_storage___temperature = 0;
	for (unsigned int i = 0; i < n_pipes / 2; i++)
	{
		new_avg_collector_temperature += new_avg_pipe_temperature[i] / (n_pipes / 2);
		old_avg_collector_temperature += old_avg_pipe_temperature[i] / (n_pipes / 2);
		new_avg_storage___temperature += new_avg_pipe_temperature[i + n_pipes / 2] / (n_pipes / 2);
		old_avg_storage___temperature += old_avg_pipe_temperature[i + n_pipes / 2] / (n_pipes / 2);
	}
	/*
	 * Here we populate the inlet and outlet temperatures for collector and storage pipes.
	 * We assume having 40 pipes.
	 * The logic is:
	 * collector inlet pipes (00-09) receive water from the storage header
	 * collector outlet pipes (10-19) receive water from the collector header (!check this!)
	 * storage inlet pipes (20-29) receive water from the collector header
	 * storage outlet pipes (30-39) receive water from the storage header (!check this!)
	 */
	for (unsigned int i = 0; i < n_pipes; i++)
	{
		if (i < n_pipes / 2)
		{
			if (i < n_pipes / 4)
				current_new_inlet__temperatures_pipes[i] = new_avg_storage___temperature;
			else
				current_new_inlet__temperatures_pipes[i] = new_avg_collector_temperature;
		}
		else
		{
			if (i < (n_pipes / 2 + n_pipes / 4))
				current_new_inlet__temperatures_pipes[i] = new_avg_collector_temperature;
			else
				current_new_inlet__temperatures_pipes[i] = new_avg_storage___temperature;
		}
	}
	/*
	 * The water takes appoximately 60 s to flow in each pipe segment, ~240 s for
	 * the whole system. Depending on the timestep choosen, the water will
	 * complete a certain number of cycles. In each cycle the water carries heat
	 * between the collector and storage pipes. This is taken into account in the
	 * next loop.
	 */
	for (unsigned int j = 0; j < (unsigned int)number_of_loops; j++)
	{
		/*
		 * Print content of average pipe temperature vectors for debugging.
		 * Likewise for inlet and outlet temperatures for each loop.
		 */
		if (print_info == true)
		{
			if (j == 0)
			{
				std::cout << "\t\tCollector avg soil: ";
				for (unsigned int i = 0; i < 20; i++)
					std::cout << (0.5 * new_avg_pipe_temperature[i] + 0.5 * old_avg_pipe_temperature[i]) << "\t";
				std::cout << "\n";
				std::cout << "\t\tStorage avg soil  : ";
				for (unsigned int i = 0; i < 20; i++)
					std::cout << (0.5 * new_avg_pipe_temperature[20 + i] + 0.5 * old_avg_pipe_temperature[20 + i]) << "\t";
				std::cout << "\n";
			}

			std::cout << "\tloop: " << j << "\n";
			std::cout << "\t\tCollector inlets : ";
			for (unsigned int i = 0; i < 20; i++)
				std::cout << current_new_inlet__temperatures_pipes[i] << "\t";
			std::cout << "\n";
			std::cout << "\t\tCollector outlets: ";
			for (unsigned int i = 0; i < 20; i++)
				std::cout << current_new_outlet_temperatures_pipes[i] << "\t";
			std::cout << "\n";
			std::cout << "\t\tStorage inlets   : ";
			for (unsigned int i = 0; i < 20; i++)
				std::cout << current_new_inlet__temperatures_pipes[20 + i] << "\t";
			std::cout << "\n";
			std::cout << "\t\tStorage outlets  : ";
			for (unsigned int i = 0; i < 20; i++)
				std::cout << current_new_outlet_temperatures_pipes[20 + i] << "\t";
			std::cout << "\n";
		}

		for (unsigned int i = 0; i < n_pipes; i++)
		{
			/*
			 * Calculate the pipe efficiency, for this we need to calculate
			 * first the overall heat transfer coefficient. This depends on
			 * the pipes convective heat transfer coefficient (which in turn
			 * depends on the temperature of the water) and the transfer
			 * coefficient on the pipe's wall (actually, this term is the
			 * bottle neck).
			 * Obviously, the temperature of the water varies along the pipe.
			 * We set it to the average of the pipe's inlet and the surrounding
			 * soil.
			 */
			ConvectiveCoefficient pipe_convective_coefficient(pipe_external_radius, pipe_flow_rate, 10 /*%*/);
			if (i < n_pipes / 2)
				pipe_convective_coefficient
					.EG_temperature(0.5 * new_avg_collector_temperature +
									0.5 * current_new_inlet__temperatures_pipes[i]);
			else
				pipe_convective_coefficient
					.EG_temperature(0.5 * new_avg_storage___temperature +
									0.5 * current_new_inlet__temperatures_pipes[i]);
			/**
			 * Units
			 * pipe_heat_coefficient         - (W/m2K)
			 * pipe_overall_heat_coefficient - (W/m2K)
			 * collector_UA                  - (W/K)
			 * mC                            - (kg/s)(J/kgK)=(W/K)
			 * pipe_efficiency               - ()
			 * heat_flux_in_current_pipe     - (W)
			 */
			double pipe_heat_coefficient = pipe_convective_coefficient.get_convective_coefficient();
			double pipe_overall_heat_coefficient = 1. / ((1. / pipe_wall_heat_coefficient) + (1. / pipe_heat_coefficient));
			double collector_UA = pipe_overall_heat_coefficient * pipe_internal_area;
			double mC = pipe_flow_rate * water_specific_heat_capacity;
			double pipe_efficiency = 1. - exp(-1. * collector_UA / mC);
			double heat_flux_in_current_pipe = 0.;

			if (i < n_pipes / 2)
				heat_flux_in_current_pipe = pipe_efficiency * mC *
											(1.0 * new_avg_collector_temperature + 0.0 * old_avg_collector_temperature - current_new_inlet__temperatures_pipes[i]);
			else
				heat_flux_in_current_pipe = pipe_efficiency * mC *
											(1.0 * new_avg_storage___temperature + 0.0 * old_avg_storage___temperature - current_new_inlet__temperatures_pipes[i]);

			current_new_pipe_heat_flux[i] += heat_flux_in_current_pipe; // (W)
			current_new_outlet_temperatures_pipes[i] = (heat_flux_in_current_pipe / mC) + current_new_inlet__temperatures_pipes[i]; // (K)
		}
		/*
		 * Re-define the pipe inlet temperature. 0-9 are flow pipes
		 * and 10-19 are return pipes (collector). Likewise, 20-29
		 * are flow pipes and 30-39 are return pipes (storage).
		 * The assignation is
		 * P10_i <- P09_o
		 * P11_i <- P08_o
		 * .
		 * .
		 * .
		 * P19_i <- P00_o
		 * and so on, the same is done for the storage pipes
		 */
		for (unsigned int i = 0; i < 10; i++)
		{
			current_new_inlet__temperatures_pipes[10 + i] = current_new_outlet_temperatures_pipes[9 - i];
			current_new_inlet__temperatures_pipes[30 + i] = current_new_outlet_temperatures_pipes[29 - i];
		}
		/*
		 * The inlet temperatures for P00...P09 and P20...P29 are
		 * defined as the average temperatures from P10 to P19 and
		 * from P30 to P39 respectively
		 */
		double collector_inlet_storage___outlet = 0.;
		double storage___inlet_collector_outlet = 0.;
		for (unsigned int i = 0; i < 10; i++)
		{
			collector_inlet_storage___outlet += current_new_outlet_temperatures_pipes[30 + i];
			storage___inlet_collector_outlet += current_new_outlet_temperatures_pipes[10 + i];
		}
		for (unsigned int i = 0; i < 10; i++)
		{
			current_new_inlet__temperatures_pipes[i] = collector_inlet_storage___outlet / 10.;
			current_new_inlet__temperatures_pipes[20 + i] = storage___inlet_collector_outlet / 10.;
		}
	}
	/*
	 * This is the average heat flux for all loops in the current iteration of
	 * the current timestep
	 */
	for (unsigned int i = 0; i < n_pipes; i++)
	{
		if (i < n_pipes / 2)
			current_new_pipe_heat_flux[i] =
				current_new_pipe_heat_flux[i] / (int)number_of_loops;
		else
			current_new_pipe_heat_flux[i] =
				current_new_pipe_heat_flux[i] / (int)number_of_loops;
	}
}