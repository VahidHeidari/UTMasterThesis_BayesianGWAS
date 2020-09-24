#include <cctype>
#include <iostream>

#include "logger.h"
#include "population.h"
#include "sim-pars.h"

using namespace std;



std::string config_path = DEF_CONFIG_PATH;
SimPars sim_pars;
Population population;
Population admix_population;

void RunSimulation(int iters, Population& population)
{
	for (int i = 0; i < iters; ++i) {					// Run simulation.
		population.StepSimulation();
		population.PrintNumAdmixIndivs();

		if (i == 0 || i + 1 % 100 == 0 || i + 1 == iters)		// Print simulation progress.
			logger << "Simulation step #" << i + 1 << " of " << iters << std::endl;
	}
}

void InitAndRunAdmixture()
{
	logger << "Run admixture simulation . . ." << std::endl;
	admix_population.InitAdmix(sim_pars, population);		// Admix populations.
	admix_population.PrintAdmixProportions("Admix Init");
	admix_population.PrintNumAdmixIndivs();
	RunSimulation(sim_pars.admix_iters, admix_population);
}

int main(int argc, char** argv)
{
	// Log start time.
	logger << endl << endl << endl;
	logger << "-------- DATAGEN START ----------" << endl;
	logger << Time <<std::endl;
	logger << "---------------------------------" << endl;

	if (argc > 1)
		config_path = argv[1];								// Get user config file path, instead of default path.

	if (!sim_pars.Init(config_path))						// Read simulation parameters.
		return 1;

	sim_pars.LogInputParameters();							// Log simulation parameters.
	population.Init(sim_pars);								// Initialize random popuation.

	RunSimulation(sim_pars.iters, population);				// Run simulation.
	if (sim_pars.IsSimModeAdmixture())
		InitAndRunAdmixture();								// Run admixture simulation.

	logger << " ----- RESULTS OF SIMULATIONS -----" << std::endl;
	//population.PrintIndividuals();						// Log last generation individuals.
	population.PrintClustersFrequencies(sim_pars);			// Log last generation allele frequencies.

	if (sim_pars.IsSimModeCaseControl()) {
		logger << "Infect disease . . ." << std::endl;
		population.InfectDisease(sim_pars);
		population.CreateExhMotahariInput(sim_pars);
		population.CreateBEAM3Input();
	} else if (sim_pars.IsSimModeAdmixture()) {
		population.PrintAdmixProportions("No Admix End");
		logger << std::endl << std::endl;
		logger << " ----- RESULTS OF ADMIXTURE SIMULATION -----" << std::endl;
		admix_population.PrintIndividuals();
		admix_population.PrintClustersFrequencies(sim_pars);
		admix_population.PrintAdmixProportions("Result");
		admix_population.PrintNumAdmixIndivs();
		admix_population.CreateSTRUCTUREInput();
		admix_population.CreateMySTRUCTUREAdmixInput();
	} else {
		logger << "Create Structure input configs . . ." << std::endl;
		population.CreateSTRUCTUREInput();
		population.CreateMySTRUCTUREInput();
	}

	// Log end time.
	logger << "--------- DATAGEN END -----------" << endl;
	logger << Time << endl;
	logger << "---------------------------------" << endl;
	return 0;
}
