#pragma once

#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <algorithm>

#include "EquationOfState.h"
#include "rootfinder.h"

class Integrator {
public:
	Integrator(EoS* eq, bool verbose = false, bool integrate = true) : eqs(eq), verbose(verbose) { if (integrate) { generateFunction(.001, 100.); }; };
	void writeToFile(std::ofstream* output);
	void generateFunction(double dt, double tmax);
	std::vector<double> getDataLast() { return data.back(); };
	std::vector<std::vector<double>> getData() { return data; };
	double getTLast() { return tvec.back(); };
	std::vector<double> RK4Step(double t, double dt, std::vector<double> state);

protected:
	EoS* eqs;
	void printState(std::vector<double> state);
	std::vector<std::vector<double>> data;
	std::vector<double> tvec;
	bool verbose;
};

void Integrator::writeToFile(std::ofstream* output)
{
	for (int i = 0; i < data.size(); i++) {
		*output << tvec[i] << " ";
		for (int j = 0; j < data[i].size(); j++) {
			*output << data[i][j] << " ";
		}
		*output << "\n";
	}
}

void Integrator::printState(std::vector<double> state)
{
	for (int i = 0; i < state.size(); i++) {
		std::cout << state[i] << ", ";
	};
	std::cout << std::endl;
}

void Integrator::generateFunction(double dt0 = .001, double tmax = 100.) {
	std::vector<double> state = eqs->getInitialState();
	double t = 0;
	double dt = dt0;
	data.clear();
	data.push_back(state);
	tvec.clear();
	tvec.push_back(t);
	
	double dtout = 2 * dt;
	double tout = t - dtout;

	std::vector<double> stateLast(state.size());
	int it = 0;
	while (eqs->stopCondition(state) > 0.0 && t < tmax) {
		std::vector<double> step = RK4Step(t, dt, state);
		double deltaMax = 0.;
		for (int i = 0; i < state.size(); i++) {
			deltaMax = std::max(deltaMax, fabs(step[i]) / std::max(1., fabs(state[i])));
		};
		
		if (deltaMax < 1.e-2) {
			//std::cout << "took step" << std::endl;
			stateLast = state;
			if (t>=tout) {
				data.push_back(state);
				tvec.push_back(t);
				tout += dtout;
			}
			if (verbose) {
				std::cout << "dt=" << dt << ", deltaMax=" << deltaMax << std::endl;
				std::cout << "State: ";
				printState(state);
				std::cout << "Step: ";
				printState(step);
				
			}
			t += dt;
			it += 1;
			for (int i = 0; i < state.size(); ++i) {
				state[i] += step[i];
			}
			if (verbose) std::cout << "StopCondition after step: " << eqs->stopCondition(state) << std::endl;
			if (deltaMax < 1.e-3) { dt = dt * 1.5; if (verbose) std::cout << "Step too small" << std::endl; };
		}
		else {
			dt = dt / 2.;
		};
	};
	if (t >= tmax) { std::cout << "Max iteration count exceeded in generateFunction" << std::endl; };
	auto zeroFunc = [this, stateLast, t](double dt) {//this changes sign when the RK step overruns the stop condition
		auto step = RK4Step(t, dt, stateLast);		  //the root of this is the closest point to the stop condition
		for (int i = 0; i < stateLast.size(); ++i) {
			step[i] += stateLast[i];
		}
		return eqs->stopCondition(step);
	};
	//get as close to the stop condition as possible for the last point
	double dtmin = Find1DRoot(0.0, dt, zeroFunc);
	std::vector<double> step = RK4Step(t, dtmin, stateLast);
	for (int i = 0; i < state.size(); ++i) {
		state[i] = stateLast[i] + step[i];
	}
	t = t + dtmin;
	data.push_back(state);
	tvec.push_back(t);
	//std::cout << eqs->stopCondition(state) << std::endl;
};

std::vector<double> Integrator::RK4Step(double t, double dt, std::vector<double> state) {
	std::vector<double> step(state.size());

	// The first RK4 step 
	std::vector<double> k1 = eqs->getDiffs(t,state);
	for (unsigned int i = 0; i < state.size(); ++i)
		step[i] = state[i] + 0.5 * dt * k1[i];

	// The second RK4 step 
	std::vector<double> k2 = eqs->getDiffs(t + 0.5 * dt, step);
	for (unsigned int i = 0; i < state.size(); ++i)
		step[i] = state[i] + 0.5 * dt * k2[i];

	// The third RK4 step 
	std::vector<double> k3 = eqs->getDiffs(t + 0.5 * dt, step);
	for (unsigned int i = 0; i < state.size(); ++i)
		step[i] = state[i] + dt * k3[i];

	// The third RK4 step 
	std::vector<double> k4 = eqs->getDiffs(t + dt, step);

	// Combine all of the substeps to get fourth order accurate update
	std::vector<double> stepFinal(state.size());
	for (unsigned int i = 0; i < state.size(); ++i)
		stepFinal[i] = dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;

	return stepFinal;
}
