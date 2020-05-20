#include <iostream> 
#include <fstream> 
#include <stdexcept>
#include <vector> 
#include <math.h>
#include <algorithm>
#include <ctime>

#include "integrator.h"
#include "rootfinder.h"
#include "EquationOfState.h"

class GammaEquations : public EoS {
public:
	GammaEquations(double a, double gamma, double v0, bool critStop = false) : a(a), gamma(gamma), v0(v0), rmax(100.), critStop(critStop) {}
	int size() { return 3; }
	std::vector<double> getInitialState() { return { 0.0,log(v0),0.0 }; }

	std::vector<double> getDiffs(double t, std::vector<double> state) {
		if (critStop) {
			return {
				1. - exp(2. * state[1] - state[2]),
				a * exp(-state[0] - state[2]) - 2.,
				-(gamma - 1) * (2. * (1. - exp(2. * state[1] - state[2])) + (a * exp(-state[0] - state[2]) - 2.))
			};
		};
		return { 
			fabs(1. - exp(2 * state[1]-state[2])),
			fabs(a * exp(-state[0]-state[2]) - 2.), 
			-(gamma - 1)*(2.*(1. - exp(2. * state[1] - state[2])) + (a * exp(-state[0] - state[2]) - 2.))
		};
	}

	double stopCondition(std::vector<double> state) {
		if (critStop) { return std::min(1. - exp(2. * state[1]-state[2]), a * exp(-state[0]-state[2]) - 2.); };
		return rmax - exp(state[0]);
	}

	double ZeroFunction(double t, std::vector<double> state) {//return f1-f2, used to search for critical point
		return 1. - exp(2. * state[1] - state[2]) - a * exp(-state[0] - state[2]) + 2.;
	}

protected:
	double a;
	double gamma;
	double v0;
	double rmax;
	bool critStop;
};

int main() {
	double a = 10.;
	double g = 1.01;
	std::string fileName = "gamma.txt";
	std::ofstream output(fileName);
	for (double v0 = .001; v0 < .01; v0 += .001) {
		GammaEquations eqs(a, g, v0);
		Integrator solver(&eqs);
		solver.writeToFile(&output);
	};
	//output.close();

	auto zeroFunc = [a,g](double v0) {
		GammaEquations eqs(a, g, v0, true);
		Integrator solver(&eqs);
		double z = eqs.ZeroFunction(solver.getTLast(), solver.getDataLast());
		return z;
	};
	clock_t timereq = clock();
	double v0crit = Find1DRoot(1.e-10, .9, zeroFunc, 1.e-11);
	timereq = clock() - timereq;
	std::cout << "v0=" << v0crit << ", time required: " << (float)timereq / CLOCKS_PER_SEC << std::endl;
	GammaEquations eqs(a, g, v0crit, true);
	Integrator solver(&eqs);
	std::cout << eqs.ZeroFunction(solver.getTLast(), solver.getDataLast()) << std::endl;
	GammaEquations eqs2(a, g, v0crit);
	Integrator solver2(&eqs2);
	solver2.writeToFile(&output);
	output.close();
}

