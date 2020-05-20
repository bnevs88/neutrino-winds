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

class IsothermalEquations : public EoS {
public:
	IsothermalEquations(double a, double v0, bool critStop = false) : a(a), v0(v0), rmax(100.), critStop(critStop) {}
	int size() { return 2; }
	std::vector<double> getInitialState() { return { 0.0,log(v0) }; }
	
	std::vector<double> getDiffs(double t, std::vector<double> state) {
		if (critStop) return { 1. - exp(2. * state[1]), a* exp(-state[0]) - 2. };
		return { fabs(1. - exp(2 * state[1])),fabs(a * exp(-state[0]) - 2.) };
	}

	double stopCondition(std::vector<double> state) {
		if (critStop) { return std::min(1. - exp(2. * state[1]), a * exp(-state[0]) - 2.); };
		return rmax - exp(state[0]);
	}

	double ZeroFunction(double t, std::vector<double> state) {//return f1-f2, used to search for critical point
		return 1. - exp(2. * state[1]) - a * exp(-state[0]) + 2.;
	}	

protected:
	double a;
	double v0;
	double rmax;
	bool critStop;
};

int main() {
	double a = 10.;
	std::string fileName = "isothermal.txt";
	std::ofstream output(fileName);
	for (double v0 = .001; v0 < .01; v0+=.001) {
		IsothermalEquations eqs(a, v0);
		Integrator solver(&eqs);
		solver.writeToFile(&output);
	};
	//output.close();

	auto zeroFunc = [a](double v0) {
		IsothermalEquations eqs(a, v0, true);
		Integrator solver(&eqs);
		double g = eqs.ZeroFunction(solver.getTLast(), solver.getDataLast());
		return g;
	};
	clock_t timereq = clock();
	double v0crit = Find1DRoot(.00001, .9, zeroFunc, 1.e-8);
	timereq = clock() - timereq;
	std::cout << "v0=" << v0crit << ", time required: " << (float)timereq / CLOCKS_PER_SEC << std::endl;
	std::cout << zeroFunc(v0crit) << std::endl;
	IsothermalEquations eqs(a, v0crit);
	Integrator solver(&eqs);
	solver.writeToFile(&output);
	output.close();
}

