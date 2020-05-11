#include<iostream>
#include<math.h>
#include<vector>
#include<fstream>
#include<ctime>
#include "Equations.h"
#include "Integrator.h"
using namespace std;

int main() {
	//cout << "Hello world" << endl;
	Equations eqs(1.01,10,200);
	//cout << "Compiled equations" << endl;
	//Equations eqs;
	//cout << eqs.dw(new double[4] { 1,1,1,1 }) << endl;
	Integrator In(eqs);
	//cout << "Compiled Integrator" << endl;
	//double* teststep = In.coupledRKStep(.1, new double[4]{ 1,1,1,1 });
	//for (int i = 0; i < 4; i++) { cout << teststep[i] << endl; };
	//In.generateFunction(new double[4]{ 0,0,log(.1),0 });
	//cout << "Generated functions" << endl;
	//In.findZeros(new double[4]{ 1,1,1,1 });
	//cout << "found zeros" << endl;
	//In.findV0();
	//cout << "found v0" << endl;
	double x = exp(100000);
	double y = exp(x);
	if (x > 1000) { cout << x+7 << endl; };
	clock_t timereq = clock();
	double v0 = In.findV0();
	timereq = clock() - timereq;
	cout << "v0 = " << v0 << ", time elapsed = " << (float)timereq/CLOCKS_PER_SEC << endl;
	In.generateFunction(new double[4]{ 0,0,.001,0 }, 50);
	In.writeToFile("./gammaheating2.txt");
	//cout << In.coupledRKStep(.001, new double[4]{ .00875, });
	//vector<double*> data = In.generateFunction(new double[4] { 0,0,log(v0),0 });
	//cout << "Final state: " << endl;
	//for (int i = 0; i < 4; i++) {
	//	cout << data.back()[i] << endl;
	//}
	//In.writeToFile("./output.txt");
	In.scan("./gammaheating.txt",1E-10,.5,10);
	cout << exp(numeric_limits<double>::max());
	return 0;
}

//trying to integrate without taking the absolute value results in an infinite loop that I haven't tracked down yet