#include<iostream>
#include<math.h>
#include<vector>
#include<fstream>
#include<ctime>
#include "Equations.h"
#include "Integrator.h"
using namespace std;

int main() {
	Equations eqs(1.01,10,200,1);
	//cout << "Compiled equations" << endl;
	//Equations eqs;
	Integrator In(eqs);

	//In.generateFunction(new double[4]{ 0,0,log(1E-8),0 }, 50);
	//In.writeToFile("./output.txt");
	//cout << In.findZeros(new double[4]{ 0,0,log(.000353531 - 2.98023E-8),0 }) << endl;
	//cout << In.findZeros(new double[4]{ 0,0,log(.000353531),0 }) << endl;
	
	//cout << "Compiled Integrator" << endl;
	//double* teststep = In.coupledRKStep(.1, new double[4]{ 1,1,1,1 });
	//for (int i = 0; i < 4; i++) { cout << teststep[i] << endl; };
	//In.generateFunction(new double[4]{ 0,0,log(.1),0 });
	//cout << "Generated functions" << endl;
	//In.findZeros(new double[4]{ 1,1,1,1 });
	//cout << "found zeros" << endl;
	//In.findV0();
	//cout << "found v0" << endl;
	//cout << "finding v0" << endl;
	clock_t timereq = clock();
	double v0 = In.findV0();
	timereq = clock() - timereq;
	cout << "v0 = " << v0 << ", time elapsed = " << (float)timereq/CLOCKS_PER_SEC << endl;
	//cout << In.findZeros(new double[4]{ 0,0,log(v0),0 }) << endl;
	//vector<double*> data = In.generateFunction(new double[4] { 0,0,log(v0),0 });
	//cout << "Final state: " << endl;
	//for (int i = 0; i < 4; i++) {
	//	cout << data.back()[i] << endl;
	//}
	//In.generateFunction(new double[4]{ 0,0,log(.0003755),0 },50);
	//In.generateFunction(new double[4]{ 0,0,log(v0),0 }, 50);
	//In.writeToFile("./output.txt");
	//cout << In.findZeros(new double[4]{ 0,0,log(.0003755),0 }) << endl;
	//cout << "scanning" << endl;
	In.scan("./gammahc.txt",.0003,.0004,10);
	return 0;
}

//trying to integrate without taking the absolute value results in an infinite loop that I haven't tracked down yet