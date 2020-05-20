#pragma once

class EoS {
public:
	virtual int size() { return 0; };
	//return differentials of state quantities
	virtual std::vector<double> getDiffs(double t, std::vector<double> vals) = 0; 
	virtual std::vector<double> getInitialState() = 0;
	virtual double stopCondition(std::vector<double> state) = 0;
};