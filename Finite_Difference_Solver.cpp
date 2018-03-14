#ifndef MVECTOR_H // the 'include guard'
#define MVECTOR_H // see C++ Primer Sec. 2.9.2

#include <thread>
#include <sstream>
#include <fstream>
#include <vector>
#include <iostream>
#include <cmath>

// Class that represents a mathematical vector
class MVector
{
public:
	// constructors
	MVector() {}
	explicit MVector(int n) : v(n) {}
	MVector(int n, double x) : v(n, x) {}

	// access element (lvalue)
	double &operator[](int index) { return v[index]; }

	// access element (rvalue)
	double operator[](int index) const { return v[index]; }

	int size() const { return v.size(); } // number of elements

private:
	std::vector<double> v;
};

//Operator overloads

//a function to compare MVector sizes
void MVector_Size_Comparison(MVector rhs, MVector lhs)
{
	if (lhs.size() != rhs.size())
	{
		std::cout << "error: MVectors are of different sizes :" << rhs.size() << " & " << " " << lhs.size() << std::endl;
		exit(0);
	}
}

// Operator overload for "scalar*vector" 
inline MVector operator*(const double& lhs, const MVector& rhs)
{
	MVector returnVal(rhs);
	for (int i = 0; i < rhs.size(); i++) returnVal[i] = lhs*rhs[i];
	return returnVal;
}

// Operator overload for "vector+vector" 
inline MVector operator+(const MVector& lhs, const MVector& rhs)
{
	MVector_Size_Comparison(rhs, lhs);
	MVector returnVal(rhs);
	for (int i = 0; i < rhs.size(); i++) returnVal[i] = lhs[i] + rhs[i];
	return returnVal;
}

// Operator overload for "vector-vector" 
inline MVector operator-(const MVector& lhs, const MVector& rhs)
{
	MVector_Size_Comparison(rhs, lhs);
	MVector returnVal(rhs);
	for (int i = 0; i < rhs.size(); i++) returnVal[i] = lhs[i] - rhs[i];
	return returnVal;
}

// Operator overload for "vector*vector" 
inline MVector operator*(const MVector& lhs, const MVector& rhs)
{
	MVector_Size_Comparison(rhs, lhs);
	MVector returnVal(rhs);
	for (int i = 0; i < rhs.size(); i++) returnVal[i] = lhs[i] * rhs[i];
	return returnVal;
}

// Operator overload for "vector/scalar" 
inline MVector operator/(const MVector& lhs, const double& rhs)
{
	MVector returnVal(lhs);
	for (int i = 0; i < lhs.size(); i++) returnVal[i] = lhs[i] / rhs;
	return returnVal;
}

//Operator overload for throwing MVector terms
std::ostream &operator<<(std::ostream &fss, const MVector &m)
{
	fss << m[0];
	for (int i = 1; i < m.size(); i++)
	{
		fss << "," << m[i];
	}
	return fss;
}

struct MFunction
{
	virtual MVector operator()(const double& x, const MVector& y) = 0;
};

double power(double x, int p) //a function to compute powers
{
	double temp = 1;
	if (p >= 0) { for (p; p > 0; p--) { temp = temp*x; } return temp; }
	else { p = -p; for (p; p > 0; p--) { temp = temp*x; } return 1 / temp; }
}

// An Euler scheme ODE solver, performs one itteration of the Euler method. specify a file name and it will create a .csv
int EulerSolve(int steps, double a, double b, MVector &y, MFunction &f, std::string fname = "")
{
	double h = ((b - a) / steps);
	if (h < 0) return 1; //[a,b] is not a valid domain

	if (fname != "") //the version that creates a .csv
	{
		std::ofstream file;
		file.open(fname + ".csv");
		if (!file) { std::cout << "error, file could not open" << std::endl; return 1; }
		file.precision(17);
		file << a << "," << y << std::endl;
		for (double x = a; x < b + h / 2; x += h)
		{
			y = y + h*f(x, y);
			file.precision(17);
			file << x << "," << y << std::endl;
		}
		file.close();

	}
	else //the regular version of the code
	{
		for (double x = a; x < b + h / 2; x += h)
		{
			y = y + h*f(x, y);
		}
	}
	return 0;
}

// A Midpoint method ODE solver, specify a file name and it will create a .csv
int MidpointSolve(int steps, double a, double b, MVector &y, MFunction &f, std::string fname = "")
{
	double h = ((b - a) / steps);
	if (h < 0) return 1; //[a,b] is not a valid domain

	if (fname != "") //the version that creates a .csv
	{
		std::ofstream file;
		file.open(fname + ".csv");
		if (!file) { std::cout << "error, file could not open" << std::endl; return 1; }
		file.precision(17);
		file << a << "," << y << std::endl;

		for (double x = a; x <= b + h / 2; x += h)
		{
			y = y + h*f(x + h / 2, y + h*f(x, y) / 2);
			file.precision(17);
			file << x << "," << y << std::endl;
		}
		file.close();
	}
	else //the regular version of the code
	{
		for (double x = b; x <= b + h / 2; x += h)
		{
			y = y + h*f(x + h / 2, y + h*f(x, y) / 2);
		}
	}
	return 0;
}
// A Runge-Kutta ODE solver, specify a file name and it will create a .csv
int RKSolve(int steps, double a, double b, MVector &y, MFunction &f, std::string fname = "")
{
	double h = ((b - a) / steps);
	if (h < 0) return 1; //[a,b] is not a valid domain
	MVector k1, k2, k3, k4;

	if (fname != "") //the version that creates a .csv
	{
		std::ofstream file;
		file.open(fname + ".csv");
		if (!file) { std::cout << "error, file could not open" << std::endl; return 1; }
		file.precision(17);
		file << a << "," << y << std::endl;

		for (double x = a; x <= b+h/2; x += h) //we add h/2 to our upper bound to compensate for floating point error
		{
			k1 = f(x, y);
			k2 = f(x + h / 2, y + h*k1 / 2);
			k3 = f(x + h / 2, y + h*k2 / 2);
			k4 = f(x + h, y + h*k3);
			y = y + h*(k1 + 2 * k2 + 2 * k3 + k4) / 6;
			file.precision(17);
			file << x << "," << y << std::endl;
		}
		file.close();
	}
	else //the regular version of the code
	{
		for (double x = a; x <= b + h / 2; x += h) //we add h/2 to our upper bound to compensate for floating point error
		{
			k1 = f(x, y);
			k2 = f(x + h / 2, y + h*k1 / 2);
			k3 = f(x + h / 2, y + h*k2 / 2);
			k4 = f(x + h, y + h*k3);
			y = y + h*(k1 + 2 * k2 + 2 * k3 + k4) / 6;
		}
	}
	return 0;
}



class F2 : public MFunction
{
public:
	virtual MVector operator()(const double& x, const MVector& y)
	{
		MVector temp(2);
		temp[0] = x;
		temp[1] = y[1];
		return temp;
	}
};

class ODE1 : public MFunction
{
public:
	virtual MVector operator()(const double& x, const MVector& y)
	{
		MVector temp(2);
		temp[0] = y[1];
		temp[1] = (32 + 2 * power(x, 3) - y[0] * y[1]) / 8.0;
		return temp;
	}
};

class Falkner_Skan : public MFunction
{
public:
	Falkner_Skan(double beta_ = 0.5) { beta = beta_; }
	virtual MVector operator()(const double& x, const MVector& y)
	{
		MVector temp(6);
		temp[0] = y[1];//df
		temp[1] = y[2]; //df'
		temp[2] = beta*(power(y[1], 2) - 1.0) - y[0] * y[2]; //df''
		temp[3] = y[4]; //dz1
		temp[4] = y[5]; //dz2
		temp[5] = beta * 2.0 * y[1] * y[4] - y[3] * y[2] - y[5] * y[0]; //dz3
		return temp;
	}
	void SetBeta(double beta_) { beta = beta_; }
private:
	double beta;
};

class BVP1 : public MFunction
{
public:
	MVector operator()(const double& x, const MVector& y)
	{
		MVector temp(4);
		temp[0] = y[1]; //Y terms
		temp[1] = (32 + 2 * power(x, 3) - y[0] * y[1]) / 8.0;
		temp[2] = y[3]; //Z terms
		temp[3] = -y[0] * y[3] /8.0 - y[1] * y[2]/8.0;
		return temp;
	}
};


double BVPSolve2(MFunction &f, int maxNewtonSteps, double a, double b, double y0, double guess, double tol = 1e-8)
{
	MVector y(4);
	double phi, phidash;
	for (int i = 0; i < maxNewtonSteps; i++)
	{
		// y[0] = y, y[1] = y', y[2] = Z_1, y[3] = Z_2,
		y[0] = y0; y[1] = guess; y[2] = 0; y[3] = 1; 
		RKSolve(100, a, b, y, f); // solve IVP
		phi = y[0] - 43/3.0; // calculate residual
		phidash = y[2]; // 'Jacobian' phidash = Z_1
		if (std::abs(phi) < tol) return guess; // exit if converged
		guess -= phi / phidash; // apply newton step
	}
	std::cout << std::endl << "Guess did not converge to within tolerance, ";
	RKSolve(100, a, b, y, f);
	if (y[0] - 43 / 3.0 <= phi)
	{
		std::cout << "guess is convering. choosing current guess : " << guess << std::endl;
		std::cout << "Current guess has absolute error: " << std::abs(y[1] - 1) << std::endl << std::endl;
		return guess;
	}
	else { std::cout << "guess is not converging. choosing guess = 1 :" << std::endl; return 1; }
}

int BVPSolve(MFunction &f, int maxNewtonSteps, double a, double b, double y0, double y1, double farBC, double &guess, double tol = 1e-8)
{
	MVector y(6);
	double phi, phidash;
	for (int i = 0; i < maxNewtonSteps; i++)
	{
		// y[0] = y, y[1] = y', y[2] = y'', y[3] = Z_1, y[4] = Z_2, y[5]= Z_3
		y[0] = y0; y[1] = y1; y[2] = guess; y[3] = 0.0; y[4] = 0.0; y[5] = 1.0;
		RKSolve(100, a, b, y, f); // solve IVP
		phi = y[1] - farBC; // calculate residual
		phidash = y[4]; // 'Jacobian' phidash = Z_1(nu - > infty)
		if (std::abs(phi) < tol) return 0; // exit if converged
		guess -= phi / phidash; // apply newton step
	}
	//what happens if the guess doenst converge within the max steps
	std::cout << std::endl << "Guess did not converge to within tolerance, ";
	RKSolve(100, a, b, y, f);
	if (y[1] - farBC <= phi) //cheking for convergence, making do if we have it
	{
		std::cout << "guess is convering. choosing current guess : " << guess << std::endl;
		std::cout << "Current guess has absolute error: " << std::abs(y[1] - farBC) << std::endl << std::endl;
		return 0;
	}
	else { std::cout << "guess is not converging. Current Guess is: " << guess << std::endl; return 1; } //return an error if we have divergence.
}
#endif

//Problem specific functions

void EvaluateFS(int farbound, double guess = 1.2, double beta = 0.5)
{
	Falkner_Skan f(beta); //ODE
	double a = 0.0; //domain
	std::cout << "beta = " << beta << std::endl;
	MVector Y(6);
	Y[0] = 0; //initial conditions of y
	Y[1] = 0; //initial conditions of y'
	Y[2] = guess;
	if (BVPSolve(f, 100, a, farbound, Y[0], Y[1], 1, Y[2])) //itterating the NR method
	{
		std::cout << "Could not compute Y[2]" <<std::endl; 
		exit(0); //canceling if we cannot compute a guess
	}
	std::cout << "Y[2] = " << Y[2] << std::endl;
	std::ostringstream convert; //dynamicly naming the file the results are saved to
	convert << "Flknrskn_b=" << beta;
	RKSolve(100, a, farbound, Y, f, convert.str()); //writing to file
	std::cout << "Y= " << Y[0] << " , " << Y[1] << " , " << Y[2] << std::endl << std::endl;
}


double NegBetaEval(double guess, double beta)
{
	Falkner_Skan f(beta); //ODE
	double a = 0.0, b = 10.0; //domain
	std::cout << "beta = " << beta << std::endl;
	MVector Y(6);
	Y[0] = 0; //initial conditions of y
	Y[1] = 0; //initial conditions of y'
	Y[2] = guess;
	if (BVPSolve(f, 100, a, b, Y[0], Y[1], 1, Y[2])) //itterating the NR method
	{
		return 0; //Mapping to zero if convergence fails
	}
	return Y[2]; 
}


//A function to test for valid boundarys of the falkner skan equations
int BoundaryConvergence(double guess, double tolerance = 1e-7,double beta = 0.5)
{
	Falkner_Skan f(beta); //ODE
	double a = 0.0, temp = 10.0; //temp is a record of the previous guess

	std::ofstream file;
	std::ostringstream convert; //dynamicly naming the file the results are saved to
	convert << "Convergence testing b = " << beta <<".csv";
	file.open(convert.str());
	if (!file) { std::cout << "error, file could not open" << std::endl; return 1; } //standard failure to open error

	for (double b = 1.0; b <= 15; b++) //varying the upper bound
	{
		std::cout << "Farbound is: " << b << std::endl;
		MVector Y(6);
		Y[0] = 0.0; //initial conditions of y
		Y[1] = 0.0; //initial conditions of y'
		Y[2] = guess;
		if (BVPSolve(f, 100, a, b, Y[0], Y[1], 1, guess)) //itterating the NR method
		{
			std::cout << "Could not compute Y[2]" << std::endl;
			exit(0);
		}
		std::cout.precision(17);
		std::cout << "Guess = " << Y[2] << std::endl << std::endl;
		file.precision(17);
		file << b << "," << Y[2];
		if (std::abs(Y[2] - temp) < tolerance) //when the change in the guess is less than the tolerance, we finish
		{
			std::cout << "Guess has converged, farbound is: "<< b << std::endl << std::endl; 
			break;
		}
		temp = Y[2];
	}
	file.close();
	return 0;
}

int main()
{
	//---------------------------------------------------------------------------------
	BVP1 f;
	double a = 1.0, b = 3.0;
	MVector Y(4);
	for (int j = 1; j <= 10; j++)
	{
		double i = power(j, 2);

		std::ostringstream EulerName, MidpointName, RKName; //dynamicly naming the file the results are saved to
		EulerName << "EulerData" << j;
		MidpointName << "MidpointData" << j;
		RKName << "RKData" << j;

		Y[0] = 17; //resetting the IC's
		Y[1] = 1;

		EulerSolve(i, a, b, Y, f, EulerName.str());

		Y[0] = 17;
		Y[1] = 1;

		MidpointSolve(i, a, b, Y, f, MidpointName.str());

		Y[0] = 17;
		Y[1] = 1;

		RKSolve(i, a, b, Y, f, RKName.str());

	}

	//---------------------------------------------------------------------------------
	//EvaluateFS(10,1.2);
	//---------------------------------------------------------------------------------
	/* Finding Boundary Convergence
	BoundaryConvergence(1.2, 1e-7, 1);
	BoundaryConvergence(1.2, 1e-8, 0.5);
	BoundaryConvergence(1.2, 1e-8, 0);
	*/
	//---------------------------------------------------------------------------------
	/*//attampted multithreading
	std::ofstream file;
	file.open("BetaPlot.csv");
	if (!file) { std::cout << "error, file could not open" << std::endl; return 1; }

	double a = 0.0, b = 10.0; //domain
	Falkner_Skan f; //function

	int BetaSteps = 100;

	int nThreads = 4; //multithreading
	std::vector<std::thread> threads(nThreads);

	for (int i = 1; i <= BetaSteps; i += nThreads)
	{
		MVector cache(nThreads, 1.2); //a queue of guess values to be written to file, 1.2 is the initial guess

		for (int j = 0; j < nThreads; j++)
		{
			if (j + i > BetaSteps) break;
			f.SetBeta((j+i)/BetaSteps);
			threads[j] = std::thread(BVPSolve, f, 100, a, b, 0, 0, 1, cache[j]);
		}

		for (int j = 0; j < nThreads; j++)
		{
			if (j + i > 10) break;
			threads[j].join();
			file.precision(17);
			file << cache[j] << std::endl;
		}
	}
	file.close();
	*/
	//---------------------------------------------------------------------------------
	/*
	std::ofstream file;
	file.open("BetaPlot.csv");
	if (!file) { std::cout << "error, file could not open" << std::endl; return 1; }

	double a = 0.0, b = 10.0; //domain
	Falkner_Skan f; //function

	double BetaSteps = 100;
	double guess = 1.2;

	MVector v(BetaSteps+1, guess);

	for (int i = 0; i <= BetaSteps+0.5; i++) //+0.5 because BetaSteps is a double
	{
		double k = i / BetaSteps; //BetaSteps is a double so that k is non-zero for i > 0
		std::cout << k << std::endl;
		f.SetBeta(k);
		if (BVPSolve(f, 100, a, b, 0, 0, 1, v[i]))
		{
			file.close();
			return 1;
		}
		file << v[i] << std::endl;
	}
	file.close();
	*/
	//---------------------------------------------------------------------------------
	/*
	//Computing Negative Beta
	std::ofstream file;
	file.open("NegativeBeta.csv");
	if (!file) { std::cout << "error, file could not open" << std::endl; return 1; }

	for (double i = -9; i <= 1; i++)
	{
		if (i > -8 && i < 1) continue; //these values do not give usfull results, so they are ommited
		double guess = i / 100.0;
		for (int j = 0; j <= 400; j++)
		{
			double beta = -j / 2000.0; //trial and error has shown that beta above this point fail, to verify this, increase the bound on j
			file.precision(17);
			file << NegBetaEval(guess, beta) << ",";
		}
		file << 0 << std::endl; //to complete the .csv, insignificant because NegBetaEval maps to zero for further terms anyway
	}
	file.close();
	*/
	//---------------------------------------------------------------------------------
	/*
	double b = 1;
	BVP1 f;
	MVector Y(4);
	Y[0] = 17; 
	Y[1] = BVPSolve2(f, 100, 1, 3, Y[0], 1);
	std::cout << "Guess: " << Y[1] << std::endl;
	RKSolve(100, 1, 3, Y, f);
	std::cout << "Y(3)= " << Y[0] << ", Y'(3)= " << Y[1] << std::endl << std::endl;
	*/
	
	/*
	//We are solving falkner skan with beta values ranging from 0 to 1
	Falkner_Skan f; //ODE
	double a = 0.0, b = 10.0; //domain
	for (int i = 0; i <= 50; i++)
	{
		double B = i / 50.0;
		f.SetBeta(B);
		std::cout << "beta = " << B << std::endl;
		MVector Y(6);
		Y[0] = 0; //initial conditions of y
		Y[1] = 0; //initial conditions of y'
		Y[2] = BVPSolve(f, 100, a, b, Y[0], Y[1], 1, 0.927684);
		std::cout << "Y[2] = " << Y[2] << std::endl;
		std::ostringstream convert; //dynamicly naming the file the results are saved to
		convert << "Flknr_scn B=" << B;
		RKSolve(100, a, b, Y, f,convert.str());
		std::cout << "Y= " << Y[0] << " , " << Y[1] << " , " << Y[2] << std::endl << std::endl;
	}
	std::cout << "FINISHED" << std::endl;

	/*while (true)
	{
		Falkner_Skan f; //ODE
		MVector Y(6);
		Y[0] = 0; //initial conditions of y
		Y[1] = 0; //initial conditions of y'
		Y[2] = BVPSolve(f, 100, a, b, Y[0], Y[1], 1.0); //initial conditions of y''
		std::cout << "f'' initial condition: " << Y[2] << std::endl;
		RKSolve(100, a, b, Y, f);
		std::cout << "Y= " << Y[0] << " , " << Y[1] << " , " << Y[2] << std::endl << std::endl;
		std::cout << std::endl;
		b += 0.5;
		if (b > 15.0) break;
	}*/
	return 0;
}
