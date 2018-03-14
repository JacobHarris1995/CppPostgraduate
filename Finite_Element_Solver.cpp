#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

#include <thread>
//a class for defining the initial profiles
struct function 
{
	virtual double operator()(const double& x) = 0;
};

class sinfunct : public function
{
public:
	virtual double operator()(const double& x){ return 1.5+std::sin(x); }
};

class squaresinfunct : public function
{
public:
	virtual double operator()(const double& x)
	{
		if (x <= 1)
			return 1;
		else
			return 0;
	}
};

//Base class for Advection and burgers elements
class UElement
{
public:
	// Pointer to the left neighbour
	UElement *Left_neighbour_pt;
	// Pointer to the right neighbour
	UElement *Right_neighbour_pt;
	// Storage for the coordinates
	std::vector<double> X;
	// Storage for the unknowns
	std::vector<double> U;
	// Constructor: initialising the vectors to hold two entries.
	UElement() 
	{	
		std::vector<double> X_(2), U_(2);
		X = X_;
		U = U_;
	}
	// Returning the value of the coordinate at local coordinate s using
	// equation (1.2)
	virtual double interpolated_x(double s) { return ((1.0 - s)*X[0] / 2.0 + (1.0 + s)*X[1] / 2.0); }
	// Returning the value of the unknown at local coordinate s using
	// equation (1.4)
	virtual double interpolated_u(double s) { return ((1.0 - s)*U[0] / 2.0 + (1.0 + s)*U[1] / 2.0); }
	//Calculate the	flux
	virtual double flux(const double u) = 0; //these are arbitrary and are defined in the subclasses
	virtual double dflux(const double u) = 0;
	//intergrate the flux
	virtual double integrate_flux()
	{
		return flux(interpolated_u(-1.0 / std::sqrt(3))) + flux(interpolated_u(1.0 / std::sqrt(3)));
	}
	//calculate h
	virtual double h(double a, double b)
	{
		//approximate the max of dflux
		double step = std::abs(b - a) / 100.0, max = 0.0, temp = 0.0;
		for (int i = 0; i < 100; i++)
		{
			temp = std::abs(dflux(interpolated_u(a + i*step)));
			if (temp > max) { max = temp; }
		}
		
		return 0.5*(flux(a) + flux(b) - max*(b-a));
	}
}; //End of the class definition

class AdvectionElement : public UElement
{
	virtual double flux(const double u) { return u; }
	virtual double dflux(const double u) { return 1; }
};


class BurgersElement : public UElement
{
	virtual double flux(const double u) { return u*u*0.5; }
	virtual double dflux(const double u) { return u; }
};

//the base class for Advection vectors and Burgers vectors
class AdvectionVector
{
public:
	AdvectionVector(function &f, int N_, double a, double b)
	{
		N = N_;
		std::vector<AdvectionElement> V_(N);
		V = V_;
		double interval = (b - a) / N;
		//seting up pointers
		V[0].Left_neighbour_pt = &V[N - 1];
		V[0].Right_neighbour_pt = &V[1];
		//interval
		V[0].X[0] = a;
		V[0].X[1] = V[0].X[0] + interval;
		//initial guesses
		V[0].U[0] = f(V[0].X[0]);
		V[0].U[1] = f(V[0].X[1]);
		//the body of elements
		for (int el = 1; el < N - 1; el++)
		{
			//seting up pointers
			V[el].Left_neighbour_pt = &V[el - 1];
			V[el].Right_neighbour_pt = &V[el + 1];
			//interval
			V[el].X[0] = V[el - 1].X[1];
			V[el].X[1] = V[el].X[0] + interval;
			//initial guesses
			V[el].U[0] = f(V[el].X[0]);
			V[el].U[1] = f(V[el].X[1]);
		}
		//the final element
		//seting up pointers
		V[N - 1].Left_neighbour_pt = &V[N - 2];
		V[N - 1].Right_neighbour_pt = &V[0];
		//interval
		V[N - 1].X[0] = V[N - 2].X[1];
		V[N - 1].X[1] = b;
		//initial guesses
		V[N - 1].U[0] = f(V[N - 1].X[0]);
		V[N - 1].U[1] = f(V[N - 1].X[1]);
	}
	//operator overload for []
	AdvectionElement operator[](int i) const { return V[i]; }
	//return the size of the advection vector
	virtual int size() const { return N; }
	//itterate the timestep
	void timestep(double dt)
	{
		std::vector<double> u0(N), u1(N);
		double U1, U0, h0, h1;
		for (int i = 0; i < N; i++)
		{		

			U1 = (*V[i].Left_neighbour_pt).U[1];
			U0 = (*V[i].Right_neighbour_pt).U[0];

			h0 = V[i].h(U1, V[i].U[0]);
			h1 = V[i].h(V[i].U[1], U0);
			//iterating the algorithm
			u0[i] = 2.0 * (-1.5*V[i].integrate_flux() + 2.0*h0 + h1) / (V[i].X[1] - V[i].X[0]);
			u1[i] = 2.0 * (1.5*V[i].integrate_flux() - h0 - 2.0*h1) / (V[i].X[1] - V[i].X[0]);
		}
		//updating the stored U values
		for (int i = 0; i < N; i++)
		{
			V[i].U[0] += dt*u0[i];
			V[i].U[1] += dt*u1[i];
		}
	}
private:
	int N;
	std::vector<AdvectionElement> V;
};

class BurgersVector
{
public:
	BurgersVector(function &f, int N_, double a, double b)
	{
		N = N_;
		std::vector<BurgersElement> V_(N);
		V = V_;
		double interval = (b - a) / N;
		//seting up pointers
		V[0].Left_neighbour_pt = &V[N - 1];
		V[0].Right_neighbour_pt = &V[1];
		//interval
		V[0].X[0] = a;
		V[0].X[1] = V[0].X[0] + interval;
		//initial guesses
		V[0].U[0] = f(V[0].X[0]);
		V[0].U[1] = f(V[0].X[1]);
		//the body of elements
		for (int el = 1; el < N - 1; el++)
		{
			//seting up pointers
			V[el].Left_neighbour_pt = &V[el - 1];
			V[el].Right_neighbour_pt = &V[el + 1];
			//interval
			V[el].X[0] = V[el - 1].X[1];
			V[el].X[1] = V[el].X[0] + interval;
			//initial guesses
			V[el].U[0] = f(V[el].X[0]);
			V[el].U[1] = f(V[el].X[1]);
		}
		//the final element
		//seting up pointers
		V[N - 1].Left_neighbour_pt = &V[N - 2];
		V[N - 1].Right_neighbour_pt = &V[0];
		//interval
		V[N - 1].X[0] = V[N - 2].X[1];
		V[N - 1].X[1] = b;
		//initial guesses
		V[N - 1].U[0] = f(V[N - 1].X[0]);
		V[N - 1].U[1] = f(V[N - 1].X[1]);
	}
	//operator overload for []
	BurgersElement operator[](int i) const { return V[i]; }
	virtual int size() const { return N; }
	//itterate the timestep
	void timestep(double dt)
	{
		std::vector<double> u0(N), u1(N);
		double U1, U0, h0, h1;
		for (int i = 0; i < N; i++)
		{

			U1 = (*V[i].Left_neighbour_pt).U[1];
			U0 = (*V[i].Right_neighbour_pt).U[0];

			h0 = V[i].h(U1, V[i].U[0]);
			h1 = V[i].h(V[i].U[1], U0);
			//iterating the algorithm
			u0[i] = 2.0 * (-1.5*V[i].integrate_flux() + 2.0*h0 + h1) / (V[i].X[1] - V[i].X[0]);
			u1[i] = 2.0 * (1.5*V[i].integrate_flux() - h0 - 2.0*h1) / (V[i].X[1] - V[i].X[0]);
		}
		//updating the stored U values
		for (int i = 0; i < N; i++)
		{
			V[i].U[0] += dt*u0[i];
			V[i].U[1] += dt*u1[i];
		}
	}
private:
	int N;
	std::vector<BurgersElement> V;
};


//operator overload for displaying AdvectionVectors
std::ostream &operator<<(std::ostream &oss, const AdvectionVector &UV)
{
	for (int i = 0; i < UV.size(); i++)
	{
		oss << UV[i].interpolated_x(0) << "," << (UV[i]).interpolated_u(0) << std::endl;
	}
	return oss;
}

//operator overload for displaying BurgersVectors
std::ostream &operator<<(std::ostream &oss, const BurgersVector &UV)
{
	for (int i = 0; i < UV.size(); i++)
	{
		oss << UV[i].interpolated_x(0) << "," << (UV[i]).interpolated_u(0) << std::endl;
	}
	return oss;
}

//we create a record of all x at all time steps from t=0 to t = max
//this lets us plot an animation of the wave moving.
int main() 
{
	squaresinfunct g; //initial profile
	sinfunct f;

	int interval = 100; //number of elements
	double a = 0.0, b = 2.0 * 3.14159265358979323846; //domain
	double MTime = 10.0; //Max time
	double dt = 0.4/(interval); //timestep, based on cfl condition
	double Max = MTime/dt; //how many itterations

	BurgersVector Data(f, interval, a, b);
	//AdvectionVector Data(f, interval, a, b);
	std::ofstream file;
	file.open("File_Name.csv");
	if (!file) { std::cout << "Could not open file for editing" << std::endl; return 1; }
	
	file << Data << std::endl; //saving the initial profile
	for (int j = 1; j < Max; j++)
	{
		Data.timestep(dt);
		file << Data << std::endl;
		std::cout << "Progress at: " << 100 * (j + 1) / Max << "%" << std::endl;
	}
	std::cout << "Finished, generated up to time: " << floor(Max)*dt << std::endl;
	file.close();
	return 0;
}
