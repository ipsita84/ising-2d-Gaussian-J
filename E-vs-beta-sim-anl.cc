// g++ -Wall -O3 E-vs-beta-sim-anl.cc -o normal
// Run with command line arguments, e.g. ./testo betamin betamax delbeta
// Considering 2d ising model in zero magnetic field with random J sign
//warming up system for first N_mc/10 loops
//averaging energy for the next N_mc updates
//incorporating SIMULATED ANNEALING

#include <iostream>
#include <fstream>
#include <ctime>
#include <boost/multi_array.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/lexical_cast.hpp>

// gen is a variable name
// Its data-type is boost::random::mt19937
boost::random::mt19937 gen;
using boost::lexical_cast;
using boost::bad_lexical_cast;
using namespace std;

typedef
boost::multi_array < int, 2 > array_2d;
// typedef keyword allows you to create an alias fo a data type

// Magnitude of J
double J = 1.0;
const unsigned int axis1 = 8;
const unsigned int axis2 = axis1;
// above assigns length along each dimension of the 2d configuration

//No.of Monte Carlo updates we want
unsigned int N_mc = 1e6;

//Function templates
int roll_coin(int a, int b);
double random_real(int a, int b);
double energy_tot(array_2d sitespin, array_2d J_x, array_2d J_y);
double nn_energy(array_2d sitespin,  array_2d J_x, array_2d J_y, unsigned int row, unsigned int col);

int main(int argc, char const * argv[])
{
	if (argc != 4)
	{
		cout << "Expecting three inputs: beta_min, beta_max, del_beta."
		     << endl << "Got " << argc - 1 << endl;
		return 1;
	}

	array_2d J_x(boost::extents[axis1][axis2]);
	array_2d J_y(boost::extents[axis1][axis2]);
	//Read the random signed bonds for a particular stored realization
	ifstream gxin("Jx.dat");
	ifstream gyin("Jy.dat");

	for (unsigned int i = 0; i < axis1; ++i)
	{
		for (unsigned int j = 0; j < axis2; ++j)
		{
			gxin>>J_x[i][j];
			gyin>>J_y[i][j];
		}
	}

	gxin.close();
	gyin.close();

	double beta_min(0), beta_max(0), del_beta(0);

	try
	{
		beta_min = lexical_cast<double>(argv[1]);
		beta_max = lexical_cast<double>(argv[2]);
		del_beta = lexical_cast<double>(argv[3]);
	}
	catch (const bad_lexical_cast & x)
	{
		cout << "Cannot convert input to double" << endl;
		return 2;
	}

//	cout << "Enter minimum beta" << endl;
//	cin >> beta_min;
//	cout << "Enter maximum beta" << endl;
//	cin >> beta_max;
//	cout << "Enter increment of beta" << endl;
//	cin >> del_beta;
	ofstream fout("E.dat");	// Opens a file for output
	ofstream gout("EA.dat");

//      Create a 2d array that is axis1 * axis2
	array_2d sitespin(boost::extents[axis1][axis2]);
//      stores the spin configuration of the system

	array_2d sitespinsum(boost::extents[axis1][axis2]);


//      initial state chosen by random no. generator above
	for (unsigned int i = 0; i < axis1; ++i)
		for (unsigned int j = 0; j < axis2; ++j)
			sitespin[i][j] = 2 * roll_coin(0, 1) - 1;

	double energy = energy_tot(sitespin, J_x, J_y);
	
	fout << 0 << '\t' << 0 << endl;
	gout << 0 << '\t' << 0 << endl;

	for (double beta =beta_min + del_beta; beta < beta_max + del_beta; beta += del_beta)
	{
		double en_sum(0), EA(0);
		unsigned int sys_size = axis1 * axis2;

		for (unsigned int k = 0; k < axis1; ++k)
			for (unsigned int l = 0; l < axis2; ++l)
			sitespinsum[k][l] = 0;

		for (unsigned int i = 1; i <=1e5+N_mc; ++i)
		{
			for (unsigned int j = 1; j <= sys_size; ++j)
			{
				
//				Now choose a random spin site with site no.=label
				unsigned int label, row, col ;
				label = roll_coin(1, sys_size);

				if (label % axis2 == 0)
				{
					row = (label / axis2) - 1;
					col = axis2 -1 ;
				}
				else
				{
					col = label % axis2 - 1;
					row = (label-col-1)/axis2;
				}

				double energy_diff =-2 * nn_energy(sitespin,J_x, J_y, row, col);
				//Generate a random no. r such that 0 < r < 1
				double r = random_real(0, 1);
				double acc_ratio = exp(-1.0 * energy_diff* beta);

				//Spin flipped if r <= acceptance ratio
				if (r <= acc_ratio)
				{
					sitespin[row][col] *= -1;
					energy += energy_diff;
				}
			}

			if (i > 1e5) 
			{	en_sum += energy;
				for (unsigned int k = 0; k < axis1; ++k)
					for (unsigned int l = 0; l < axis2; ++l)
					sitespinsum[k][l] += sitespin[k][l] ;
			}
		}

		fout << beta << '\t' << en_sum / N_mc << endl;

		for (unsigned int i = 0; i < axis1; ++i)
			for (unsigned int j = 0; j < axis2; ++j)
			EA += sitespinsum[i][j]*sitespinsum[i][j] / (N_mc*N_mc);
//prints Edwards-Anderson order parameter for each beta value
		gout << beta << '\t' << EA / sys_size << endl;
	}

	fout.close();
	gout.close();
	return 0;
}

//function to generate random integer
// between 2 integers a & b, including a & b
int roll_coin(int a, int b)
{
	boost::random::uniform_int_distribution <> dist(a, b);
	return dist(gen);
}

//function to generate random real no.
// between 2 integers a & b, including a & excluding b

double random_real(int a, int b)
{
	boost::random::uniform_real_distribution <> dist(a, b);
	// uniform_real_distribution: continuous uniform distribution
	//on some range [min, max) of real number
	return dist(gen);
}

//function to calculate total energy
//for a given spin configuration
//with periodic boundary conditions

double energy_tot(array_2d sitespin, array_2d J_x, array_2d J_y)
{
	double energy = 0;

	for (unsigned int i = 0; i < axis1 - 1; ++i)
	{
		for (unsigned int j = 0; j < axis2 - 1; ++j)
		{
			energy += J_x[i][j]*sitespin[i][j]*sitespin[i+1][j];
			energy += J_y[i][j]*sitespin[i][j]*sitespin[i][j+1];
		}
	}

	//periodic boundary conditions
	for (unsigned int j = 0; j < axis2-1; ++j) // for i=axis1-1
	 {energy += J_x[axis1-1][j]*sitespin[axis1-1][j] * sitespin[0][j];
	  energy += J_y[axis1-1][j]*sitespin[axis1-1][j]*sitespin[axis1-1][j+1];
	 }

	for (unsigned int i = 0; i < axis1-1; ++i) // for j=axis2-1
	 {energy += J_y[i][axis2-1]*sitespin[i][axis2-1] * sitespin[i][0];
	  energy += J_x[i][axis2-1]*sitespin[i][axis2-1]*sitespin[i+1][axis2-1];
	 }

	energy += J_x[axis1-1][axis2-1]*sitespin[axis1-1][axis2-1]*sitespin[0][axis2-1];
	energy += J_y[axis1-1][axis2-1]*sitespin[axis1-1][axis2-1]*sitespin[axis1-1][0];
	
	

	return energy;
}

//Calculating interaction energy change for spin
//at random site->(row,col) with its nearest neighbours
double nn_energy(array_2d sitespin,  array_2d J_x, array_2d J_y, unsigned int row, unsigned int col)
{
	double nn_en = 0;

	if (row > 0 && row < axis1 - 1)
	{
		nn_en +=J*J_x[row-1][col] * sitespin[row][col] * sitespin[row-1][col];
		nn_en +=J*J_x[row][col]*sitespin[row][col] * sitespin[row+1][col];
	}

	if (col > 0 && col < axis2 - 1)
	{
		nn_en +=J*J_y[row][col-1]*sitespin[row][col] * sitespin[row][col-1];
		nn_en +=J*J_y[row][col]*sitespin[row][col] * sitespin[row][col+1];
	}

	if (row == 0)
	{
		nn_en +=J*J_x[axis1-1][col]*sitespin[0][col] * sitespin[axis1-1][col];
		nn_en += J *J_x[0][col]* sitespin[0][col] * sitespin[1][col];
	}

	if (row == axis1-1)
	{
		nn_en +=J*J_x[axis1-2][col] * sitespin[axis1 - 1][col] * sitespin[axis1-2][col];
		nn_en +=J*J_x[axis1-1][col] * sitespin[axis1-1][col] * sitespin[0][col];
	}

	if (col == 0)
	{
		nn_en +=J*J_y[row][axis2-1] * sitespin[row][0] * sitespin[row][axis2-1];
		nn_en += J*J_y[row][0] * sitespin[row][0] * sitespin[row][1];
	}

	if (col == axis2-1)
	{
		nn_en +=J*J_y[row][axis2-2] * sitespin[row][axis2-1] * sitespin[row][axis2-2];
		nn_en +=J*J_y[row][axis2-1] * sitespin[row][axis2-1] * sitespin[row][0];
	}

	return nn_en;
}
