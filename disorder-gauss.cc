// g++ -Wall -O3 disorder-realzn-GAUSS.cc -o testo

#include <iostream>
#include <fstream>
#include <ctime>
#include <boost/multi_array.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>

// gen is a variable name
// Its data-type is boost::random::mt19937
boost::random::mt19937 gen(std::time(0));
using namespace std;

typedef
boost::multi_array < double, 2 > array_2d;
// typedef keyword allows you to create an alias fo a data type

const unsigned int axis1 = 8 ;
const unsigned int axis2 = axis1;
// above assigns length along each dimension of the 2d configuration

//Function templates
double gaussian(double a, double b);

int main()
{
	array_2d J_x(boost::extents[axis1][axis2]);
	array_2d J_y(boost::extents[axis1][axis2]);
	//Assign random sign to each NN bond & store in an array
	ofstream fout("Jx.dat");	// Opens a file for output
	ofstream gout("Jy.dat");

	for (unsigned int i = 0; i < axis1; ++i)
	{
		for (unsigned int j = 0; j < axis2; ++j)
		{
			J_x[i][j] = gaussian(0, 1) ;
			J_y[i][j] = gaussian(0, 1) ;
			fout<< J_x[i][j] << endl;
			gout<< J_y[i][j] << endl;
		}
	}

	fout.close();
	gout.close();
	return 0;
}

//function to generate random integer
// from a gaussian distribution with mean-a, variance=b
double gaussian(double a, double b)
{
	boost::random::normal_distribution <> dist(a, b);
	return dist(gen);
}

// http://www.boost.org/doc/libs/1_56_0/doc/html/boost/random/normal_distribution.html#idp93575424-bb
