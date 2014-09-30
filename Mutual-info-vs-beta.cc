//  g++ -Wall -O3 `pkg-config --cflags --libs gsl tabdatrw-0.3 interp2dpp` Mutual-info-vs-beta.cc -o testo
// 2nd Renyi entropy for classical 2d Ising model in zero magnetic field
//Metropolis algorithm employed

#include <iostream>
#include <fstream>
#include <ctime>
#include <boost/multi_array.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <gsl/gsl_integration.h>
#include <boost/lexical_cast.hpp>
#include <interp2d.hpp> // For interpolation
#include <tabdatrw.hpp> // For tabdatr and tabdatw


// gen is a variable name
// Its data-type is boost::random::mt19937
boost::random::mt19937 gen;
using namespace std;
using boost::lexical_cast;
using boost::bad_lexical_cast;

unsigned int axis2 = 8;

//Function templates
double f (double , void * params);
double g (double , void * params);

int main(int argc, char const * argv[])
{
	double beta_min(0), beta_max(0), del_beta(0);

	if (argc != 4)
	{
		cout << "Expecting three inputs: beta_min, beta_max, del_beta."
		     << endl << "Got " << argc - 1 << endl;
		return 1;
	}

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

	beta_max += del_beta;
//	cout << "Enter minimum beta" << endl;
//	cin >> beta_min;
//	cout << "Enter maximum beta" << endl;
//	cin >> beta_max;
//	cout << "Enter increment of beta at each step" << endl;
//	cin >> del_beta;
	double mut_info(0); //mutual information I_2
	ofstream fout("I2-vs-beta-8.dat"); // Opens a file for output
	vvdouble vm = tabdatr("Em-8.dat", 2);//modified energy data
	interp_data idm(vm,1);
	vvdouble vn = tabdatr("E-8.dat", 2);//normal energy data
	interp_data idn(vn,1);
//	gsl_integration_workspace * w
//          = gsl_integration_workspace_alloc (1000);
	gsl_integration_cquad_workspace * w
	    = gsl_integration_cquad_workspace_alloc (100);
	size_t nevals = 1e3;

	for (double beta = beta_min; beta < beta_max + del_beta; beta += del_beta)
	{
		mut_info = 0 ;
		double term1(0), term2(0), term3(0), abs_error(0);
		gsl_function F;
		F.function = &f;
		F.params = &idn;
//Function: int gsl_integration_qags (const gsl_function * f,double a,double b, double epsabs, double epsrel,size_t limit,gsl_integration_workspace * workspace,double * result, double *abserr)
//gsl_integration_cquad (const gsl_function * f, double a, double b, double epsabs, double epsrel, gsl_integration_cquad_workspace * workspace, double * result, double * abserr, size_t * nevals)
		gsl_integration_cquad (&F, 0,     beta, 1e-6, 1e-4,
		                       w, &term2, &abs_error, &nevals);
		gsl_integration_cquad (&F, 0, 2.0*beta, 1e-6, 1e-4,
		                       w, &term3, &abs_error, &nevals);
		F.function = &g;
		F.params = &idm;
		gsl_integration_cquad (&F, 0, beta, 1e-6, 1e-4, w, &term1,
		                       &abs_error, &nevals);
		mut_info =2.0*term1 -2.0* term2 - term3;
		fout << beta << '\t' << mut_info /axis2 << endl;
	}

	fout.close();
	return 0;
}

//performing numerical integration using gsl
double f (double beta, void * params)
{
	interp_data p = *(interp_data *) params;
	return p.interp_akima(beta);
}

double g (double beta, void * params)
{
	interp_data p = *(interp_data *) params;
	return p.interp_akima(beta);
}



