#include "wfg.h"

#include <cstring>
#include <string>
using std::string;
#include <vector>
using std::vector;
#include "Toolkit/ExampleProblems.h"
using namespace WFG::Toolkit;
using namespace WFG::Toolkit::Examples;
#include "Toolkit/TransFunctions.h"
#include "assert.h"


namespace
{
vector< double > problem_calc_fitness
(
		const vector< double >& z,
		const int k,
		const int M,
        const char* fn
)
{
	if ( strcmp("WFG1",fn) == 0 )
	{
		return Problems::WFG1( z, k, M );
	}	
	else if ( strcmp("WFG2",fn)== 0  )
	{
		return Problems::WFG2( z, k, M );
	}
	else if ( strcmp("WFG3",fn)== 0  )
	{
		return Problems::WFG3( z, k, M );
	}
	else if ( strcmp("WFG4",fn) == 0 )
	{
		return Problems::WFG4( z, k, M );
	}
	else if ( strcmp("WFG5",fn) == 0 )
	{
		return Problems::WFG5( z, k, M );
	}
	else if ( strcmp("WFG6",fn) == 0)
	{
		return Problems::WFG6( z, k, M );
	}
	else if ( strcmp("WFG7",fn) == 0 )
	{
		return Problems::WFG7( z, k, M );
	}
	else if ( strcmp("WFG8",fn) == 0 )
	{
		return Problems::WFG8( z, k, M );
	}
	else if ( strcmp("WFG9",fn) == 0 )
	{
		return Problems::WFG9( z, k, M );
	}
	else if ( strcmp("I1",fn) == 0 )
	{
		return Problems::I1( z, k, M );
	}
	else if ( strcmp("I2",fn) == 0 )
	{
		return Problems::I2( z, k, M );
	}
	else if ( strcmp("I3",fn) == 0 )
	{
		return Problems::I3( z, k, M );
	}
	else if ( strcmp("I4",fn) == 0 )
	{
		return Problems::I4( z, k, M );
	}
	else if ( strcmp("I5",fn) == 0)
	{
		return Problems::I5( z, k, M );
	}
	else
	{
		fprintf(stderr,"Unknown problem %s\n", fn );
		assert( false );
		return vector< double >();
	}
}
}

int wfg_eval(const double* x, int n, int k, int M, const char* problem, double* fit )
{
	vector< double > z;  // the decison vector
	vector< double > f;  // the fitness vector
	for( int i = 0; i < n; i++ )
	{
		z.push_back( x[i]*(2*(i+1)) );
	}	
	f = problem_calc_fitness( z, k, M, problem);
	for( int i = 0; i < (int)f.size(); i++ )
		fit[i] = f[i]; 
	
	return 0;	
}
