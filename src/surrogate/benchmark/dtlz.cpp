#include "dtlz.h"

namespace dtlz
{
#define PI 3.1415926535897932384626433832795

void DTLZ1(const double *x, double* f, const unsigned no_vars, const unsigned no_objs)
{
    int i = 0;
    int j = 0;
    int n = no_vars;
    int k = n - no_objs + 1;

    double g = 0;
    for (i = n - k + 1; i <= n; i++)
    {
	g += pow(x[i-1]-0.5,2) - cos(20 * PI * (x[i-1]-0.5));
    }
    g = 100 * (k + g);

    for (i = 1; i <= no_objs; i++)
    {
	double fi = 0.5 * (1 + g);
	for (j = no_objs - i; j >= 1; j--)
	{
	    fi *= x[j-1];
	}
	if (i > 1)
	{
	    fi *= 1 - x[(no_objs - i + 1) - 1];
	}

	f[i-1] = fi;
    }
}

void DTLZ2(const double *x, double* f, const unsigned no_vars, const unsigned no_objs)
{
    int i = 0;
    int j = 0;
    int n = no_vars;
    int k = n - no_objs + 1;

    double g = 0;
    for (i = n - k + 1; i <= n; i++)
    {
	g += pow(x[i-1]-0.5,2);
    }

    for (i = 1; i <= no_objs; i++)
    {
	double fi = (1 + g);
	for (j = no_objs - i; j >= 1; j--)
	{
	    fi *= cos(x[j-1] * PI / 2);
	}
	if (i > 1)
	{
	    fi *= sin(x[(no_objs - i + 1) - 1] * PI / 2);
	}

	f[i-1] = fi;
    }
}



void DTLZ3(const double *x, double* f, const unsigned no_vars, const unsigned no_objs)
{
    int i = 0;
    int j = 0;
    int n = no_vars;
    int k = n - no_objs + 1;

    double g = 0;
    for (i = n - k + 1; i <= n; i++)
    {
	g += pow(x[i-1]-0.5,2) - cos(20 * PI * (x[i-1]-0.5));
    }
    g = 100 * (k + g);

    for (i = 1; i <= no_objs; i++)
    {
	double fi = (1 + g);
	for (j = no_objs - i; j >= 1; j--)
	{
	    fi *= cos(x[j-1] * PI / 2);
	}
	if (i > 1)
	{
	    fi *= sin(x[(no_objs - i + 1) - 1] * PI / 2);
	}

	f[i-1] = fi;
    }
}

void DTLZ4(const double *x, double* f, const unsigned no_vars, const unsigned no_objs)
{
    int i = 0;
    int j = 0;
    double alpha = 100;
    int n = no_vars;
    int k = n - no_objs + 1;

    double g = 0;
    for (i = n - k + 1; i <= n; i++)
    {
	g += pow(x[i-1]-0.5,2);
    }

    for (i = 1; i <= no_objs; i++)
    {
	double fi = (1 + g);
	for (j = no_objs - i; j >= 1; j--)
	{
	    fi *= cos(pow(x[j-1],alpha) * PI / 2);
	}
	if (i > 1)
	{
	    fi *= sin(pow(x[(no_objs - i + 1) - 1],alpha) * PI / 2);
	}

	f[i-1] = fi;
    }
}

void DTLZ5(const double *x, double* f, const unsigned no_vars, const unsigned no_objs)
{
    int i = 0;
    int j = 0;
    int n = no_vars;
    int k = n - no_objs + 1;
    double *theta = new double[no_objs];
    double t = 0;
    double g = 0;

    for (i = n - k + 1; i <= n; i++)
    {
	g += pow(x[i-1] - 0.5, 2);
    }

    t = PI / (4 * (1 + g));
    theta[0] = x[0] * PI / 2;
    for (i = 2; i <= no_objs - 1; i++)
    {
	theta[i-1] = t * (1 + 2 * g * x[i-1]);
    }

    for (i = 1; i <= no_objs; i++)
    {
	double fi = (1 + g);
	for (j = no_objs - i; j >= 1; j--)
	{
	    fi *= cos(theta[j-1]);
	}
	if (i > 1)
	{
	    fi *= sin(theta[(no_objs - i + 1) - 1]);
	}

	f[i-1] = fi;
    }

    delete theta;
}

void DTLZ6(const double *x, double* f, const unsigned no_vars, const unsigned no_objs)
{
    int i = 0;
    int j = 0;
    int n = no_vars;
    int k = n - no_objs + 1;
    double *theta = new double[no_objs];
    double t = 0;
    double g = 0;

    for (i = n - k + 1; i <= n; i++)
    {
	g += pow(x[i-1], 0.1);
    }

    t = PI / (4 * (1 + g));
    theta[0] = x[0] * PI / 2;
    for (i = 2; i <= no_objs - 1; i++)
    {
	theta[i-1] = t * (1 + 2 * g * x[i-1]);
    }

    for (i = 1; i <= no_objs; i++)
    {
	double fi = (1 + g);
	for (j = no_objs - i; j >= 1; j--)
	{
	    fi *= cos(theta[j-1]);
	}
	if (i > 1)
	{
	    fi *= sin(theta[(no_objs - i + 1) - 1]);
	}

	f[i-1] = fi;
    }

    delete theta;
}

void DTLZ7(const double *x, double* f, const unsigned no_vars, const unsigned no_objs)
{
    int i = 0;
    int j = 0;
    int n = no_vars;
    int k = n - no_objs + 1;
    double g = 0;
    double h = 0;

    for (i = n - k + 1; i <= n; i++)
    {
	g += x[i-1];
    }
    g = 1 + 9 * g / k;

    for (i = 1; i <= no_objs - 1; i++)
    {
	f[i-1] = x[i-1];
    }

    for (j = 1; j <= no_objs - 1; j++)
    {
	h += x[j-1] / (1 + g) * (1 + sin(3 * PI * x[j-1]));
    }
    h = no_objs - h;
    f[no_objs - 1] = (1 + g) * h;
}

} /* dtlz */

