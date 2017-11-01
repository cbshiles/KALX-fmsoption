// normalexp.h - normal in the middle, exponential at both ends
#pragma once
#include "../fmsdual/dual.h"

#pragma warning(push)
#pragma warning(disable: 4100) // unreferrence formal parameter

// Area under g(x) = f(a) exp(alpha (x - a), x < a; f(x), a < x < b; f(b) exp(beta (x - b)), x > b,
// where f is the standard normal density function
inline
double exp_normal_exp(double a, double alpha, double b, double beta)
{
	ensure (a < 0);
	ensure (alpha > 0);
	ensure (b > 0);
	ensure (beta < 0);

	return normal_pdf(a)/alpha + normal_cdf(b) - normal_cdf(a) - normal_pdf(b)/beta;
}
// h(x) = g(x)/exp_normal_exp
// cumulative distribution function of h
template<class T>
inline 
T exp_normal_exp_cdf(T x, double a, double alpha, double b, double beta)
{
	T ans = -1.0;
	if (x < a)
		ans = -normal_pdf(a)*exp( -a*(x-a) )/( a*exp_normal_exp(a,alpha,b,beta) );
	else if (x < b)
		ans = (-normal_pdf(a)/a + normal_cdf(x) - normal_cdf(a) )/exp_normal_exp(a,alpha,b,beta);
	else
		ans = (-normal_pdf(a)/a + normal_cdf(b) - normal_cdf(a) - normal_pdf(b)*exp( -b*(x-b) )/b + normal_pdf(b)/b)/exp_normal_exp(a,alpha,b,beta);
	return ans;
}

// integral to x of exp(sigma x)h(x) dx
// assumptions: assumption of a<0<b and a = -alpha and b = beta
template<class T>
inline
T normalexp_exp(T x, double sigma, double a, double alpha, double b, double beta)
{
	T ans = -1.0;
	if (x < a)
		ans = normal_pdf(a)*exp( (sigma - a)*x + a*a )/(exp_normal_exp(a,alpha,b,beta) * (sigma - a));
	else if (x < b) {
		ans = normal_pdf(a)*exp(a*sigma)/(sigma-a) + exp(.5*sigma*sigma)*(normal_cdf(x-sigma) - normal_cdf(a-sigma));
		ans = ans/exp_normal_exp(a,alpha,b,beta);
	}
	else {
		ans = normal_pdf(a)/(sigma-a);
		ans = ans + exp(.5*sigma*sigma)*(normal_cdf(b-sigma) - normal_cdf(a-sigma));
		ans = ans + normal_pdf(b)*(exp((sigma-b)*x + b*b)  - exp(sigma*b))/(sigma-b);
		ans = ans/exp_normal_exp(a,alpha,b,beta);
	}
	return ans;
}

// If F = f exp(ct + sigma sqrt(t) X), E[max{k - S}, F] = ? normalexp_cdf - ? normalexp_exp
// Choose c so thate E[S] = f.

#pragma warning(pop)