// toption_ene.cpp
// Copyright (c) 2011 KALX, LLC. All rights reserved. No warranty is made.
#include "../fmsdual/dual.h"
#include "option.h"
#include "exp_normal_exp.h"

using namespace option;

#ifdef max
#undef max
#endif

void
fms_test_option_ene(void)
{
	double eps = std::numeric_limits<double>::epsilon();

	double f = 100;
	double s = 0.2;
	double k = 100;
	double t = 0.25;

	double a = -2.0;
	double alpha = 2.0;
	double b = 2.0;
	double beta = 2.0;

	double p = black<>::value(f, s, -k, exp_normal_exp<>(t, a, alpha, b, beta));
	double c = black<>::value(f, s, +k, exp_normal_exp<>(t, a, alpha, b, beta));

	ensure (f - k == c - p);

	double srt = s*sqrt(t);
	double d1 = log(f/k)/srt + srt/2;
	double d2 = d1 - srt;
/*
	ensure (c == f*normal_cdf(d1) - k*normal_cdf(d2));
	ensure (p == k*normal_cdf(-d2) - f*normal_cdf(-d1));

	typedef dual::number<double,3> d3;
	d3 F(f,1), V;

	V = black<d3>::value<d3>(F, s, -k, exp_normal_exp<d3>(t, a, alpha, b, beta));
	double delta = -normal_cdf(-d1); // dv/df
	ensure (fabs(V.derivative(1) - delta) <= eps);					// TO DO change to printing difference
	ensure (fabs(V.derivative(1) - black<>::delta<>(f, s, -k, exp_normal_exp<>(t, a, alpha, b, beta))) <= eps);
	double gamma = normal_pdf(d1)/(f*srt); // d^2v/df^2
	ensure (fabs(V.derivative(2) - gamma) <= eps);

	V = black<d3>::value<d3>(F, s, k, exp_normal_exp<d3>(t, a, alpha, b, beta));
	delta = normal_cdf(d1);
	ensure (fabs(V.derivative(1) - delta) <= eps);
	ensure (fabs(V.derivative(1) - black<>::delta<>(f, s, k, exp_normal_exp<>(t, a, alpha, b, beta))) <= eps);
	ensure (fabs(V.derivative(2) - gamma) <= eps); // same gamma

	d3 S(s,1);
	V = black<d3>::value<d3>(f, S, -k, exp_normal_exp<d3>(t, a, alpha, b, beta));
	double vega = f*srt*normal_pdf(d1)/s; // dv/ds
	ensure (fabs(V.derivative(1) - vega) <= eps);

	d3 T(t,1);
	V = black<d3>::value<d3>(f, s, k, exp_normal_exp<d3,d3>(T, a, alpha, b, beta));
	double theta = f*srt*normal_pdf(d1)/(2*t); // dv/dt
	ensure (fabs(V.derivative(1) - theta) <= theta*eps);

	d3 K(-k,1);
	V = black<d3>::value<d3>(f, s, K, exp_normal_exp<d3>(t, a, alpha, b, beta));	
	double kappa = -black<>::binary<>(f, s, -k, exp_normal_exp<>(t, a, alpha, b, beta)); // dv/dk
	ensure (fabs(V.derivative(1) - kappa) <= eps);


	distribution_parameters_base<>* pd = new exp_normal_exp_parameters<>(t, a, alpha, b, beta);
	p = black<>::value(f, s, -k, pd->operator*());
*/
}

