// toption.cpp
// Copyright (c) 2011 KALX, LLC. All rights reserved. No warranty is made.
#include <ctime>
#include <functional>
#include "../fmsdual/dual.h"
#include "option.h"
#include "standard_brownian.h"

using namespace option;

#ifdef max
#undef max
#endif

void
fms_test_brownian(void)
{
	double eps = std::numeric_limits<double>::epsilon();

	double f = 100;
	double s = 0.2;
	double k = 100;
	double t = 0.25;
	double p = black<>::value<>(f, s, -k, standard_brownian<>(t));
	double c = black<>::value<>(f, s, +k, standard_brownian<>(t));

	ensure (f - k == c - p);

	double srt = s*sqrt(t);
	double d1 = log(f/k)/srt + srt/2;
	double d2 = d1 - srt;

	ensure (c == f*normal_cdf(d1) - k*normal_cdf(d2));
	ensure (p == k*normal_cdf(-d2) - f*normal_cdf(-d1));

	typedef dual::number<double,3> d3;
	d3 F(f,1), V;

	V = black<d3>::value<>(F, s, -k, standard_brownian<d3,double,double,d3>(t));
	double delta = -normal_cdf(-d1); // dv/df
	ensure (fabs(V.derivative(1) - delta) <= eps);
	ensure (fabs(V.derivative(1) - black<>::delta<>(f, s, -k, standard_brownian<>(t))) <= eps);
	double gamma = normal_pdf(d1)/(f*srt); // d^2v/df^2
	ensure (fabs(V.derivative(2) - gamma) <= eps);


	V = black<d3>::value<d3>(F, s, k, standard_brownian<d3,double,double,d3>(t));
	delta = normal_cdf(d1);
	ensure (fabs(V.derivative(1) - delta) <= eps);
	ensure (fabs(V.derivative(1) - black<>::delta<>(f, s, k, standard_brownian<>(t))) <= eps);
	ensure (fabs(V.derivative(2) - gamma) <= eps); // same gamma

	d3 S(s,1);
	V = black<d3>::value<d3>(f, S, -k, standard_brownian<d3,double,d3,d3>(t));
	double vega = f*srt*normal_pdf(d1)/s; // dv/ds
	ensure (fabs(V.derivative(1) - vega) <= f*eps);

	d3 T(t,1);
	V = black<d3>::value<d3>(f, s, k, standard_brownian<d3,d3,double,d3>(T));
	double theta = f*srt*normal_pdf(d1)/(2*t); // dv/dt
	ensure (fabs(V.derivative(1) - theta) <= theta*eps);

	d3 K(-k,1);
	V = black<d3>::value<d3>(f, s, K, standard_brownian<d3,double,double,d3>(t));	
	double kappa = -black<>::binary<>(f, s, -k, standard_brownian<>(t)); // dv/dk
	ensure (fabs(V.derivative(1) - kappa) <= eps);


	const distribution_base<>* pd = new standard_brownian<>(t);
	p = black<>::value(f, s, -k, *pd);

	V = black<d3>::value(d3(f,1), d3(s), d3(k), standard_brownian<d3,d3,d3,d3>(d3(t)));
	ensure (V._(0) == c);
	ensure (fabs(V._(1) - delta) <= eps);
	ensure (fabs(V._(2) - gamma) <= eps);

	V = black<d3>::value(d3(f), d3(s,1), d3(k), standard_brownian<d3,d3,d3,d3>(d3(t)));
	ensure (V._(0) == c);
	ensure (fabs(V._(1) - vega) <= f*eps);

	V = black<d3>::value(d3(f), d3(s), d3(-k,1), standard_brownian<d3,d3,d3,d3>(d3(t)));
	ensure (V._(0) == c);
	ensure (fabs(V._(1) - kappa) <= eps);

	V = black<d3>::value(d3(f), d3(s), d3(k), standard_brownian<d3,d3,d3,d3>(d3(t,1)));
	ensure (V._(0) == c);
	ensure (fabs(V._(1) - theta) <= theta*eps);
/*
	double sigma;
	sigma = black<>::implied(f, p, -k, standard_brownian<d3,double,d3,d3>(t));
	double s_ = .4;
	p = black<>::value<>(f, s_, -k, standard_brownian<>(t));
	sigma = black<>::implied(f, p, -k, standard_brownian<d3,double,d3,d3>(t));
	c = black<>::value<>(f, s_, +k, standard_brownian<>(t));
	sigma = black<>::implied(f, c, k, standard_brownian<d3,double,d3,d3>(t));
*/
}

inline double 
timer(size_t n, const std::function<void(void)>& f)
{
	clock_t b = clock();
	while (n--) f();
	clock_t e = clock();

	return 1.*(e - b)/CLOCKS_PER_SEC;
}

void
fms_option_timing(void)
{
	double dt;
	double f = 100;
	double s = .2;
	double k = 100;
	double t = 0.25;

	size_t n = 1;
	
	dt = timer(n, [f,s,k,t](void) { black<>::value(f, s, k, standard_brownian<>(t)); });
	dt = dt;
	typedef dual::number<double,2> d2;
	dt = timer(n, [f,s,k,t](void) { black<d2>::value(d2(f), d2(s), d2(k), standard_brownian<d2,d2,d2,d2>(d2(t))); });
	dt = timer(n, [f,s,k,t](void) { black<d2>::value(d2(f), d2(s), d2(k), standard_brownian<d2,d2,d2,d2>(d2(t))); });
}
