// toption.cpp
// Copyright (c) 2011 KALX, LLC. All rights reserved. No warranty is made.
#include <ctime>
#include <functional>
#include "../fmsdual/dual.h"
#include "option.h"
#include "standard_brownian.h"
#include "gamma.h"

using namespace option;

template<class D>
void
fms_test_option_(const D& t)
{
	double eps = std::numeric_limits<double>::epsilon();

	double f = 100;
	double s = 0.2;
	double k = 100;
	double p = black<>::value(f, s, -k, t);
	double c = black<>::value(f, s, +k, t);

	ensure (f - k == c - p);
}

void
fms_test_option(void)
{
	double t = 1;
	double a = 1, b = 1;

	fms_test_option_<>(standard_brownian<>(t));
//	fms_test_option_<>(gamma<>(t, a, b));
}

