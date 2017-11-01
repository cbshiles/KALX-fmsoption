// fms_option_test.cpp - test option valuation
#include "../fmsdual/dual.h"
#include "daly.h"
#include "brownian.h"
#include "option.h"

using namespace fms;

template<class I, class T>
void fms_option_test_(void)
{
	T f = 100, s = 0.2, k = 100, t = 0.25, v1, v2;

	v1 = option::value(f, s, k, t, brownian<I>());
	v2 = option::call(f, s, k, t, brownian<I>());
	ensure (v1 == v2);

	v2 = option::put(f, s, k, t, brownian<I>());
	ensure (v1 == v2);

	v1 = option::value(f, s, -k, t, brownian<I>());
	ensure (v1 == v2);

	k = 90;
	v1 = option::call(f, s, k, t, brownian<I>());
	v2 = option::put(f, s, k, t, brownian<I>());
	ensure (f - k == v1 - v2);

	v2 = option::call(f, s, k, 2*t, brownian<I>());
	ensure (v2 > v1);

	v2 = option::call(f, s, k, t/2, brownian<I>());
	ensure (v2 < v1);

	dual::number<T,2> F(f,1), V;
	V = option::call(F, s, k, t, brownian<I,dual::number<T,2>>());
}

void fms_option_test(void)
{
	fms_option_test_<daly, double>();
}