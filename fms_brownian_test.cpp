// fms_normal_test - test fms::normal functions
#ifndef ensure
#include <cassert>
#define ensure(x) assert(x)
#endif
#include "daly.h"
#include "brownian.h"

using namespace fms;

template<class I, class T>
void fms_brownian_test_(void)
{
	T t;
	t = brownian<daly,T,T,T>::kappa(.1,.2);
	ensure (t == .1*.1*.2/2);
}

void fms_brownian_test(void)
{
	fms_brownian_test_<daly, double>();
}