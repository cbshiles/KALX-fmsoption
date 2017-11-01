// fms_normal_test.cpp - test normal distirbutions
#ifndef ensure
#include <cassert>
#define ensure(x) assert(x)
#endif

#include "daly.h"
using namespace fms;

template<class I, class Z>
void fms_normal_test_(void)
{
	ensure (0.5 == normal<I>::cdf<Z>(0));
}

void fms_normal_test(void)
{
	fms_normal_test_<daly, double>();
}