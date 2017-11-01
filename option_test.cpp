// main.cpp - Test option pricing routines
// Copyright (c) 2011 KALX, LLC. All rights reserved. No warranty is made.
#include <iostream>

void fms_normal_test(void);
void fms_brownian_test(void);
void fms_option_test(void);

int main()
{
	try {
		fms_normal_test();
		fms_brownian_test();
		fms_option_test();
	}
	catch (const std::exception& ex) {
		std::cerr << ex.what() << std::endl;

		return -1;
	}

	return 0;
}