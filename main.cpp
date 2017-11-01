// main.cpp - Test option pricing routines
// Copyright (c) 2011 KALX, LLC. All rights reserved. No warranty is made.
#include <iostream>

void fms_test_option(void);
void fms_test_brownian(void);
void fms_option_timing(void);
void fms_test_option_ene(void);

int main()
{
	try {
		fms_test_option();

		fms_test_brownian();
		fms_option_timing();
		fms_test_option_ene();
	}
	catch (const std::exception& ex) {
		std::cerr << ex.what() << std::endl;

		return -1;
	}

	return 0;
}