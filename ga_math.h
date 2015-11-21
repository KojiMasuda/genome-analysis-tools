#ifndef _GA_MATH_H_
#define _GA_MATH_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float ga_mean (float arr[], unsigned long n);
float ga_var (float arr[], unsigned long n);
float ga_covar (float arr1[], float arr2[], unsigned long n);
float ga_ustd (float arr[], unsigned long n);
float ga_t_table (unsigned long dof);

#endif
