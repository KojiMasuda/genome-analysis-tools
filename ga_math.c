#include "ga_math.h"

#define LOG(m) \
  fprintf(stderr, \
  "%s:line%d:%s(): " m "\n", \
  __FILE__, __LINE__, __FUNCTION__)

const float t_table[][2] = {{1, 12.70620},{2, 4.302653},{3, 3.182446},{4, 2.776445},{5, 2.570582},{6, 2.446912},{7, 2.364624},{8, 2.306004},{9, 2.262157},{10, 2.228139},{11, 2.200985},{12, 2.178813},{13, 2.160369},{14, 2.144787},{15, 2.131450},{16, 2.119905},{17, 2.109816},{18, 2.100922},{19, 2.093024},{20, 2.085963},{21, 2.079614},{22, 2.073873},{23, 2.068658},{24, 2.063899},{25, 2.059539},{26, 2.055529},{27, 2.051831},{28, 2.048407},{29, 2.045230},{30, 2.042272},{31, 2.039513},{32, 2.036933},{33, 2.034515},{34, 2.032245},{35, 2.030108},{36, 2.028094},{37, 2.026192},{38, 2.024394},{39, 2.022691},{40, 2.021075},{41, 2.019541},{42, 2.018082},{43, 2.016692},{44, 2.015368},{45, 2.014103},{46, 2.012896},{47, 2.011741},{48, 2.010635},{49, 2.009575},{50, 2.008559},{51, 2.007584},{52, 2.006647},{53, 2.005746},{54, 2.004879},{55, 2.004045},{56, 2.003241},{57, 2.002465},{58, 2.001717},{59, 2.000995},{60, 2.000298},{61, 1.999624},{62, 1.998972},{63, 1.998341},{64, 1.997730},{65, 1.997138},{66, 1.996564},{67, 1.996008},{68, 1.995469},{69, 1.994945},{70, 1.994437},{72, 1.993464},{74, 1.992543},{76, 1.991673},{78, 1.990847},{80, 1.990063},{85, 1.988268},{90, 1.986675},{95, 1.985251},{100, 1.983972},{110, 1.981765},{120, 1.979930},{130, 1.978380},{140, 1.977054},{150, 1.975905},{160, 1.974902},{170, 1.974017},{180, 1.973231},{200, 1.971896},{250, 1.969498},{300, 1.967903},{400, 1.965912},{500, 1.964720},{1000, 1.962339},{5000, 1.960439},{10000, 1.960201},{50000, 1.960011},{100000, 1.959988},{500000, 1.959969},{1000000, 1.959966}}; //97.5 percent point for t-distribution and dof

float ga_mean (float arr[], unsigned long n);
float ga_var (float arr[], unsigned long n);
float ga_covar (float arr1[], float arr2[], unsigned long n);
float ga_ustd (float arr[], unsigned long n);
float ga_t_table (unsigned long dof);

/*
 * This returns mean value from array.
 * arr[]: array
 * n    : length of array
 */
float ga_mean (float arr[], unsigned long n){
  unsigned long i;
  float sum = 0.0;

  for (i = 0; i < n; i++) {
    sum += arr[i];
  }

  return sum / n;
}

/*
 * This returns variance from array.
 * arr: array
 * n  : length of array
 */
float ga_var (float arr[], unsigned long n){
  unsigned long i;
  float dev = 0.0, ave;

  ave = ga_mean(arr, n);
  for (i = 0; i < n; i++) {
    dev += (ave - arr[i]) * (ave - arr[i]);
  }
  return dev / n ;
}

/*
 * This returns covariance from two arrays.
 * arr1: array1
 * arr2: array2
 * n   : length of array
 */
float ga_covar (float arr1[], float arr2[], unsigned long n){
  unsigned long i;
  float dev = 0.0, ave1, ave2;

  ave1 = ga_mean(arr1, n);
  ave2 = ga_mean(arr2, n);
  for (i = 0; i < n; i++) {
    dev += (ave1 - arr1[i]) * (ave2 - arr2[i]);
  }
  return dev / n ;
}

/*
 * This returns unbiased standard deviation from array.
 * arr: array
 * n  : length of array
 */
float ga_ustd (float arr[], unsigned long n){
  unsigned long i;
  float dev = 0.0, ave;
  if (n < 2) {
    LOG("The array length is less than two.");
    return 0.0;
  }

  ave = ga_mean(arr, n);
  
  for (i = 0; i < n; i++) {
    dev += (ave - arr[i]) * (ave - arr[i]);
  }
  return sqrt(dev / (n - 1));
}

/*
 * This returns 97.5 percent point for t-distribution with specified degree of freedum.
 * This doesn't calculate anything, just returns the record of percent points.
 * dof: degree of freedum
 */
float ga_t_table (unsigned long dof) { //97.5 percent point for t-distribution
  int i = 0;
  if (dof < 1) {
    LOG("Degree of freedum must be more than or equal to one.");
    return 0;
  }
  else if (dof > 1000000) {
    return 1.959964;
  }
  else if (dof == 1000000) {
    return 1.959966;
  }
  while(1){
    if (dof >= t_table[i][0] && dof < t_table[i+1][0]) break;
    i++;
  }
  return t_table[i][1];
}
