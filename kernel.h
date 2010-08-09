/*
 kernel.cc header file
 
 (c) Tim van Werkhoven <t.i.m.vanwerkhoven@xs4all.nl> 2010
 
 This work is licensed under the Creative Commons Attribution-NonCommercial-
 ShareAlike 3.0 Netherlands License. To view a copy of this license, visit 
 http://creativecommons.org/licenses/by-nc-sa/3.0/nl/ or send a letter to 
 Creative Commons, 171 Second Street, Suite 300, San Francisco, California, 
 94105, USA.
 */

#ifndef HAVE_KERNEL_H
#define HAVE_KERNEL_H

#include <gsl/gsl_matrix.h>

gsl_matrix *calc_kernel(int q, int order, int ngrid, gsl_matrix **weightMat);
gsl_matrix *calc_weights(int order, int ngrid);
int newtoncotes_factors(int order, double factors[]);
double calc_Sq(int q, int i, int j, int ngrid, double *A1, double A2);
double calc_Sq_f1(double thpp, void *params);
double calc_Sq_f1c(double thpp, void *params);

#endif // HAVE_KERNEL_H
