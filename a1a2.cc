/*
 Calculate function A1(r) and constant A2 from Dai (2001).
 
 A1(r) = \Int_0^1 r``dr`` \Int_0^2{\pi} (r^2 + r``^2 - 2 r r`` \cos(\theta``) dtheta``
 A2 = 2 \Int_0^1 r A1(r) dr
 
 A2 should be approximately 2.991717027

 (c) Tim van Werkhoven <t.i.m.vanwerkhoven@xs4all.nl> 2010
 
 This work is licensed under the Creative Commons Attribution-NonCommercial-
 ShareAlike 3.0 Netherlands License. To view a copy of this license, visit 
 http://creativecommons.org/licenses/by-nc-sa/3.0/nl/ or send a letter to 
 Creative Commons, 171 Second Street, Suite 300, San Francisco, California, 
 94105, USA.
 */

#include <math.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_integration.h>

#include "a1a2.h"
#include "klbasis.h"

// Calculate A1(r)
double calc_A1(double r) {
	double result, error;
	static gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	static gsl_function F;
	F.function = &calc_A1_f1;
	F.params = &r;
	
	gsl_integration_qags (&F, 0, 1, ABSERR, RELERR, 1000,
						  w, &result, &error);
	
	return result;
}

// The first integral of A1(r)
double calc_A1_f1(double rpp, void *params) {
	double result, error;
	double r = ((double *) params)[0];
	double inparms[2] = {r, rpp};
	static gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	static gsl_function F;
	F.function = &calc_A1_f2;
	F.params = inparms;
	
	gsl_integration_qags (&F, 0, 2 * M_PI, ABSERR, RELERR, 1000,
						  w, &result, &error); 
	
	return rpp * result;
}

// The second integral of A1(r)
double calc_A1_f2(double thpp, void *params) {
	double r = ((double *) params)[0];
	double rpp = ((double *) params)[1];
	return pow(r*r + rpp*rpp - 2*r*rpp*gsl_sf_cos(thpp), 5./6.);
}


// Calculate A2
double calc_A2() {
	double result, error;
	
	static gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	static gsl_function F;
	F.function = &calc_A2_f1;
	
	gsl_integration_qags (&F, 0, 1, ABSERR, RELERR, 1000,
						  w, &result, &error);
	
	return 2*result;
}

// Integral function of A2
double calc_A2_f1(double r, void *params) {
	return r * calc_A1(r);
}

