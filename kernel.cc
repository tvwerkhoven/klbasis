/*
 Calculate the integration kernel S_q(i,j)

 (c) Tim van Werkhoven <t.i.m.vanwerkhoven@xs4all.nl> 2010
 
 This work is licensed under the Creative Commons Attribution-NonCommercial-
 ShareAlike 3.0 Netherlands License. To view a copy of this license, visit 
 http://creativecommons.org/licenses/by-nc-sa/3.0/nl/ or send a letter to 
 Creative Commons, 171 Second Street, Suite 300, San Francisco, California, 
 94105, USA.
*/

#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include "kernel.h"
#include "klbasis.h"

gsl_matrix *calc_kernel0(int order, int ngrid, double *A1, double A2, gsl_matrix **weightMat) {
	double tmp;
	// Allocate memory for kernel and weight matrix
	gsl_matrix *kernelS = gsl_matrix_alloc(ngrid, ngrid);
	
	// Construct weight matrix
	*weightMat = calc_weights(order, ngrid);
	
	// Calculate kernel matrix
	for (int i=0; i<ngrid; i++) {
		for (int j=i; j<ngrid; j++) {
			// Kernel is symmetric, so we only need to calculate one half
			tmp = calc_Sq(i*1./(ngrid-1), j*1./(ngrid-1), q);
			gsl_matrix_set(kernelS, i, j, tmp);
			gsl_matrix_set(kernelS, j, i, tmp);
		}
	}
	
	// Apply weight to matrix
	gsl_matrix_mul_elements(kernelS, *weightMat);
	
	return kernelS;
}

gsl_matrix *calc_kernelq(int q, int order, int ngrid, gsl_matrix **weightMat) {
	double tmp;
	// Allocate memory for kernel and weight matrix
	gsl_matrix *kernelS = gsl_matrix_alloc(ngrid, ngrid);
	
	// Construct weight matrix
	*weightMat = calc_weights(order, ngrid);
	
	// Calculate kernel matrix
	for (int i=0; i<ngrid; i++) {
		for (int j=i; j<ngrid; j++) {
			// Kernel is symmetric, so we only need to calculate one half
			tmp = calc_Sq(i*1./(ngrid-1), j*1./(ngrid-1), q);
			gsl_matrix_set(kernelS, i, j, tmp);
			gsl_matrix_set(kernelS, j, i, tmp);
		}
	}
	
	// Apply weight to matrix
	gsl_matrix_mul_elements(kernelS, *weightMat);
	
	return kernelS;
}

// Calculate weights to be applied to kernel matrix, i.e. sqrt(w_i * w_k) 
// from Dai (2001), Eqn. 41
gsl_matrix *calc_weights(int order, int ngrid) {
	double weights[order+1];

	gsl_matrix * weightMat = gsl_matrix_alloc(ngrid, ngrid);
	gsl_matrix * weightVec = gsl_matrix_alloc(ngrid, 1);
	
	if (newtoncotes_factors(order, weights))
		return NULL;
	
	// Multiply first and last weight by two to get the extended coefficients.
	weights[0] *= 2;
	weights[order] *= 2;
	
	// Fill vector with weights
	for (int i=0; i<ngrid; i++)
		gsl_matrix_set(weightVec, i, 0, sqrt(weights[i % order]));
	
	// Correct first and last terms
	gsl_matrix_set(weightVec, 0, 0, sqrt(weights[0]/2.));
	gsl_matrix_set(weightVec, ngrid-1, 0, sqrt(weights[0]/2.));
	
	// Multiply vector with itself to produce weight matrix
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, weightVec, weightVec, 0.0, weightMat);
	
	// Scale with stepsize
	gsl_matrix_scale(weightMat, 1./(ngrid-1));

	return weightMat;
}

// Calculate Newton-Cotes factors for a specific order. Note that these 
// coefficients are *not* the extended version.
int newtoncotes_factors(int order, double factors[]) {
	if (order == 1) {
		factors[0] = factors[order] = 0.5;
	}
	else if (order == 2) {
		factors[0] = factors[order] = 1./3;
		factors[1] = 1./4;
	}
	else if (order == 3) {
		factors[0] = factors[order] = 3./8;
		factors[1] = factors[order-1] = 9./8;
	}
	else if (order == 4) {
		factors[0] = factors[order] = 14./45;
		factors[1] = factors[order-1] = 64./45;
		factors[2] = 24./45;
	}
	else if (order == 5) {
		factors[0] = factors[order] = 5*19./288;
		factors[1] = factors[order-1] = 5*75./288;
		factors[2] = factors[order-2] = 5*50./288;
	}
	else if (order == 9) {
		factors[0] = factors[order] = 9*2857./89600;
		factors[1] = factors[order-1] = 9*15741./89600;
		factors[2] = factors[order-2] = 9*1080./89600;
		factors[3] = factors[order-3] = 9*19344./89600;
		factors[4] = factors[order-4] = 9*5778./89600;
	}
	else 
		return -1;
	
	return 0;
}

// Calculate S_q(r,r`) (Dai (2001), Eqn. 38)
double calc_Sq(double r, double rp, int q) {
	if (r == 0 || rp == 0) 
		return 0;
	
	// Setup integration parameters
//	static gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
//	static gsl_integration_qawo_table * qawo_table = 
//		gsl_integration_qawo_table_alloc(q, 2 * M_PI, GSL_INTEG_COSINE, 16);
	
	static gsl_integration_glfixed_table *gl_table = gsl_integration_glfixed_table_alloc (100);
	
	static gsl_function F;
	double result, error;
	
	// This should be a0 = 0.3450933461
	static const double a0 = pow(24./5 * gsl_sf_gamma(6./5), 5./6) / (pow(2., 5./3) * M_PI);
	
	// Calculate the integration part
	double inparms[3] = {r, rp, q};
//	F.function = &calc_Sq_f1;
	F.function = &calc_Sq_f1c;
	F.params = inparms;
	
	result = gsl_integration_glfixed(&F, 0.0, 2 * M_PI, gl_table);
//	gsl_integration_qawo (&F, 0.0, ABSERR, RELERR, 1000, w, qawo_table, 
//						  &result, &error);
	
	return -a0 * sqrt(rp) * sqrt(r) * result;
}

// Integration function for calc_Sq() without cosine factor
double calc_Sq_f1(double thpp, void *params) {
	double r = ((double *) params)[0];
	double rp = ((double *) params)[1];

#ifdef COS_LUT_SIZE
	static int have_table = 0;
	double cos_LUT[COS_LUT_SIZE];
	if (have_table == 0) {
		printf("Using cosine LUT for calc_Sq_f1(), this should be calculated only once...\n");
		for (int i=0; i<COS_LUT_SIZE; i++)
			cos_LUT[i] = cos(2.*M_PI*i/(COS_LUT_SIZE-1.));
		have_table = 1;
	}
	return pow(r*r + rp*rp - 2*r*rp*cos_LUT[(int) (COS_LUT_SIZE * thpp/(2.*M_PI))], 5./6);
#else	
	return pow(r*r + rp*rp - 2*r*rp*cos(thpp), 5./6);
#endif
}

// Integration function for calc_Sq() with cosine factor
double calc_Sq_f1c(double thpp, void *params) {
	double r = ((double *) params)[0];
	double rp = ((double *) params)[1];
	double q = ((double *) params)[2];

#ifdef COS_LUT_SIZE
	static int have_table = 0;
	static double cos_LUT[COS_LUT_SIZE];
	if (have_table == 0) {
		printf("Using cosine LUT for calc_Sq_f1c(), this should be calculated only once...\n");
		for (int i=0; i<COS_LUT_SIZE; i++)
			cos_LUT[i] = cos(2.*M_PI*i/(COS_LUT_SIZE-1.));
		have_table = 1;
	}
	return pow(r*r + rp*rp - 2*r*rp*cos_LUT[(int) (COS_LUT_SIZE*thpp/(2.*M_PI))], 5./6) * cos_LUT[(int) (COS_LUT_SIZE * ((q * thpp)/(2.*M_PI))) % COS_LUT_SIZE];
#else	
	return pow(r*r + rp*rp - 2*r*rp*cos(thpp), 5./6) * cos(q * thpp);
#endif
}
