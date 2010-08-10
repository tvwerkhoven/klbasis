/*
 Header file for klbasis.cc
 
 (c) Tim van Werkhoven <t.i.m.vanwerkhoven@xs4all.nl> 2010
 
 This work is licensed under the Creative Commons Attribution-NonCommercial-
 ShareAlike 3.0 Netherlands License. To view a copy of this license, visit 
 http://creativecommons.org/licenses/by-nc-sa/3.0/nl/ or send a letter to 
 Creative Commons, 171 Second Street, Suite 300, San Francisco, California, 
 94105, USA.
 */

#ifndef HAVE_KLBASIS_H
#define HAVE_KLBASIS_H

#include <pthread.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

// Absolute and relative error for A1(r) and A2 integrals
#define ABSERR 0
#define RELERR 1e-9
// Allocate this much for eigenvalues and eigenfunctions at a time
#define ALLOCSIZE 5

struct _kl_config {
	int ngrid;			// Number of gridpoints to solve eigenfunctions on
	int order;			// Interpolation order to use
	int cache;			// Cache intermediate results
	int maxQ;			// Maximum value of q to try
	int minQ;			// Minimum value of q to try
	double limit;		// Cut-off value for eigenvalues
	int nthreads;		// Number of threads to use
	pthread_t *threads;
};

struct _kl_modes {
	gsl_vector *eigv;
	gsl_matrix *eigf;
	gsl_matrix_uint *pq;	
	int nm;
	int nalloc;
	pthread_mutex_t lock;	// Lock for writing to kl_modes
};

struct thr_info {
	int id;
	int minq;
	int maxq;
	struct _kl_modes *kl_modes;
	struct _kl_config *cfg;
};

void show_version();
void show_clihelp(char *execname, bool error);
void *thread_worker(void *arg);
int calc_kl(int q, struct _kl_config *cfg, struct _kl_modes* out);

#endif // HAVE_KLBASIS_H
