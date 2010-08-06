/*
 Construct basis set of KL modes from scratch. Based on numerically solving 
 the Karhunen-Loève integral. Based on Dai (2001) and Fried (1978)
 
 See also Roddier (1990), Wang et al. (1978), Noll (1976) for further 
 reference.
 
 Compile with:
	g++ -O3 -Wall -L/sw/lib -I/sw/include -lgsl -lgslcblas util.cc a1a2.cc kernel.cc klbasis.cc -o klbasis
 
 (c) Tim van Werkhoven <t.i.m.vanwerkhoven@xs4all.nl> 2010
 
 This work is licensed under the Creative Commons Attribution-NonCommercial-
 ShareAlike 3.0 Netherlands License. To view a copy of this license, visit 
 http://creativecommons.org/licenses/by-nc-sa/3.0/nl/ or send a letter to 
 Creative Commons, 171 Second Street, Suite 300, San Francisco, California, 
 94105, USA.
*/

#include <stdio.h>
#include <getopt.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>

#include "a1a2.h"
#include "kernel.h"
#include "util.h"
#include "klbasis.h"

void show_version() {
	printf("klbasis (version 0.1, built %s %s)\n", __DATE__, __TIME__);
	printf("Copyright (c) 2010 Tim van Werkhoven <T.I.M.vanWerkhoven@xs4all.nl>\n");
	printf("\nFOAM comes with ABSOLUTELY NO WARRANTY. This is free software,\n"
		   "and you are welcome to redistribute it under certain conditions;\n"
		   "see the file COPYING for details.\n");
}

void show_clihelp(char *execname, bool error = false) {
	if (error)
		fprintf(stderr, "Try '%s --help' for more information.\n", execname);
	else {
		printf("Usage: %s [option]...\n\n", execname);
		printf("  -n, --ngrid=N        Use N gridpoints for eigenfunction solutions.\n"
			   "      --maxq=N         Solve KL functions up to q=N.\n"
			   "      --minq=N         Solve KL functions from q=N.\n"
			   "  -o, --order=N        Use interpolation order N.\n"
			   "      --cache          Cache intermediate results.\n"
			   "  -h, --help           Display this help message.\n"
			   "      --version        Display version information.\n\n");
		printf("Report bugs to Tim van Werkhoven <T.I.M.vanWerkhoven@xs4all.nl>.\n");
	}
}


int main(int argc, char *argv[]) {
	// Parse command line options
	int r, option_index = 0;
	
	int ngrid=253;		// Number of gridpoints to solve eigenfunctions on
	int order=4;		// Interpolation order to use
	int cache=0;		// Cache intermediate results
	int minQ=0;			// Minimum value of q to try
	int maxQ=2;			// Maximum value of q to try
	
	static struct option const long_options[] = {
		{"ngrid", required_argument, NULL, 'n'},
		{"maxq", required_argument, NULL, 3},
		{"minq", required_argument, NULL, 4},
		{"order", required_argument, NULL, 'o'},
		{"help", no_argument, NULL, 'h'},
		{"version", no_argument, NULL, 1},
		{"cache", no_argument, NULL, 2},
		{NULL, 0, NULL, 0}
	};
	
	while((r = getopt_long(argc, argv, "c:hv", long_options, &option_index)) != EOF) {
		switch(r) {
			case 0:
				break;
			case 'n':
				ngrid = (int) atoi(optarg);
				break;
			case 'o':
				order = (int) atoi(optarg);
				break;
			case 'h':
				show_clihelp(argv[0]);
				return -1;
			case 1:
				show_version();
				return -1;
			case 2:
				cache = 1;
				break;
			case 3:
				maxQ = atoi(optarg);
				break;
			case 4:
				minQ = atoi(optarg);
				break;
			case '?':
				show_clihelp(argv[0], true);
				return -1;
			default:
				break;
		}
	}
	
	// Check options
	if ((ngrid-1)/order != (ngrid-1)*1.0/order) {
		ngrid = (ngrid-1)/order * order + order + 1;
		fprintf(stderr, "Warning, ngrid should be equal to n*order + 1. Corrected ngrid to %d.\n", ngrid);
	}
	
	// Reserve memory for output KL modes and eigenvalues
	gsl_vector *kl_eigenv = gsl_vector_alloc(ALLOCSIZE);
	gsl_matrix *kl_eigenf = gsl_matrix_alloc(ngrid, ALLOCSIZE);
	int kl_nmodes = 0;
	int kl_nalloc = ALLOCSIZE;
	
	// Tabulate function A1(r) and calculate constant A2 (Dai Eqn. 23 and 24)
	printf("Calculating A1(r) and A2...\n");
	gsl_vector *A1 = gsl_vector_alloc(ngrid);
	double A2;

	if (cache && !gsl_restore_vector("klbasis-A1.gsl", A1)) {
		printf("Restored A1 from cache.\n");
	}
	else {
		for (int r=0; r<ngrid; r++)
			gsl_vector_set(A1, r, calc_A1((double) r /(ngrid-1.0)));
	}
	
	A2 = calc_A2();
	printf("Found A2 = %.16f\n", A2);
	
	if (cache)
		gsl_store_vector("klbasis-A1", A1, true, true, "%.16f");
	
	// Loop over all q, start at maxQ because that will determine the 
	// eigenvalue cutoff for further KL modes
	int q=maxQ;
	
	gsl_matrix *kernelS, *weightMat;
	gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(ngrid);
	gsl_vector * eigenV = gsl_vector_alloc(ngrid);
	gsl_matrix * eigenF = gsl_matrix_alloc(ngrid, ngrid);

	// Try to load from cache
	if (cache &&
			!gsl_restore_matrix(format("klbasis-eigenF_q=%d.gsl", q), eigenF) &&
			!gsl_restore_vector(format("klbasis-eigenV_q=%d.gsl", q), eigenV))
		printf("Restored first set of eigenfunction & -values for q=%d from cache.\n", q);
	// Otherwise calculate...
	else {
		printf("Calculating intergration kernel S_q=%d...\n", q);
		// Build kernel matrix S_q (Dai Eqn. 42)
		kernelS = calc_kernel(q, order, ngrid, &weightMat);
		printf("Solving eigensystem...\n");
		// Solve eigensystem into eigenfunctions & -values
		gsl_eigen_symmv (kernelS, eigenV, eigenF, w);
	}
	
	// Cache to disk
	if (cache) {
		gsl_store_matrix(format("klbasis-eigenF_q=%d", q), eigenF, true, true);
		gsl_store_vector(format("klbasis-eigenV_q=%d", q), eigenV, true, true);
	}
	
	double limit = gsl_vector_get(eigenV, 1);
	int p = 1;
	printf("Keeping first eigenvalue @ q=%d, p=%d, value=%.16f\n", q, p, limit);

	// Loop over all other values of q
	for (q=maxQ-1; q>=minQ; q--) {
		if (cache &&
			!gsl_restore_matrix(format("klbasis-eigenF_q=%d.gsl", q), eigenF) &&
			!gsl_restore_vector(format("klbasis-eigenV_q=%d.gsl", q), eigenV))
			printf("Restored eigenfunction & -values for q=%d from cache.\n", q);
		else {
			printf("Calculating intergration kernel S_q=%d...\n", q);
			kernelS = calc_kernel(q, order, ngrid, &weightMat);
			printf("Solving eigensystem...\n");
			gsl_eigen_symmv (kernelS, eigenV, eigenF, w);
		}
		
		if (cache) {
			gsl_store_matrix(format("klbasis-eigenF_q=%d", q), eigenF, true, true);
			gsl_store_vector(format("klbasis-eigenV_q=%d", q), eigenV, true, true);
		}
		
		// Keep all eigenfunctions & -values >= limit
		p=1;
		for (int i=1; gsl_vector_get(eigenV, i)>=limit; i++, p++) {
			double vec = gsl_vector_get(eigenV, i);
			if (vec < 0) {
				printf("Eigenvalue @ q=%d, p=%d, value=%g is negative! Abort!\n", q, p, vec);
				exit(-1);
			}
			printf("Keeping eigenvalue @ q=%d, p=%d, value=%.16f\n", q, p, vec);
		}
		
		// Normalize and convert to Karhunen-Loève functions
		
	}
	
	return 0;
}
