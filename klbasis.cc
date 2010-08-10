/*
 Construct basis set of KL modes from scratch. Based on numerically solving 
 the Karhunen-Loève integral. Based on Dai (2001) and Fried (1978)
 
 See also Roddier (1990), Wang et al. (1978), Noll (1976) for further 
 reference.
 
 (c) Tim van Werkhoven <t.i.m.vanwerkhoven@xs4all.nl> 2010
 
 This work is licensed under the Creative Commons Attribution-NonCommercial-
 ShareAlike 3.0 Netherlands License. To view a copy of this license, visit 
 http://creativecommons.org/licenses/by-nc-sa/3.0/nl/ or send a letter to 
 Creative Commons, 171 Second Street, Suite 300, San Francisco, California, 
 94105, USA.
*/

#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <pthread.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_permutation.h>

#include "a1a2.h"
#include "kernel.h"
#include "util.h"
#include "klbasis.h"

void show_version() {
	printf("%s (version %s, built %s %s)\n", PACKAGE_NAME, PACKAGE_VERSION, __DATE__, __TIME__);
	printf("Copyright (c) 2010 Tim van Werkhoven <%s>\n", PACKAGE_BUGREPORT);
	printf("\n%s comes with ABSOLUTELY NO WARRANTY. This is free software,\n"
		   "and you are welcome to redistribute it under certain conditions;\n"
		   "see the file COPYING for details.\n", PACKAGE_NAME);
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
			   "  -l, --limit=L        Use this limit as eigenvalue cutoff.\n"
			   "  -t, --threads=N      Number of threads to use.\n"
			   "      --cache          Cache intermediate results or re-use them.\n"
			   "  -h, --help           Display this help message.\n"
			   "      --version        Display version information.\n\n");
		printf("Report bugs to Tim van Werkhoven <%s>.\n", PACKAGE_BUGREPORT);
	}
}

void *thread_worker(void *arg) {
	thr_info *info = (thr_info *) arg;
	int id = info->id;
	struct _kl_modes *kl_modes = info->kl_modes;
	struct _kl_config *cfg = info->cfg;
	
	// Now calculate the rest, loop over all other values of q
	for (int q=id; ; q += cfg->nthreads) {
		printf("Thread %d working on Q=%d.\n", id, q);
		int n = calc_kl(q, cfg, kl_modes);
		
		// If maxQ is set, stop when we reach this q
		if (cfg->maxQ != -1 && q>=cfg->maxQ)
			break;
		// If no more eigenvalues are above the cutoff limit, we're done
		else if (n == 0)
			break;
	}
	printf("Thread %d stopped @ Q=%d.\n", id, q);
	return NULL;
}

int calc_kl(int q, struct _kl_config *cfg, struct _kl_modes* out) {
	//printf("Calculating KL modes for q=%d. order=%d, ngrid=%d\n", q, cfg->order, cfg->ngrid);
	// TODO: Update ngrid dynamically
	int p=0;
	gsl_matrix *kernelS, *weightMat;	// Integration kernel to be inverted
	gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(cfg->ngrid);
	gsl_vector *tmp_col = gsl_vector_alloc(cfg->ngrid);
	gsl_vector *eigenV = gsl_vector_alloc(cfg->ngrid); // Store eigenvectors for specific Q here
	gsl_matrix *eigenF = gsl_matrix_alloc(cfg->ngrid, cfg->ngrid); // Store eigenfunctions for specific Q here
	
	// Try to load from cache
	if (cfg->cache &&
		!gsl_restore_matrix(format("klcache-eigenF_q=%d.gsl", q), eigenF) &&
		!gsl_restore_vector(format("klcache-eigenV_q=%d.gsl", q), eigenV))
		printf("Restored first set of eigenfunction & -values for q=%d from cache.\n", q);
	// Otherwise calculate...
	else {
		//printf("Calculating intergration kernel S_q=%d...\n", q);
		// Build kernel matrix S_q (Dai Eqn. 42)
		kernelS = calc_kernel(q, cfg->order, cfg->ngrid, &weightMat);
		//printf("Solving eigensystem...\n");
		// Solve eigensystem into eigenfunctions & -values
		gsl_eigen_symmv (kernelS, eigenV, eigenF, w);
	}
	
	// First set gets special treatment
	if (cfg->limit == 0.0) {
		cfg->limit = gsl_vector_get(eigenV, 1);
		p++;
		printf("Keeping first eigenvalue @ q=%d, p=%d, value=%.16f\n", q, p, cfg->limit);
		
		if (out->nm == out->nalloc) {
			fprintf(stderr, "Please increase ALLOCSIZE...\n");
			exit(-1);
		}

		gsl_vector_set(out->eigv, out->nm, cfg->limit);
		gsl_matrix_get_col(tmp_col, eigenF, 1);
		gsl_matrix_set_col(out->eigf, out->nm, tmp_col);
		gsl_matrix_uint_set(out->pq, 0, out->nm, p);
		gsl_matrix_uint_set(out->pq, 1, out->nm, q);
		out->nm++;
	}
	// Store regular modes & values here
	else {
		for (int i=1; i<cfg->ngrid && gsl_vector_get(eigenV, i)>=cfg->limit; i++) {
			p++;
			double vec = gsl_vector_get(eigenV, i);
			if (vec < 0) {
				fprintf(stderr, "Eigenvalue @ q=%d, p=%d, value=%g is negative! Abort!\n", q, p, vec);
				return -1;
			}
			//printf("Keeping eigenvalue @ q=%d, p=%d, value=%.16f\n", q, p, vec);
			
			pthread_mutex_lock(&(out->lock));
			if (out->nm == out->nalloc) {
				// Allocate new bigger memory block
				gsl_vector *tmp_eigv = gsl_vector_calloc(out->nalloc + ALLOCSIZE);
				gsl_matrix *tmp_eigf = gsl_matrix_alloc(cfg->ngrid, out->nalloc + ALLOCSIZE);
				gsl_matrix_uint *tmp_pq = gsl_matrix_uint_alloc(2, out->nalloc + ALLOCSIZE);
				out->nalloc += ALLOCSIZE;
				
				// Copy from old to new
				for (int i=0; i<out->nm; i++) {
					gsl_vector_set(tmp_eigv, i, gsl_vector_get(out->eigv, i));
					gsl_matrix_get_col(tmp_col, out->eigf, i);
					gsl_matrix_set_col(tmp_eigf, i, tmp_col);
					int p = gsl_matrix_uint_get(out->pq, 0, i);
					int q = gsl_matrix_uint_get(out->pq, 1, i);
					gsl_matrix_uint_set(tmp_pq, 0, i, p);
					gsl_matrix_uint_set(tmp_pq, 1, i, q);
				}
				
				gsl_vector_free(out->eigv);
				gsl_matrix_free(out->eigf);
				gsl_matrix_uint_free(out->pq);
				out->eigv = tmp_eigv;
				out->eigf = tmp_eigf;
				out->pq = tmp_pq;
			}
			gsl_vector_set(out->eigv, out->nm, vec);
			gsl_matrix_get_col(tmp_col, eigenF, i);
			gsl_matrix_set_col(out->eigf, out->nm, tmp_col);
			gsl_matrix_uint_set(out->pq, 0, out->nm, p);
			gsl_matrix_uint_set(out->pq, 1, out->nm, q);
			out->nm++;
			pthread_mutex_unlock(&(out->lock));
		}
	}
	
	// TODO: Normalize and convert to Karhunen-Loève functions
	
	// Cache to disk
	if (cfg->cache) {
		gsl_store_matrix(format("klcache-eigenF_q=%d", q), eigenF, true, true);
		gsl_store_vector(format("klcache-eigenV_q=%d", q), eigenV, true, true);
	}
	
	return p;
}


int main(int argc, char *argv[]) {
	// Parse command line options
	int r, option_index = 0;
	struct _kl_config cfg;
	
	cfg.ngrid = 253;
	cfg.order = 4;
	cfg.cache = 0;
	cfg.minQ = 0;
	cfg.maxQ = -1;
	cfg.limit = 0.0;
	cfg.nthreads = 2;

	static struct option const long_options[] = {
		{"ngrid", required_argument, NULL, 'n'},
		{"threads", required_argument, NULL, 't'},
		{"maxq", required_argument, NULL, 3},
		{"minq", required_argument, NULL, 4},
		{"order", required_argument, NULL, 'o'},
		{"limit", required_argument, NULL, 'l'},
		{"help", no_argument, NULL, 'h'},
		{"version", no_argument, NULL, 1},
		{"cache", no_argument, NULL, 2},
		{NULL, 0, NULL, 0}
	};
	
	while((r = getopt_long(argc, argv, "c:hvt:n:o:l:", long_options, &option_index)) != EOF) {
		switch(r) {
			case 0:
				break;
			case 'n':
				cfg.ngrid = (int) atoi(optarg);
				break;
			case 'o':
				cfg.order = (int) atoi(optarg);
				break;
			case 'h':
				show_clihelp(argv[0]);
				return -1;
			case 1:
				show_version();
				return -1;
			case 2:
				cfg.cache = 1;
				break;
			case 't':
				cfg.nthreads = atoi(optarg);
				break;
			case 'l':
				cfg.limit = atof(optarg);
				break;
			case 3:
				cfg.maxQ = atoi(optarg);
				break;
			case 4:
				cfg.minQ = atoi(optarg);
				break;
			case '?':
				show_clihelp(argv[0], true);
				return -1;
			default:
				break;
		}
	}
	
	// Check options
	if (cfg.maxQ == -1 && cfg.limit == 0) {
		fprintf(stderr, "Error, need either maxQ set or a limit of eigenvalues to use.\n");
		exit(-1);
	}
	else if (cfg.maxQ != -1 && cfg.limit != 0) {
		fprintf(stderr, "Error, cannot use maxQ and limit simultaneously.\n");
		exit(-1);
	}
	
	if ((cfg.ngrid-1)/cfg.order != 
			(cfg.ngrid-1)*1.0/cfg.order) {
		cfg.ngrid = (cfg.ngrid-1)/cfg.order * cfg.order + cfg.order + 1;
		fprintf(stderr, "Warning, ngrid should be equal to n*order + 1. Corrected ngrid to %d.\n", cfg.ngrid);
	}
	
	printf("%s starting. ngrid=%d, order=%d, cache=%d, q=%d--%d, limit=%g.\n", \
		   PACKAGE_NAME, cfg.ngrid, cfg.order, cfg.cache, cfg.minQ, cfg.maxQ, cfg.limit);
	
	// Reserve memory for output KL modes and eigenvalues
	struct _kl_modes kl_modes;
	
	kl_modes.eigv = gsl_vector_calloc(ALLOCSIZE);
	kl_modes.eigf = gsl_matrix_alloc(cfg.ngrid, ALLOCSIZE);
	kl_modes.pq = gsl_matrix_uint_alloc(2, ALLOCSIZE);
	kl_modes.nm = 0;
	kl_modes.nalloc = ALLOCSIZE;
	pthread_mutex_init(&(kl_modes.lock), NULL);
	
	if (cfg.maxQ != -1) {
		// If maxQ is set, calculate KL modes once for this Q to determine 
		// the cutoff value (limit) for eigenvalues
		calc_kl(cfg.maxQ, &cfg, &kl_modes);
		cfg.maxQ--;
	}
	
	// Setup threads
	struct thr_info threads[cfg.nthreads];
	cfg.threads = (pthread_t *) malloc((cfg.nthreads) * (sizeof(pthread_t)));
	
	for (int t=0; t<cfg.nthreads; t++) {
		threads[t].id = t;
		threads[t].cfg = &cfg;
		threads[t].kl_modes = &kl_modes;

		pthread_create(&cfg.threads[t], NULL, thread_worker, (void *)(threads+t));
	}
	
	// Wait for threads here
	for (int t=0; t<cfg.nthreads; t++) {
		pthread_join(cfg.threads[t], NULL);
	}
	
	// Sort vectors
	gsl_permutation * p = gsl_permutation_alloc(kl_modes.nalloc);
	gsl_sort_vector_index (p, kl_modes.eigv);
	
	// Print list of eigenvalues in descending order
	for (int n=kl_modes.nalloc-1; n>=(kl_modes.nalloc-kl_modes.nm); n--) {
		int idx = gsl_permutation_get(p, n);
		unsigned int p = gsl_matrix_uint_get(kl_modes.pq, 0, idx);
		unsigned int q =  gsl_matrix_uint_get(kl_modes.pq, 1, idx);
		double val = gsl_vector_get(kl_modes.eigv, idx);
		printf("Mode p=%d, q=%d:  %g\n", p, q, val);
		if (q != 0)
			printf("Mode p=%d, q=%d: %g\n", p, -q, val);
			
	}
	
	// Store final eigenfunctions and eigenvalues
	printf("Storing results to disk.\n");
	
	gsl_store_matrix(format("klbasis-eigenmodes", cfg.maxQ), kl_modes.eigf, true, true);
	gsl_store_vector(format("klbasis-eigenvalues", cfg.maxQ), kl_modes.eigv, true, true);

	return 0;
}
