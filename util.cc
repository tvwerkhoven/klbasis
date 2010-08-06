/*
 Some utility functions.
 
 (c) Tim van Werkhoven <t.i.m.vanwerkhoven@xs4all.nl> 2010
 
 This work is licensed under the Creative Commons Attribution-NonCommercial-
 ShareAlike 3.0 Netherlands License. To view a copy of this license, visit 
 http://creativecommons.org/licenses/by-nc-sa/3.0/nl/ or send a letter to 
 Creative Commons, 171 Second Street, Suite 300, San Francisco, California, 
 94105, USA.
 */

#include <cstdio>
#include <cstdarg>
#include <cstdlib>
#include <string>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "util.h"

int gsl_store_vector(const std::string file, const gsl_vector *vec, bool ascsv, bool asgsl, const std::string csvfmt) {
	FILE *fd;
	int ret=0;
	// Store csv
	if (ascsv) {
		fd = fopen((file + ".csv").c_str(), "w");
		ret += gsl_vector_fprintf(fd, vec, csvfmt.c_str());
		fclose(fd);
	}
	// Store GSL
	if (asgsl) {
		fd = fopen((file + ".gsl").c_str(), "w");
		ret += gsl_vector_fwrite(fd, vec);
		fclose(fd);
	}
	return ret;
}

int gsl_store_matrix(const std::string file, const gsl_matrix *mat, bool ascsv, bool asgsl, const std::string csvfmt) {
	FILE *fd;
	int ret=0;
	// Store csv
	if (ascsv) {
		fd = fopen((file + ".csv").c_str(), "w");
		ret += gsl_matrix_fprintf(fd, mat, csvfmt.c_str());
		fclose(fd);
	}
	// Store gsl
	if (asgsl) {
		fd = fopen((file + ".gsl").c_str(), "w");
		ret += gsl_matrix_fwrite(fd, mat);
		fclose(fd);
	}
	return ret;
}

int gsl_restore_vector(const std::string file, gsl_vector *vec) {
	FILE *fd;
	int ret=0;
	
	fd = fopen(file.c_str(), "r");
	if (!fd) return -1;

	ret = gsl_vector_fread(fd, vec);
	fclose(fd);
	
	return ret;
}

int gsl_restore_matrix(const std::string file, gsl_matrix *mat) {
	FILE *fd;
	int ret=0;
	
	fd = fopen(file.c_str(), "r");
	if (!fd) return -1;
	
	ret = gsl_matrix_fread(fd, mat);
	fclose(fd);
	
	return ret;
}
