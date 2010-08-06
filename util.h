/*
 Some utility functions -- header file
 
 (c) Tim van Werkhoven <t.i.m.vanwerkhoven@xs4all.nl> 2010
 
 This work is licensed under the Creative Commons Attribution-NonCommercial-
 ShareAlike 3.0 Netherlands License. To view a copy of this license, visit 
 http://creativecommons.org/licenses/by-nc-sa/3.0/nl/ or send a letter to 
 Creative Commons, 171 Second Street, Suite 300, San Francisco, California, 
 94105, USA.
 */

#ifndef HAVE_UTIL_H
#define HAVE_UTIL_H

#include <cstdio>
#include <cstdarg>
#include <cstdlib>
#include <string>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

static inline std::string vformat(const char *format, va_list va) {
	char buf[4096];
	vsnprintf(buf, sizeof buf, format, va);
	return std::string(buf);
}

static inline std::string format(const char *format, ...) {
	va_list va;
	va_start(va, format);
	std::string result = vformat(format, va);
	va_end(va);
	return result;
}

int gsl_store_vector(const std::string file, const gsl_vector *vec, bool ascsv, bool asgsl, const std::string csvfmt = "%.16f");
int gsl_store_matrix(const std::string file, const gsl_matrix *mat, bool ascsv, bool asgsl, const std::string csvfmt = "%.16f");
int gsl_restore_vector(const std::string file, gsl_vector *vec);
int gsl_restore_matrix(const std::string file, gsl_matrix *mat);

#endif // HAVE_UTIL_H
