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

// Absolute and relative error for A1(r) and A2 integrals
#define ABSERR 0
#define RELERR 1e-9
// Allocate this much for eigenvalues and eigenfunctions at a time
#define ALLOCSIZE 50

void show_version();
void show_clihelp(char *execname, bool error);

#endif // HAVE_KLBASIS_H
