/*
 Header file for a1a2.cc
 
 (c) Tim van Werkhoven <t.i.m.vanwerkhoven@xs4all.nl> 2010
 
 This work is licensed under the Creative Commons Attribution-NonCommercial-
 ShareAlike 3.0 Netherlands License. To view a copy of this license, visit 
 http://creativecommons.org/licenses/by-nc-sa/3.0/nl/ or send a letter to 
 Creative Commons, 171 Second Street, Suite 300, San Francisco, California, 
 94105, USA.
 */

#ifndef HAVE_A1A2_H
#define HAVE_A1A2_H

double calc_A1(double r);
double calc_A1_f1(double rpp, void *params);
double calc_A1_f2(double thpp, void *params);

double calc_A2();
double calc_A2_f1(double r, void *params);

#endif // HAVE_A1A2_H
