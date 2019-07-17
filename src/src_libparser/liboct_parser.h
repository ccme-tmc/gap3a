/*
 Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2, or (at your option)
 any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 02111-1307, USA.
*/

#ifndef _LIB_OCT_H
#define _LIB_OCT_H

#include "gsl_complex.h"

int parse_init(char *file_out);
int parse_input(char *file_in);
void parse_end();

int parse_isdef(char *name);

int parse_int(char *name, int def);
double parse_double(char *name, double def);
gsl_complex parse_complex(char *name, gsl_complex def);
char *parse_string(char *name, char *def);

int parse_block_n(char *name);
int parse_block_cols(char *name, int l);
int parse_block_int(char *name, int l, int col, int *r);
int parse_block_double(char *name, int l, int col, double *r);
int parse_block_complex(char *name, int l, int col, gsl_complex *r);
int parse_block_string(char *name, int l, int col, char **r);

/* from parse_exp.c */
typedef struct parse_result{
  union {
	  gsl_complex c;
		char *s;
	} value;
	enum {PR_CMPLX, PR_STR} type;
} parse_result;

int parse_exp(char *exp, parse_result *t);

void parse_putsym_int(char *s, int i);
void parse_putsym_double(char *s, double d);
void parse_putsym_complex(char *s, gsl_complex c);

#endif
