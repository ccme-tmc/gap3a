include ../make.inc
LIB= ../libparser.a

OBJS = gsl_complex_math.o  parse.o parse_exp.o parser_f.o string_f.o symbols.o  liboct_parser.o

$(LIB):$(OBJS)
	$(AR) $(ARFLAGS) $(LIB) $(OBJS) 

liboct_parser.o:liboct_parser.f90
	$(FC) $(FFLAGS) -c liboct_parser.f90
gsl_complex_math.o: gsl_complex.h gsl_complex_math.h
	$(CC) $(CFLAGS) -c gsl_complex_math.c

parse.o: parse.c liboct_parser.h symbols.h
	$(CC) $(CFLAGS) -c parse.c

parse_exp.o: parse_exp.c gsl_complex.h gsl_complex_math.h liboct_parser.h symbols.h 
	$(CC) $(CFLAGS) -c parse_exp.c

parser_f.o: parser_f.c symbols.h liboct_parser.h string_f.h
	$(CC) $(CFLAGS) -c parser_f.c

string_f.o: string_f.c string_f.h
	$(CC) $(CFLAGS) -c string_f.c

symbols.o: symbols.c gsl_complex_math.h symbols.h
	$(CC) $(CFLAGS) -c symbols.c

#grammar.c: grammar.y
#	$(YACC) grammar.y
#	mv y.tab.c grammar.c

clean:
	rm -f $(OBJS) *.mod

cleanall:
	rm -f $(OBJS) *.mod $(LIB) 

