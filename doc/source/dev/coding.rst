Programming Conventions
=======================

Fortran
-------

#. All internal units should be **atomic**, or clearly stated otherwise (with exception of WIEN2k output files used as input).
#. Formatting

   * Strict ANSI Fortran90 should be used. Features marked as obsolescent in F90/95 should be avoided.
   * Whenever possible, *free-form* source format should be used. In this case, code should start from the first column. 
   * Length of each line should be limited to *80* characters, using ``&`` character for line continuation. The line truncation should take at logical place. 
   * Any code should be written in *lower-case* form.
   * ``use`` statements should include the ``only`` option and the corresponding list of global variables used by the subroutine, unless all the variables in the module are used.

#. Indentation:

   * Extra indentation of *2* columns should be added inside each loop level.
   * For non-procedure, continued lines should indent with double the number of spaces that each block is indented.
   * For procedures, indentation should take place at 2 columns before the left parentheses, followed by ``& `` and then continued parameters.

#. Variables:

   * Declaration should be of the form ``datatype(N)``. ``double precision`` or ``double complex`` should be avoided.
   * Global variables should be declared in modules.
   * Use of ``implicit none`` is mandatory.
   * Each variable should be declared separately

#. Comments and documenting

   * Functions and subroutines should be documented 
   * Purpose of each variable should be described in a short comment on the same line.
   * Subroutines should be "plentifully" commented. If you are not sure whether or not to add comment somewhere, **do it**.

#. When coding procedures

   * *DRY*: avoid redundant or repeated code: check to see if the routine you need already exists before writing a new one.
   * Defensive programming is strongly recommended.
   * Every function or subroutine, no matter how small, should be in its own file named ``routine.f90``, where ``routine`` is the function or subroutine name.
     It is recommended that the routines are named so as to make their purpose apparent from the name alone.
   * Keep usage of ``goto`` statements to a minimum, only when it is impossible to avoid.
     They should be used for exiting loops only and always point to a continue statement.
   * All called procedures (intrinsic or external) within the subroutines should be explicitly declared.

#. Debugging

#. Safety

   * Local allocatable arrays must be deallocated on exit of the routine to prevent memory leakage.
   * Use corresponding generator and destructor subroutines to allocate memory to data within derived types.
   * Each routine should terminate the program when given improper input.
   * Report errors prior to termination with a short desription using the ``outerr`` subroutine.


Python
------

#. All python scripts should work with Python 2.7 and 3.6 or newer versions.
#. Use "if-main" scheme for all Python scripts.

