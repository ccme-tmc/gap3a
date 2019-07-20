=======================
Programming Conventions
=======================

-------
Fortran
-------

^^^^^^
Styles
^^^^^^

* Strict ANSI Fortran90 should be used. Features marked as obsolescent in F90/95 should be avoided.
* Whenever possible, *free-form* source format should be used. In this case, code should start from the **1st** column. 
* Any code should be written in **lower-case** form.
* ``use`` statements should include the ``only`` option and the corresponding list of global variables used by the subroutine, unless all the variables in the module are used.
* Extra indentation of *2* columns should be added inside each loop level.
* Length of each line should be limited to *80* characters, using ``&`` character for line continuation.
* The line truncation should take at logical place, and the line continuing should take place such that 
  
  #. Procedures: parameters are aligned right after the left parenthese.
  #. Assignments: right values are aligned along the very first value after ``=``
  #. Strings: characters are aligned right after the left quotation mark.

^^^^^^^^^
Variables
^^^^^^^^^

* Declaration of double of complex variable should use the form ``datatype(N)``.
  ``double precision`` or ``double complex`` should be avoided.
* Global variables should be declared in modules.
* Subroutine parameters must be declared with ``intent``.
* Initialize all variables.
* Each variable should be declared separately.
* Fortran keywords should not be used as variable names.
* Use ``constant`` module to define and extract physical constants.
  Physical constants should never be hardwired into the executable portion of the code.

^^^^^^^^^^^^^^^^^^
Writing procedures
^^^^^^^^^^^^^^^^^^

* Use of ``implicit none`` is mandatory.
* **DRY**: avoid redundant or repeated code: check to see if the routine you need already exists before writing a new one.
* Every function or subroutine, no matter how small, should be in its own file named ``routine.f90``,
  where ``routine`` is the function or subroutine name.
  It is recommended that the routines are named so as to make their purpose apparent from the name alone.
* Check consistency between input parameters whenever possible.
* All internal units should be **atomic**, or clearly stated otherwise (with exception of WIEN2k output files used as input).
* All called procedures (intrinsic or external) within the subroutines should be explicitly declared.
* Keep usage of ``goto`` statements to a minimum, only when it is impossible to avoid.
  They should be used for exiting loops only and always point to a continue statement.

^^^^^^^^^^^^^^^^^^
Inputs and outputs
^^^^^^^^^^^^^^^^^^

* All units of input and output should be encapsulated in a derived type variable named as ``io``.
* The following information should be dumped to ``io%ou6``:

  #. System information
  #. Calculation results
  #. Warnings
  #. Errors

* The following information should be dumped to ``io%ou0``:
  
  #. Debugging information

* Output with negative unit should be suppressed.
* All I/O statements should appropriately contain the status specifier, ``err``, ``iostat`` and ``end``.
* Errors should be reported prior to termination with a short desription using the ``outerr`` subroutine.

^^^^^^^^^^^^^^^^^^^^^
Profiling and timings
^^^^^^^^^^^^^^^^^^^^^

* Timing of each functionality should utilize a variable of derived type ``timer`` and its methods,
  implemented in ``profiling`` module .

^^^^^^^^^^^^^^^^^^^^^^^^
Comments and documenting
^^^^^^^^^^^^^^^^^^^^^^^^

* API should be sufficiently documented in the source by using the reStructedText markup language,
  in order to utilize the documentation generator Sphinx and Sphinx-fortran extension.
* The function, basic logic and the algorithm used within should be described in the header,
  i.e. the lines after the declaration of procedure and before the first line of code.
* Purpose of each variable and each field of derived type should be described in a short comment on the same line.
* Subroutines should be "plentifully" commented. If you are not sure whether or not to add comment somewhere, **do it**.
* A tentative template of documenting subroutine is followed as below

  .. code-block:: fortran

     subroutine subname(a, b)
     !  <short description>
     !
     !  <Long description describing the function, equations and logic, as well as references>
     !
     !  The subroutine solves
     !
     !  .. math::
     !
     !     <latex math>
     !
     !  See also
     !
     !  <References>
     !
     !  Revision history:
     !  
     !  * YYYY.MM.DD <NAME> <Description of activity>
     !
     !(source code following)
       integer, intent(in)  :: a ! explaination of a
       real(8), intent(out) :: b ! explaination of b
       integer              :: i ! counter
     !(source code following)
     end subroutine subname

^^^^^^
Safety
^^^^^^

* Local allocatable arrays must be deallocated on exit of the routine to prevent memory leakage.
* Use corresponding generator and destructor subroutines to allocate memory to data within derived types.
* Each routine should terminate the program when given improper input, and report an error.

------
Python
------

#. All python scripts should work with Python 2.7 and 3.6 or newer versions.
#. Use "if-main" scheme for all Python scripts.

