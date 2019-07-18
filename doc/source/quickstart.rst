===========
Quick Start
===========

--------------
GW Calculation
--------------

^^^^^^^^^^^^^^
Basic Workflow
^^^^^^^^^^^^^^

The basic workflow for the GW calculation is as follows.

1. Run a full SCF calculation using WIEN2K.

2. Prepare GW input files by running the shell script, ``init_gap3a``

   .. code-block:: shell

      $ init_gap<ver> -f casename -d gwdir -nkp nkp

   where

      * "casename": the case name of the WIEN2k calculation
      * "gwdir": the working diretory for gw calculation. It is recommended that gwdir be different from the working directory for SCF.
      * "nkp": the number of k-points to use in the GW.
   
   More detailed information about the options for ``init_gap3a`` can be got by running

   .. code-block:: shell

      $ init_gap<ver> -h 

3. Modify "gwdir/gw.inp" if necessary 

4. Start the GW calculation sequentially by running in the command line

   .. code-block:: shell

      $ gap<ver>.x

   or parallely by

   .. code-block:: shell

      $ mpirun -np gap<ver>-mpi.x

   within the "gwdir". 

5. To restart a GW calculation, set in the "gw.inp"

   .. code-block:: 

      Restart = T

   otherwise the calculation always starts from the beginning, no matter whether "casename.selfx/c" exists or not.

