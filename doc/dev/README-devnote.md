# Development note

## 2018-10-24 (zmy)

  1. Backup old `make.inc`s. TODO: make.inc-Linux options

## 2018-11-03 (zmy)

  1. Remove `Lebedev-laikov.f90` in `src_util` due to duplication
  2. Add `veclen.f90` and `vecprojlen.f90` in `src_util`
  3. Add auxiliary vectors and tensors for anisotropy in `dielmat.f90`
  4. Standarize directory structure.
  
## 2018-11-04 (zmy)

  1. Change version name to `3.0a`

## 2018-11-05 (zmy)

  1. convert to `complex(8)` instead of default `complex(4)` when using `cmplx` function
    - `freq_factor` in `bzint`
    - `sub_bzintq_0`
    - `freq_intpl_ac`
    - `scgw_herm_sxc`
    - `task_acont`
    - `acpatrd`, `calcacfreq`, `getsac`,`setsac` in `src_acfreq`
  2. change default condition for `iop_coul` from 0 to -1 in `calcminm`,`calcmicm`,`calcminc`,**`calcmwm`**,**`calcselfx`**.
  3. move anisotropy-related quantites to a new module file `anisotropy.f90`
  4. create `ten_rvctrv` to calculate $q\cdot T\cdot q$ in `src_util`
  5. move initialization of head in `calceps` before the `isp` loop
  6. In `calchead`, correctly calculate head at `q0_eps` with `ten_p_aniso` for `iop_aniso.ne.-1`

## 2018-11-06 (zmy)

  1. In `calceps`, 
    - correctly calculate wings at `q0_eps` with `vec_u_aniso` and `vec_t_aniso` for `iop_aniso.ne.-1`
    - add `coef_coul` to make the 4\\pi coefficient possible to be dimension-dependent (TODO later)
  2. Initialize cutoff length in `readingw`
  3. create `coul_coef` to calculate |q|-dependent coefficient for head and wing calculations
  4. In `barcoul`, `smallq` for isotropic dieletric function, since `q0_eps` only specifies direction (TODO)
  5. Separate inversion of body `eps` from the inversion of whole dielectric matrix
  6. Try inverting the dielectric matrix by tensor A and vector a,b,u,t but failed. (TODO)

## 2018-11-07 (zmy)

  1. Disable the use of `zhemv` in calulating `epsw1` and `epsw2` when inverting the dielectric matrix for imaginary frequency.
     When using `zhemv`, calculating `head` with `epsw2+bw1` and `w2b+epsw1` gives different results,
	 although they should be identical due to the analytic expression.

## 2018-11-08 (zmy)

  1. Fix the `zgemm` in calculating `ten_a_aniso`.
  2. Rearrange some `iop_aniso` if condition
  3. Calcualte dielectric matrix on `q0_sph` except for the body part.
     Memory corruption happens for `nq0>=14`. 
