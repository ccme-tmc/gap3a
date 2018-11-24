# Development note

## 2018-10-24 (zmy)

  1. Backup old `make.inc`s. TODO: make.inc-Linux options

## 2018-11-03 (zmy)

  1. Remove duplicate `Lebedev-laikov.f90` in `src_util`
  2. Add `veclen.f90` and `vecprojlen.f90` in `src_util`
  3. Add auxiliary vectors and tensors for anisotropy in `dielmat` module
  4. Standardize directory structure.
  
## 2018-11-04 (zmy)

  1. Change version name to `3a`

## 2018-11-05 (zmy)

  1. Convert to `complex(8)` explictly instead of default `complex(4)` when using `cmplx` function
    - `freq_factor` in `bzint`
    - `sub_bzintq_0`
    - `freq_intpl_ac`
    - `scgw_herm_sxc`
    - `task_acont`
    - `acpatrd`, `calcacfreq`, `getsac`,`setsac` in `src_acfreq`
  2. Change default condition for `iop_coul` from 0 to -1 in `calcminm`,`calcmicm`,`calcminc`,**`calcmwm`**,**`calcselfx`**.
  3. Move anisotropy-related quantites to a new module file `anisotropy.f90`
  4. Create `ten_rvctrv` to calculate $q\cdot T\cdot q$ in `src_util`
  5. Move initialization of head in `calceps` before the `isp` loop
  6. In `calchead`, correctly calculate head at `q0_eps` with `ten_p_aniso` for `iop_aniso.ne.-1`

## 2018-11-06 (zmy)

  1. In `calceps`, 
    - correctly calculate wings at `q0_eps` with `vec_u_aniso` and `vec_t_aniso` for `iop_aniso.ne.-1`
    - add `coef_coul` to make the 4\\pi coefficient possible to be dimension-dependent (TODO later)
  2. Initialize cutoff length in `readingw`
  3. Create `coul_coef` to calculate |q|-dependent coefficient for head and wing calculations
  4. In `barcoul`, `smallq` for isotropic dieletric function, since `q0_eps` only specifies direction (TODO)
  5. Separate inversion of body `eps` from the inversion of whole dielectric matrix
  6. Try inverting the dielectric matrix by tensor A and vector a,b,u,t but failed. (TODO)

## 2018-11-07 (zmy)

  1. Fix inversion of dielectric matrix by disabling the use of `zhemv` in 
  calulating `epsw1` and `epsw2` when inverting the dielectric matrix for imaginary frequency.
  When using `zhemv`, calculating `head` with `epsw2+bw1` and `w2b+epsw1` gives different results,
  although they should be identical due to the analytic expression.

## 2018-11-08 (zmy)

  1. Fix the `zgemm` in calculating `ten_a_aniso`.
  2. Rearrange some `iop_aniso` if condition
  3. Calcualte dielectric matrix on `q0_sph` except for the body part. Memory corruption happens for `nq0>=14`

## 2018-11-09 (zmy)

  1. `im_g0` defaults to 1 in `barcoul` module. (12.273 to 12.275 for LiF-nk1 GW0 gap)
  2. Fix memory corruption by fixing bugs when using `DOPTS` in makefile as `FFLAGS`
  3. Use `ALLOCATABLE` instead of `POINTER` to store tensors and vectors in `anisotropy`
  4. Add `debug' macro in `make.inc` to easily use `DOPTS` for tests
  5. Update `README.md` in the root directory
  6. Add `smallq`,`smallq_div`,`qmax_q0` in `anisotropy` to define the proximity around Gamma point for integration
  7. Subroutine `init_smallq` to decide `qmax` along `q0_sph` on the defined `q0` region
  8. Subroutine `angint_eps_sph` in `anisotropy` to average over Gamma proximity. Done for head.

## 2018-11-10 (zmy)

  1. Subroutine `angint_eps_sph` for averaging wings and body over Gamma proximity.

## 2018-11-13 (zmy)

  1. Add time counting for anisotropy utilies (Not working for MPI, TODO)
  2. Explicitly write all used variables from `anisotropy` module in `calceps` subroutine
  3. Fix MPI bug by if condition in `calceps`

## 2018-11-18 (zmy)

  1. Test data should that only averaging the dielectric matrix does not fully remove the singularity at Gamma
  for anisotropic dielectric screening: the band gap depends linearly on the number of k point along z direction
  for hBN test case.
  2. Try to use `singc1ex,singc2ex` and `singc1co,singc2co` to represent the coefficient used to deal with q-1 
  and q-2 singularities in exchange and correlation self-energy, respectively. `singc1` and `singc2` are not 
  deprecated yet for the sake of stability. But they will be removed later.
  The reason why the diffferentiation is necessary is that to deal with anisotropy in dielectric function, 
  the auxiliary function Fs, which is used to calculate `singc1/2` only adds to the complexity to 
  integrate over q0, since the integration should be performed explicitly to account for anisotropy.
  However, they are still useful in calculating exchange self-energy.
  3. Figure out that `calcmwm` instead of `calc_mwm3` is used to calculate M*W*M for correlationi selfenergy.
  Hence non-analyticity for head should be considered in this part.(TODO)

## 2018-11-19 (zmy)

  1. Move the factor `4pi` in Lebedev-Laikov integration to weight `wt_q0_sph`.
  This adds missing `4pi` factor in `bodyinv` correction in `angint_eps_sph`

## 2018-11-20 (zmy)

  1. Check expansion on spherical harmonics

## 2018-11-21 (zmy)

  1. Remove body correction for the current stage.
  2. Use the scheme described as Eq.(45) in Friedrich, et al PRB 81,125102(2010).

## 2018-11-22 (zmy)

  1. Use the scheme described as Eq.(36) in Freysoldt, et al CPC 176,1(2007).

## 2018-11-23 (zmy)

  1. Merge code for generation of angular grids in `anisotropy` and `bzinteg`

