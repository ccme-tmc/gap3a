# Development note

## 2018.10.24 (zmy)
  1. Backup old `make.inc`s. TODO: make.inc-Linux options

## 2018.11.3 (zmy)
  1. Remove `Lebedev-laikov.f90` in `src_util` due to duplication
  2. Add `veclen.f90` and `vecprojlen.f90` in `src_util`
  3. Add auxiliary vectors and tensors for anisotropy in `dielmat.f90`
  4. Standarize directory structure.
  
## 2018.11.4 (zmy)
  1. Change version name to `3.0a`

## 2018.11.5 (zmy)
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
  6. correctly calculate head at `q0_eps` with `ten_p_aniso` for `iop_aniso.ne.-1`
