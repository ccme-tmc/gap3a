#!/usr/bin/env python
import sys,os,shutil,commands
from math import *

# python functions used by gap${version}_* scripts

Ry2eV = 13.6056917


s_real  = [ ['s', 0,[1.0,0.0] ] ]

p_real  = [ ['px',1,[0.70710677,0.0,-0.70710677,0.0] ],\
            ['py',1,[0.0,0.70710677, 0.0,0.70710677] ],\
            ['pz',0,[1.0,0.0] ] ]

d_real  = [ ['dz2',    0,[1.0,0.0] ],\
            ['dx2-y2', 2,[0.70710677, 0.0, 0.70710677,0.0] ],\
            ['dxy',    2,[0.0,0.70710677, 0.0,-0.70710677] ],\
            ['dyz',    1,[0.0,0.70710677, 0.0,0.70710677 ] ],\
            ['dzx',    1,[0.70710677,0.0, -0.70710677,0.0] ]]

f_real  = [ ['fxy2',      1,[ 0.70710678, 0.00000000,-0.70710678, 0.00000000] ],\
            ['fyz2',      1,[ 0.00000000, 0.70710678, 0.00000000, 0.70710678] ],\
            ['fz3',       0,[ 1.00000000, 0.00000000] ],\
            ['fx(x2-3y2)',3,[ 0.70710678, 0.00000000,-0.70710678, 0.00000000] ],\
            ['fy(3x2-y2)',3,[ 0.00000000, 0.70710678, 0.00000000, 0.70710678] ],\
            ['fz(x2-y2)', 2,[ 0.70710678, 0.00000000, 0.70710678, 0.00000000] ],\
            ['fxyz',      2,[ 0.00000000, 0.70710678, 0.00000000,-0.70710678] ] ]

def Coef_hyb(symbl):
  """
  Return the coefficients in terms of spherical harmonics for common hybrid functionals 
  in ther order of 
    [ coef_s, coef_p{-1/0/1}, coef_d_{-2,-1,0,1,2}] 
  """
  coef_all = []

  # the coefficients of real atomic orbitals in terms of spherical harmonics 
  # in the order of -l, -l+1, ..., l-1, l
  # based on https://en.wikipedia.org/wiki/Table_of_spherical_harmonics
  s2r = sqrt(2.0)/2
  s2i = complex(0.0,sqrt(2.0)/2)

  px=[ s2r, 0.0,-s2r]   
  py=[ s2i, 0.0, s2i]
  pz=[ 0.0, 1.0, 0.0]

  dz2   = [ 0.0, 0.0,  1., 0.0, 0.0]
  dx2y2 = [ s2r, 0.0,  0., 0.0, s2r]
  dxy   = [ s2i, 0.0,  0., 0.0,-s2i]
  dyz   = [ 0.0, s2i,  0., s2i, 0.0]
  dxz   = [ 0.0, s2r,  0.,-s2r, 0.0]

  zero=[]
  for i in range(9): zero.append( complex(0.0, 0.0) ) 

  if symbl == 'sp3':
    for iorb in range(4): 
      coef=zero[:]
      coef[0]= complex(0.5,0.0) 
      for i in range(3):
        if iorb == 0:
          coef[1+i] += 0.5*( px[i]+py[i]+pz[i])
        elif iorb == 1:
          coef[1+i] += 0.5*( px[i]-py[i]-pz[i])
        elif iorb == 2:
          coef[1+i] += 0.5*(-px[i]+py[i]-pz[i])
        elif iorb == 3:
          coef[1+i] += 0.5*(-px[i]-py[i]-pz[i])
      coef_all.append(coef[:]) 

  elif symbl == 'sp2':
    for iorb in range(3):
      coef=zero[:]
      coef[0] = complex(1.0/sqrt(3.0),0.0) 
      for i in range(3):
        if   iorb == 0:
          coef[1+i] += - px[i]/sqrt(6.0) + py[i]/sqrt(2.0)
        elif iorb == 1:
          coef[1+i] += - px[i]/sqrt(6.0) - py[i]/sqrt(2.0)
        elif iorb == 2:
          coef[1+i] +=   px[i]*2/sqrt(6.0)
      coef_all.append(coef[:]) 
  else:
    print "ERROR: not implemented for symbl=",symbl 
    sys.exit(1)

  return coef_all

def Write_orb_inwf(ofl,l,iat,atm,coefs=None):
  """
  Write coefficients of projected orbitals in terms of complex harmonics 
  """

  fmt    = "%2d%2d%3d%13.8f%13.8f%47s\n"
  info   ="# index of atom, L, M, coefficient (complex)"

  if l < 0: # composite projects (following the convention in Wannier90, Table 3.2 in wannier90 guide)  
    if l == -3:    # sp3 hybrid 
      symbl = 'sp3'
      norb = 4
    elif l == -2:  # sp2 hybrid 
      symbl = 'sp2'
      norb = 3 
    elif l == -1:
      symbl = -1
      norb = 2 
    elif l == -4: 
      symbl='sp3d'
      norb = 5
    elif l == -5:
      symbl='sp3d2'
      norb = 6 
    else:
      print "Error: ill-defined value for negative l=",l
      sys.exit(1) 
    coefs = Coef_hyb(symbl) 

    for iorb in range(norb): 
      coef=coefs[iorb]

      nc = 0   # count the number of compolents 
      for j in range(9): 
        if abs(coef[j]) > 0.000001: nc += 1

      ofl.write("%d        # %2s %s orbital\n"%(nc,atm,symbl) )
      
      for j in range(9): 
        if abs(coef[j]) > 0.000001:
          if j == 0:
            l=0; m=0
          elif j>= 1 and j <=3:
            l=1; m = j-2
          else:
            l=2; m = j-6
          ofl.write(fmt%(iat,l,m,coef[j].real,coef[j].imag,info))
    return

  norb = len(coefs)  
  for iorb in range(norb):
    symbl = coefs[iorb][0]
    m =     coefs[iorb][1]
    if m == 0:
      ofl.write("1         # %2s %s orbital\n"%(atm,symbl) )
      ofl.write(fmt%(iat,l,0,coefs[iorb][2][0],coefs[iorb][2][1],info))
    else:
      ofl.write("2         # %2s %s orbital\n"%(atm,symbl) )
      ofl.write(fmt%(iat,l,-m,coefs[iorb][2][0],coefs[iorb][2][1],info))
      ofl.write(fmt%(iat,l, m,coefs[iorb][2][2],coefs[iorb][2][3],info))

def Get_result(gwout,mode=0):
   """
   Use to extract the most important results from the GW output
    mode = 0 -> Eg_KS Eg_G0W0  DeltaVBM_G0W0 Eg_GW0 DeltaVBM_GW0 
   """
   Eg_PBE = f_Grep_Lines(gwout,nl,nf) 
  

def f_Getopt(opflag,nval,def_val,debug=False):
  try:
    i_op = sys.argv.index(opflag)
    del sys.argv[i_op]

    if nval == 0:
      val = True

    elif nval ==1:
      if isinstance(def_val,int):
        val = int(sys.argv[i_op]); del sys.argv[i_op]
      elif isinstance(def_val,float):
        val = float(sys.argv[i_op]); del sys.argv[i_op]
      else:
        val = sys.argv[i_op].strip(); del sys.argv[i_op]
    else:

      val=[]
      for i in range(nval):
        if isinstance(def_val[i],int):
          t = int(sys.argv[i_op]); del sys.argv[i_op]
        elif isinstance(def_val[i],float):
          t = float(sys.argv[i_op]); del sys.argv[i_op]
        else:
          t = sys.argv[i_op].strip(); del sys.argv[i_op]
        val.append(t)

  except:
    val = def_val

  if opflag != '-h' and debug : print "option " + opflag + " = ",val
  return val

def f_Check_Name(name=''):
  """
  Check whether <name> is a valid case name for the wien2k struct, i.e. check <name>.struct exists.
  If it is a null string, get the case name from the name of the current directory.
  the valid case name is returned, or if exit with an error message
  """
  if name == '' :
    cwd=os.getcwd()
    case_name = os.path.basename(cwd)
  else:
    case_name = name

  struct_file = case_name+".struct"

  if not os.path.isfile(struct_file) :
    print "ERROR: struct file " + struct_file + "not existing "
    sys.exit(1)
  return case_name

def f_Check_Complex(name):
  """
  Check whether Kohn-Sham vectores are complex, i.e. whether the system has the inversion symmetry
  """
  if os.path.isfile(name+".in1c") and not os.path.isfile(name+".in1"):
    cmplx='c'
  else:
    cmplx=''
  return cmplx

def f_Grep_Lines(file,tag,nl,nf,debug=False):
  """
  Get the "nf"-th field (separated by space) on the "nl"-th line
  that contains "tag"
  """

  if nl > 0:
    grep_cmd = "grep --no-filename '" + tag + "' " + file + "| head -n " + "%4d"%(nl)
  else:
    grep_cmd = "grep --no-filename '" + tag + "' " + file + "| tail -n " + "%4d"%(-nl)

  if debug: print "grep_cmd='",grep_cmd,"'"

  failure,output = commands.getstatusoutput( grep_cmd )
  if failure:
    print "ERROR when running " + grep_cmd
    sys.exit(1)

  if output.strip() == '':
    return None
  else:
    return output.split()[nf-1]

def f_Read_Efer(name,iop=0,sp=0,tag=''):
  """
  Read Fermi energy from case.scf/qtl(iop=0/1) 
  """
  case_name = f_Check_Name(name)

  if sp == 0:
    sptag=''
  else:
    sptag = 'up'

  if iop == 1:
    file=case_name.strip()+".scf2"+sptag+tag
    efer = float( f_Grep_Lines(file,":FER  :",1,10) )
  else:
    file=case_name.strip()+".scf"+tag
    efer = float( f_Grep_Lines(file,":FER  :",-1,10) )

  return efer


def f_Get_Natom(case_name):
  """
  Get the information about the number of atoms from the wien2k structure file,
  """
  ifile = open(case_name+".struct",'r')
  line = ifile.readline()
  line = ifile.readline()
  nat_neq = int( line[27:30] )
  ifile.close()
  return nat_neq

def f_Skip_Lines(ifile,n):
  i=0
  while i < n :
    ifile.readline()
    i += 1
  return

def f_Get_Atoms(case_name):
  """
  Get the elemental symbols of all atoms 
  """
  ifile = open(case_name+".struct",'r')

  f_Skip_Lines(ifile,1)
  line = ifile.readline()
  nat_neq = int( line[27:30] )
  f_Skip_Lines(ifile,2)

  atoms = []
  for i_neq in range(nat_neq):
    f_Skip_Lines(ifile,1)
    line = ifile.readline()
    m = int(line[15:17])
    f_Skip_Lines(ifile,m-1)
    line = ifile.readline().split()
    atoms.append(line[0])
    f_Skip_Lines(ifile,3)

  ifile.close()
  return atoms 

def f_Get_iat(atoms,at):
  """
  Get the atom index  
  """
  iat = -1
  nat = len(atoms)
  for i in range(nat):
    if len(atoms[i]) == 3 and len(at) == 2: 
      if atoms[i][0:2] == at: 
        iat = i 
        break 
    else:
      if atoms[i] == at: 
        iat = i
        break
  iat += 1 
  return iat 

def f_Get_iats(atoms,at):
  """
  Get the indices for the atom "at"  
  """
  iat_s = []
  nat = len(atoms)
  for i in range(nat):
    # removing the possible numbers in atoms symbol 
    if atoms[i][-1:].isdigit():
      atom = atoms[i][:-1]
    else:
      atom = atoms[i]  

    if atom == at: 
      iat_s.append(i) 
  return iat_s

def f_Get_Mult(case_name):
  """
  Get the information about the multiplicity of each atom
  """
  ifile = open(case_name+".struct",'r')

  f_Skip_Lines(ifile,1)
  line = ifile.readline()
  nat_neq = int( line[27:30] )
  f_Skip_Lines(ifile,2)

  nat = 0
  mult=[]
  for i_neq in range(nat_neq):
    f_Skip_Lines(ifile,1)
    line = ifile.readline()
    m = int(line[15:17])
    nat += m
    f_Skip_Lines(ifile,m+3)
    mult.append(m)

  ifile.close()
  return mult


def f_Read_Energy(fname,nat,debug=True):
  """
  Read band energies in the wien2k case.energy format 
  """
  print "Read band energies from " + fname
  ifile = open(fname,'r')

  # skip the heads related to the linearization energy 
  f_Skip_Lines(ifile,nat*2)
  ik=0
  enk_all = []
  kw_all = []
  kv_all = []
  while 1:
    line = ifile.readline()
    if not line: break # end of file readed 

    ik += 1
    # k-vectors, number of bands at each k, and the weight of this k     
    kv = [ float(line[0:19]), float(line[19:38]), float(line[38:57]) ]
    nbk = int(line[73:79])
    kw =  float(line[79:84])

    if debug:
      print "Read the band energy for ik=%5d,k=(%8.4f, %8.4f, %8.4f), nbk=%5d, wk=%5.1f"%(ik,kv[0],kv[1],kv[2],nbk,kw)
      print "%5s %12s"%("n","Enk")

    enk = []
    for ie in range(nbk):
      s_en = ifile.readline().split()
      n = int(s_en[0])
      en = float(s_en[1])
      enk.append(en)
      if debug: print "%5d %12.6f"%(n,en)

    kv_all.append(kv)
    kw_all.append(kw)
    enk_all.append(enk)

  ifile.close()

  return enk_all,kw_all,kv_all

def f_Band_Analysis(enk,efer,kvec,emin=-2.0,emax=2.0,ib0=1,debug=False):
  print \
  """
  Analyze band energies 
  """
  nsp = len(enk)
  nk  = len(kvec)
  nband = len(enk[0][0])
  print "\tNumber of spin=",nsp
  print "\tNumber of k-vector=",nk
  print "\tNumber of bands=",nband
  print "\n"

  nomax = []
  numin = []
  ikvm =  []
  ikcm =  []

  lmetal = False
  for isp in range(nsp):
    nomax.append(0)
    numin.append(nband)
    ikvm.append(0)
    ikcm.append(0)
    for ik in range(nk):
      nbk = len(enk[isp][ik])

      if debug:
        print "  %5d bands for ik=%5d"%(nbk,ik)

      if nbk < nband: nband = nbk
      nv = 0
      nc = nbk-1
      for ib in range(nbk):
        if enk[isp][ik][ib] <= efer + 0.00001:
          if ib > nv: nv = ib
        else:
          if ib < nc: nc = ib
      if nv > nomax[isp]:  nomax[isp] = nv
      if nc < numin[isp]:  numin[isp] = nc

      nv  = nomax[isp]
      nc  = numin[isp]
      ikv = ikvm[isp]
      ikc = ikcm[isp]
      if enk[isp][ik][nv] > enk[isp][ikv][nv]: ikvm[isp] = ik
      if enk[isp][ik][nc] < enk[isp][ikc][nc]: ikcm[isp] = ik
    # loop over ik

    if nomax[isp] >= numin[isp]: lmetal = True

  # end loop over isp

  # set evbm, which is used as the energy zero 
  if lmetal :
    evbm = efer
  else:
    if nsp==1:
      evbm = enk[0][ikvm[0]][nomax[0]]
    else:
      evbm = max( enk[0][ikvm[0]][nomax[0]], enk[1][ikvm[1]][nomax[1]] )

  for isp in range(nsp):
    nv = nomax[isp]
    nc = numin[isp]
    ikv = ikvm[isp]
    ikc = ikcm[isp]
    if nsp ==2: print "\nBand analysis for spin",isp
    print "Band index for VBM and CBM=",nv+ib0,nc+ib0

    if not lmetal:
      print "Insulating system:"
      egap = (enk[isp][ikc][nc] - enk[isp][ikv][nv])*Ry2eV
      if ikvm[isp] == ikcm[isp]:
        print ":BandGap(d) =%8.3f eV "%(egap)
        print "  direct gap at k=(%8.3f,%8.3f,%8.3f)"%(kvec[ikv][0],kvec[ikv][1],kvec[ikv][2])
      else:
        print ":BandGap(i) =%8.3f eV "%(egap)
        egvbm = (enk[isp][ikv][nc] - enk[isp][ikv][nv] )*Ry2eV
        egcbm = (enk[isp][ikc][nc] - enk[isp][ikc][nv] )*Ry2eV

        egm = (enk[isp][0][nc] - enk[isp][0][nv])*Ry2eV
        ikm = 0
        for ik in range(nk):
          egk = (enk[isp][ik][nc] - enk[isp][ik][nv])*Ry2eV
          if egk < egm:
            egm = egk
            ikm = ik

        print ":Eg_d(min) =%8.3f eV, at      k=(%8.3f,%8.3f,%8.3f) (ik=%4d)"%(egm,  kvec[ikm][0], kvec[ikm][1], kvec[ikm][2],ikm+1)
        print ":Eg_direct =%8.3f eV  at  VBM k=(%8.3f,%8.3f,%8.3f) (ik=%4d)"%(egvbm,kvec[ikv][0], kvec[ikv][1], kvec[ikv][2],ikv+1)
        print ":Eg_direct =%8.3f eV  at  CBM k=(%8.3f,%8.3f,%8.3f) (ik=%4d)"%(egcbm,kvec[ikc][0], kvec[ikc][1], kvec[ikc][2],ikc+1)
    else:
      print "Metallic system"
    print "Range of each band with respect to VBM (eV):"
    print "%5s%12s%12s%12s"%('n','Bottom','Top','Width')
    for ib in range(nband):

      ebmin=enk[isp][0][ib]
      ebmax=enk[isp][0][ib]

      for ik in range(nk):
        if enk[isp][ik][ib] < ebmin:
          ebmin = enk[isp][ik][ib]
        if enk[isp][ik][ib] > ebmax:
          ebmax = enk[isp][ik][ib]

      ebmin = ebmin - evbm
      ebmax = ebmax - evbm

      if ebmin > emin  and ebmax < emax:
        print "%5d%12.3f%12.3f%12.3f"%(ib+ib0,ebmin*Ry2eV,ebmax*Ry2eV,(ebmax-ebmin)*Ry2eV)
      # end loop over ib
   # end loop over isp 

  return

