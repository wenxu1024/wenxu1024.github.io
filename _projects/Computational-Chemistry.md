---
title: Some useful python and bash scripts for computational chemistry
layout: project
name: Some useful python and bash scripts for computational chemistry
---

Chemistry has evloved from "beaker and test tube" age for a long time. Paul Dirac once said  
>The underlying physical laws necessary for the mathematical theory of a large part of physics and the **whole of chemistry** are thus **completely known**, and the difficulty is only that the exact application of these laws leads to equations much too complicated to be soluble. It therefore becomes desirable that **approximate** practical methods of applying quantum mechanics should be developed, which can lead to an explanation of the main features of complex atomic systems without too much computation.
{:.text-justify}

Computational chemistry has been a very powerful and useful tool started from the last century due to the advances in computing power.

Below are three useful python or bash scripts written to calculate ion-molecule reaction rates, ion-molecule collision cross sections, and electron impact ionization cross sections based on information extracted from quantum chemistry calculations.
<!--more-->

## List of scripts

1. [Ion-molecule reaction rates calculator](#code-block-1)
2. [Ion-molecule collision cross section calculator](#code-block-2)
3. [Electron impactor cross section calculator](#code-block-3)

## Ion-molecule reaction rates calculator
{:#code-block-1}
```bash
#! /bin/bash
echo please give the reagent ion mass in amu
read amu
echo please give the temperature in Kelvin
read Temp
#echo $rmass,$T
echo species,reduced mass,alpha,ud,ud/sqrtalpha,C,k
for file in *.log; do
filename=${file%.log}
echo -n "$filename"
if grep -q "exited gracefully" $file; then
    true
else
    echo "	the Gamess job does not terminate normally, restart neccesarry"
    continue
fi
awk -v rmass="$amu" -v T="$Temp" '/ALPHA\(1,1\)/ {alpha11=$3}
				  /ALPHA\(2,2\)/ {alpha22=$3}
			          /ALPHA\(3,3\)/ {alpha33=$3}
		      	          /\/D\// {getline;dipole=$4}
				  /\$DATA/ {getline;getline;getline;
					    mass=0;
					    while($3 != "$END") {
					    	if($3 == 6) mass += 12.0;
					        if($3 == 1) mass += 1.0;
						if($3 == 8) mass += 16.0;
						if($3 == 7) mass += 14.0;
						if($3 == 16) mass += 32.0;
						if($3 == 17) mass += 35.0;
						if($3 == 15) mass += 31.0;
						if($3 == 11) mass += 23.0;
						if($3 == 19) mass += 39.0;
						getline;
	 				    }
					   }
END {
     u = rmass * mass / (rmass + mass)
     alpha=(alpha11+alpha22+alpha33)/3.0
     alpha=alpha*0.529117**3
     case=dipole/sqrt(alpha)
     if(case<0.082)
     C=0.0
     else
     { if( case<0.163)
       C=0.02
       else
       {if (case<0.250)
         C=0.045
        else
         {if (case<0.327)
          C=0.084
          else
           {if(case<0.408)
            C=0.115
            else
            { if(case<0.490)
              C=0.137
              else
              { if(case<0.572)
                C=0.156
                else
                 {if(case<0.653)
                  C=0.172
                  else
                   {if(case<0.735)
                    C=0.185
                    else
                     if(case<0.817)
                      C=0.196
                      else
                       {if(case<0.898)
                         C=0.205
                        else
                         { if(case<0.980)
                          C=0.213
                           else
                           {if(case<1.061)
                            C=0.220
                            else
                             {if(case<1.143)
                            C=0.226
                              else
                            C=0.231
                             }
                           }
                         }
                        }
                       }
                     }
              }
             }
           }
          }
         }
        }
q=1.602176487E-19*0.1*29979245800
pi=3.14159
k=2*pi*q/sqrt(u/6.0221415E23)*(sqrt(alpha*1E-24)+C*dipole*1E-18*sqrt(2/(pi*1.3806503E-16*T)))/(1E-9)
printf("%12.6f%12.6f%12.6f%12.6f%12.6f%12.2f\n",u,alpha,dipole,dipole/sqrt(alpha),C,k)}' $file
done
#k is in unit of 1E-9 cm3/(moleclue*sec)
```

## Bash script for generating \*.xyz file for mobcal to calculate ion-molecule collision cross section
{:#code-block-2}
```bash
#! /bin/bash
if grep -q "EQUILIBRIUM GEOMETRY LOCATED" $1; then
    true
else
    echo "the Gamess job type must be optimization"
    exit
fi
if grep -q "TOTAL MULLIKEN AND LOWDIN ATOMIC POPULATIONS" $1; then
    true
else
    echo "the Gamess job does not have Mulliken and lowdin population analysis"
    exit
fi
fullfilename=$1
filename=${fullfilename%.log}
filename=${filename##*/}
echo $filename
echo 1
awk ' 
    BEGIN{ }
    /EQUILIBRIUM GEOMETRY LOCATED/ {
	do {
	   getline var;
        } while(substr(var,17,3) != "ALL");
	getline;getline;
	i = 0;
	ntot = 0;
	do {
	   getline var;
	   if (var != "") {
	       S[i] = substr(var,0,2);
	       X[i] = substr(var,17,16);
	       Y[i] = substr(var,33,16);
	       Z[i] = substr(var,49,16);
	       i += 1;
	       ntot += 1;
	   }
	} while(var != "");
    }
    /TOTAL MULLIKEN AND LOWDIN ATOMIC POPULATIONS/ {
        getline;
        i = 0;
	do {
	   getline var;
	   if (var != "") {
		Q[i] = substr(var,58,13);
		i += 1;
	   }
	} while(var != "");
    }
    END {
	print ntot;
	print "ang";
	print "calc"
	print "1.0000"
	for(i = 0; i < ntot; i += 1) {
	    if(S[i] == " H") print X[i], Y[i], Z[i],"    ", 1.0, Q[i];
	    else if (S[i] == " C") print X[i], Y[i], Z[i],"   ", 12.0, Q[i];
	    else if (S[i] == " O") print X[i], Y[i], Z[i],"   ", 16.0, Q[i];
	    else if (S[i] == " N") print X[i], Y[i], Z[i],"   ", 14.0, Q[i];
	}
 }
' $1 
```

## Python script for calculating electron impact ionization cross sections
{:#code-block-3}
```python
#! /usr/bin/python
# a python script to parse the gamess output
# 1. Read 1 electron Kinetic Energy Integral
# 2. Read occupied MO coefficients and MO energies
# 3. Calculate the 1 electron kinetic energy of MO

import re
import numpy as np
# total number of cartesian gaussian basis functions pattern
pattern1 = re.compile("NUMBER OF CARTESIAN GAUSSIAN BASIS FUNCTIONS\s+=\s+(\d+)")
# total number of electrons pattern
pattern2 = re.compile("NUMBER OF ELECTRONS\s+=\s+(\d+)")
# kinetic energy integrals pattern
pattern3 = re.compile("\s+KINETIC ENERGY INTEGRALS\s+")
# MO coefficients and energies pattern
pattern4 = re.compile("\s+EIGENVECTORS\s+")
NFUNC = 0
NELEC = 0
#replace with your gamess result file name
fileName = "CS.log"
###################
filein = open(fileName, 'r+')
for i, line in enumerate(filein):
  for match in re.finditer(pattern1, line):
    if match:
      NFUNC = match.groups()[0]
    else:
      print "Gamess Execution Error"

filein = open(fileName, 'r+')
for i, line in enumerate(filein):
  for match in re.finditer(pattern2, line):
    if match:
      NELEC = match.groups()[0]
    else:
      print "Gamess Execution Error"

NFUNC = int(NFUNC)
NELEC = int(NELEC)

K_matrix = [[0 for j in range(NFUNC)] for i in range(NFUNC)]

filein = open(fileName, 'r+')
for i, line in enumerate(filein):
  for match in re.finditer(pattern3, line):
    if match == None:
      print "Gamess Execution Error, missing 1 electron kinetic energy integrals"
    else:
      filein.next()
      filein.next()
      filein.next()
      nblock = (NFUNC - 1) / 5 + 1 
      row = 0
      column = 0
      for iblock in range(nblock):
        row = 5 * iblock
        #print iblock, row, NFUNC
        for row_index in xrange(row, NFUNC):
          line = filein.next()
          if row_index == row:
            K_matrix[row_index][row_index] = float(line.split()[4])
          elif row_index == row + 1:
	    K_matrix[row_index][row] = float(line.split()[4])
	    K_matrix[row_index][row + 1] = float(line.split()[5])
          elif row_index == row + 2:
	    K_matrix[row_index][row] = float(line.split()[4])
            K_matrix[row_index][row + 1] = float(line.split()[5])
	    K_matrix[row_index][row + 2] = float(line.split()[6])
          elif row_index == row + 3:
	    K_matrix[row_index][row] = float(line.split()[4])
            K_matrix[row_index][row + 1] = float(line.split()[5])
            K_matrix[row_index][row + 2] = float(line.split()[6])
	    K_matrix[row_index][row + 3] = float(line.split()[7])
          else:
	    K_matrix[row_index][row] = float(line.split()[4])
            K_matrix[row_index][row + 1] = float(line.split()[5])
            K_matrix[row_index][row + 2] = float(line.split()[6])	
	    K_matrix[row_index][row + 3] = float(line.split()[7])
            K_matrix[row_index][row + 4] = float(line.split()[8])
        filein.next()
        filein.next()
        filein.next()

#print K_matrix
# Symmetrize K_matrix
m = len(K_matrix)
n = len(K_matrix[0])
for i in range(m):
  for j in xrange(i, n):
    K_matrix[i][j] = K_matrix[j][i]

K_matrix = np.matrix(K_matrix)
#print K_matrix
# allocate memergy for MO coefficients
C_matrix = [[0 for j in range(NFUNC)] for i in range(NFUNC)]
# allocate memergy for MO energies
E_Vector = [0 for j in range(NFUNC)]

filein = open(fileName, 'r+')
for i, line in enumerate(filein):
  for match in re.finditer(pattern4, line):
    if match == None:
      print "Gamess Execution Error missing MO coefficients and energies, possibly SCF does not Converge"
    else:
      filein.next()
      nblock = (NFUNC - 1) / 5 + 1
      col = 0
      for iblock in range(nblock):
        filein.next()
        filein.next()
        line = filein.next()
	col = 5 * iblock
	if col < NFUNC:
          E_Vector[col] = float(line.split()[0])
	if col + 1 < NFUNC:
	  E_Vector[col + 1] = float(line.split()[1])
	if col + 2 < NFUNC:
	  E_Vector[col + 2] = float(line.split()[2])
	if col + 3 < NFUNC:
	  E_Vector[col + 3] = float(line.split()[3])
	if col + 4 < NFUNC:
	  E_Vector[col + 4] = float(line.split()[3])
	filein.next()
	for row_index in range(NFUNC):
	  line = filein.next()
#	  print line.split()
          if col < NFUNC:
	    C_matrix[row_index][col] = float(line.split()[4])
	  if col + 1 < NFUNC:
            C_matrix[row_index][col + 1] = float(line.split()[5])
	  if col + 2 < NFUNC:
	    C_matrix[row_index][col + 2] = float(line.split()[6])
	  if col + 3 < NFUNC:
	    C_matrix[row_index][col + 3] = float(line.split()[7])
	  if col + 4 < NFUNC:
	    C_matrix[row_index][col + 4] = float(line.split()[8])

#print E_Vector
C_matrix = np.matrix(C_matrix)
#print C_matrix
C_Transpose = np.transpose(C_matrix)
K_temp = np.dot(np.dot(C_Transpose, K_matrix), C_matrix)

E_hartree = 27.211396132 # Hartree to eV
n = len(E_Vector)

#print K_matrix
#print C_matrix
#print E_Vector

Tk_array = [(i + 1) * 1.0 for i in range(1000)]           			#kinetic Energy of impact electron
pi = 3.14159265958
a0 = 0.5292 # Angstrom Bohr radius
R = 13.61 #eV half hartree

for Tk in Tk_array:
  sigma_tot = 0.0 # total electron impact cross section

  if NELEC % 2 == 0: # number of electrons is even. Most probably RHF.
    norb = NELEC / 2
    if norb > NFUNC:
      print "Error in Gamess SCF calculation"

    for i in range(norb):
      Neorb = 2 # electrons in single molecular orbital
      Uk = K_temp[i,i] * E_hartree
      Bind = -E_Vector[i] * E_hartree
      if Tk > Bind:
        t = Tk / Bind
        u = Uk / Bind
        S = 4.0 * pi * a0 * a0 * Neorb * R * R / Bind / Bind
        sigma_beb = S / (t + u + 1) * (np.log(t) / 2.0 * (1 - 1.0 / t / t) + 1 - 1.0 / t - np.log(t) / (t + 1))
      else:
        sigma_beb = 0
      sigma_tot += sigma_beb

  print Tk, sigma_tot
```