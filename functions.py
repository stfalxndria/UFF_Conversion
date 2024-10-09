#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 14 17:52:45 2024

@author: stefani
"""

#importing packages needed
import pandas as pd
import numpy as np
import math
import os
from datetime import datetime
from itertools import combinations as combi

from molmass import Formula


#Reading the UFF table data
_uff='data/uff_table'
with open(_uff,'r') as f:
    lines = f.readlines()
uff = []
for line in lines:
    words = [x.strip() for x in line.split(',')]
    uff.append(words)

uff_table = pd.DataFrame(uff, columns=['atoms', 'ri', 'phi', 'xi', 'di', 'psi', 'zmm', 'vsp3', 'vsp2', 'chi', 'nc'])

_angleff='data/angle.ff'
with open(_angleff,'r') as f:
    lines = f.readlines()
angle_ff = []
for line in lines:
    words = [x.strip() for x in line.split('\t')]
    angle_ff.append(words)
    
_torsionff='data/torsion.ff'
with open(_torsionff,'r') as f:
    lines = f.readlines()
torsion_ff = []
for line in lines:
    words = [x.strip() for x in line.split('\t')]
    torsion_ff.append(words)
    
_inversionff='data/inversion.ff'
with open(_inversionff,'r') as f:
    lines = f.readlines()
inversion_ff = []
for line in lines:
    words = [x.strip() for x in line.split('\t')]
    inversion_ff.append(words)

def get_hybridisation(atom):
    _hyb = None
    
    if len(atom) >=3:
        if atom[2]==1:
            _hyb = 'sp'
        elif atom[2]==2:
            _hyb = 'sp2'
        elif atom[2]==2:
            _hyb = 'aromatic'
        else:
            _hyb = 'sp3'
    elif len(atom) < 3:
        _hyb = 'sp3'
    
    return _hyb


def get_BO(atom_A,atom_B):
    _bo = None
    A_hyb=get_hybridisation(atom_A)
    B_hyb=get_hybridisation(atom_B)
                    
    if (atom_A == 'C_R' and atom_B == 'N_R') or (atom_A == 'N_R' and atom_B == 'C_R'):
        _bo = 1.41
    elif (atom_A == 'O_3z' and atom_B == 'Si3') or (atom_A == 'O_3z' and atom_B == 'Si3'):
        _bo = 1.44
    elif A_hyb == 'sp' and B_hyb == 'sp':
        _bo = 3.0
    elif A_hyb == 'sp2' and B_hyb == 'sp2':
        _bo = 2.0
    elif A_hyb == 'aromatic' and B_hyb == 'aromatic':
        _bo = 1.5
    else:
        _bo = 1.0
    return _bo

def get_natural_bond(atom_A,atom_B):
    r_ab = None
    r_i = None
    x_i = None
    r_j = None
    x_j = None
    n_bo = get_BO(atom_A, atom_B)

    for n, line in enumerate(uff_table['atoms']):
        if atom_A == uff_table['atoms'][n]:
            r_i = float(uff_table['ri'][n])
            x_i = float(uff_table['xi'][n])

    for a, ll in enumerate(uff_table['atoms']):
        if atom_B == uff_table['atoms'][a]:
            r_j = float(uff_table['ri'][a])
            x_j = float(uff_table['xi'][a])
            
            #bond order correction term
            r_bo_ij =-0.1332*(r_i+r_j)*math.log(n_bo)
            
            #equilibrium bond distance correction term
            ren_ij = r_i * r_j * ((math.sqrt(x_i) - math.sqrt(x_j))**2) / (x_i * r_i + x_j * r_j)
            
            #final equation: Natural Bond Distance
            r_ab = r_i+r_j-ren_ij+r_bo_ij
            
    return r_ab


#grouping systems needed for the separation of atoms
#checking whether the atom is a group_6 or not
g6 = ['O','S','Se','Te','Po']