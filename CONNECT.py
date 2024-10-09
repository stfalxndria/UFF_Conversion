#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 14:46:36 2024

@author: stefani

mdf to get connectivity
"""
from input import *
from functions import *
atom_charge = None

start = datetime.now()
print(start)


car = f'{structure}/{structure}.car'
mdf = f'{structure}/{structure}.mdf'
pacmof = f'{structure}/{structure}_charged.cif'

#READING PACMOF FILE
with open(pacmof,'r') as f:
    lines = f.readlines()
pac = []
for line in lines:
    word = line.split()
    pac.append(word)

pac = [x for x in pac if x != []]
                

#READING CAR FILE
with open(car,'r') as f:
    lines = f.readlines()
car_ = []
for line in lines:
    word = line.split()
    car_.append(word)
    

    
car_table = pd.DataFrame() #LIST OF ATOMS AND LABELS
atom_list_unmod = []
atom_names = []
ff_names = []
atom_elements = []

target1 = '(P1)'
target2 = 'end'
inside_target = False
for i, line in enumerate(car_):
    if target1 in line:
        line = car_[i+1]
        inside_target = True
    elif target2 in line:
        inside_target = False

    if (inside_target == True) and (line not in atom_list_unmod):
        atom_list_unmod.append(line)
        atom_names.append(line[0])
        ff_names.append(line[6])
        atom_elements.append(line[7])
        
        
car_table['index'] = np.arange(1,len(atom_names)+1,1)
car_table['labels']= atom_names
car_table['elements']=atom_elements
car_table['uff names']=ff_names

#reading Charges file
charge_list = []
for i, row in enumerate(pac):
    if len(row)==9:
        atom_label = pac[i][0]
        charge = pac[i][8]
        for a, l1 in enumerate(car_table['labels']):
            if atom_label == l1:
                charge_list.append([atom_label,charge])

############################################################################
#READING THE MDP file
with open(mdf,'r') as f:
    lines = f.readlines()
mdf_unmod = []
mdf_ = []
for line in lines:
    word = line.split()
    mdf_unmod.append(word)
    

for i, line in enumerate(mdf_unmod):
    start_ = 21
    target = '!'
    if target in line:
        end_ = i-1
        mdf_ = mdf_unmod[start_:end_]
        
mdf_ = [[s.replace("XXXX_1:", "") for s in inner_list] for inner_list in mdf_]

##############################################################################
all_bond_list = pd.DataFrame()


atoms_t = []
bonded = []
for row in mdf_:
    x = 12
    y = len(row)+1
    atoms_t.append(row[0])
    bonded.append(row[x:y])
all_bond_list['atom']=atoms_t
all_bond_list['bonded'] = bonded

all_bond_list['bonded'] = [[atom.split('%')[0] for atom in row] for row in all_bond_list['bonded']]
all_bond_list['bonded'] = [[atom.split('/')[0] for atom in row] for row in all_bond_list['bonded']]

unique_atoms = []
for i, line in enumerate(car_table['uff names']):
    atom = car_table['uff names'][i]
    elements = car_table['elements'][i]
    if atom not in [item[0] for item in unique_atoms]:
        unique_atoms.append([atom, elements])

label_list = []
for line in unique_atoms:
    label_list.append(line[0])


#BOND CONNECTIVITY
bond_connectivity = []
unique_bonds = []
for i, line in enumerate(mdf_):
    atom_A = line[0]
    x = 12
    y = len(line)
    for b in range(x,y):
        if ('%' in line[b]) and ('/' in line[b]):
            temp = line[b].split('%')[0]
            atom_B = temp.split('/')[0]
        if '%' in line[b]:
            atom_B = line[b].split('%')[0]
        elif '/' in line[b]:
            atom_B = line[b].split('/')[0]
        else:
            atom_B = line[b]
        if ([atom_A,atom_B] not in bond_connectivity) and ([atom_B,atom_A] not in bond_connectivity):
            bond_connectivity.append([atom_A,atom_B])
            
print('reading done')

#GETTING THE BONDS FOR RASPA DEF file
for b, ll in enumerate(bond_connectivity):
    uff_A = None
    uff_B = None
    atom_A = ll[0]
    atom_B = ll[1]
    for a, line in enumerate(car_table['labels']):
        if atom_A == car_table['labels'][a]:
            uff_A = car_table['uff names'][a]
    for c, l3 in enumerate(car_table['labels']):
        if atom_B == car_table['labels'][c]:
            uff_B = car_table['uff names'][c]
    
    if ([uff_A,uff_B] not in unique_bonds) and ([uff_B, uff_A] not in unique_bonds):
        unique_bonds.append([uff_A,uff_B])

print('bonds done')
            
#######################################################################
#angle connectivity
angle_connectivity = []
for i, line in enumerate(bond_connectivity):
    atom_A, atom_B = line[0], line[1]

    for a, line in enumerate(car_table['labels']):
        atom_C = car_table['labels'][a]
        
        angle1 = None
        angle2 = None
        
        if atom_B != atom_C: #checking for bond A-C, C must not be B
            if ([atom_A, atom_C] in bond_connectivity) or ([atom_C,atom_A] in bond_connectivity):
                angle1 = [atom_C, atom_A, atom_B]

        if atom_A != atom_C:#checking for bond B-C, C must not be A
            if ([atom_B, atom_C] in bond_connectivity) or ([atom_C, atom_B] in bond_connectivity):
                angle2 = [atom_A, atom_B, atom_C]

        if ([atom_C, atom_A, atom_B] not in angle_connectivity) and ([atom_B, atom_A, atom_C] not in angle_connectivity) and (angle1 != None):
            angle_connectivity.append(angle1)

        if ([atom_A, atom_B, atom_C] not in angle_connectivity) and ([atom_C, atom_B, atom_A] not in angle_connectivity) and (angle2 != None):
            angle_connectivity.append(angle2)
            
unique_angles = []

for i, line in enumerate(angle_connectivity):
    uff_A = None
    uff_B = None
    uff_C = None
    atom_A, atom_B, atom_C = line[0], line[1], line[2]
    for a, l1 in enumerate(car_table['labels']):
        if atom_A == car_table['labels'][a]:
            uff_A = car_table['uff names'][a]
    for b, l2 in enumerate(car_table['labels']):
        if atom_B == car_table['labels'][b]:
            uff_B = car_table['uff names'][b]
    for c, l3 in enumerate(car_table['labels']):
        if atom_C == car_table['labels'][c]:
            uff_C = car_table['uff names'][c]
    
    if ([uff_A,uff_B,uff_C] not in unique_angles) and ([uff_C,uff_B, uff_A] not in unique_angles):
        unique_angles.append([uff_A,uff_B,uff_C])
        
print('angles done')
        
##############################################################################
#Torsion Connectivity

torsion_connectivity=[]
for i, line in enumerate(angle_connectivity):
    atom_A, atom_B, atom_C = line[0], line[1], line[2]
    for a, line in enumerate(car_table['labels']):
        atom_D = car_table['labels'][a]
        
        tors1 = None
        tors2 = None
        
        if atom_B != atom_D:
            if ([atom_A, atom_D] in bond_connectivity) or ([atom_D, atom_A] in bond_connectivity):
                tors1 = [atom_D, atom_A, atom_B, atom_C]
        if atom_B != atom_D:
            if ([atom_C, atom_D] in bond_connectivity) or ([atom_D, atom_C] in bond_connectivity):
                tors2 = [atom_A, atom_B, atom_C, atom_D]
                
        if ([atom_D, atom_A, atom_B, atom_C] not in torsion_connectivity) and ([atom_C, atom_B, atom_A, atom_D] not in torsion_connectivity) and (tors1 != None):
            torsion_connectivity.append(tors1)
            
        if ([atom_A, atom_B, atom_C, atom_D] not in torsion_connectivity) and ([atom_D, atom_C, atom_B, atom_A] not in torsion_connectivity) and (tors2 != None):
            torsion_connectivity.append(tors2)
        
unique_torsions = []

for i, line in enumerate(torsion_connectivity):
    uff_A = None
    uff_B = None
    uff_C = None
    uff_D = None
    atom_A, atom_B, atom_C, atom_D = line[0], line[1], line[2], line[3]
    for a, l1 in enumerate(car_table['labels']):
        if atom_A == car_table['labels'][a]:
            uff_A = car_table['uff names'][a]
    for b, l2 in enumerate(car_table['labels']):
        if atom_B == car_table['labels'][b]:
            uff_B = car_table['uff names'][b]
    for c, l3 in enumerate(car_table['labels']):
        if atom_C == car_table['labels'][c]:
            uff_C = car_table['uff names'][c]
    for d, l3 in enumerate(car_table['labels']):
        if atom_D == car_table['labels'][d]:
            uff_D = car_table['uff names'][d]
    
    if ([uff_A,uff_B,uff_C, uff_D] not in unique_torsions) and ([uff_D, uff_C,uff_B, uff_A] not in unique_torsions):
        unique_torsions.append([uff_A,uff_B,uff_C,uff_D])     

print('torsion_done')

##############################################################################
#Inversion Connectivity
inversion_connectivity = []

for i, atom in enumerate(all_bond_list['atom']):
    atom_j = atom
    bonded_ = all_bond_list['bonded'][i]
    inversion_ = None
    if len(bonded_) >= 3:
        y = list(combi(bonded_, 3))
        inversion_ = [(atom_j,) + combo for combo in y]
        
        for x in inversion_:
            if x not in inversion_connectivity:
                inversion_connectivity.append(x)
                
unique_inversions = []

for i, line in enumerate(inversion_connectivity):
    uff_A = None
    uff_B = None
    uff_C = None
    uff_D = None
    atom_A, atom_B, atom_C, atom_D = line[0], line[1], line[2], line[3]
    for a, l1 in enumerate(car_table['labels']):
        if atom_A == car_table['labels'][a]:
            uff_A = car_table['uff names'][a]
    for b, l2 in enumerate(car_table['labels']):
        if atom_B == car_table['labels'][b]:
            uff_B = car_table['uff names'][b]
    for c, l3 in enumerate(car_table['labels']):
        if atom_C == car_table['labels'][c]:
            uff_C = car_table['uff names'][c]
    for d, l3 in enumerate(car_table['labels']):
        if atom_D == car_table['labels'][d]:
            uff_D = car_table['uff names'][d]
    
    if ([uff_A,uff_B,uff_C, uff_D] not in unique_inversions) and ([uff_D, uff_C,uff_B, uff_A] not in unique_inversions):
        unique_inversions.append([uff_A,uff_B,uff_C,uff_D])     

    
print('inversion_done')

    
print('file', structure,'was read')
print("--- %s time taken ---" % ((datetime.now() - start)))





    
    

        
        
