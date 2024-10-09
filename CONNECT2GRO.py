#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 14:52:17 2024

@author: stefani
UFF For GROMACS
"""
#Variables needed to be set to generate an accurate file

from input import *
from CONNECT import *

print(f"there are {len(unique_atoms)} unique atoms, {len(unique_bonds)} unique bonds, {len(unique_angles)} unique angles and {len(unique_torsions)} unique torsions observed")
start = datetime.now()
    

###############################################################################
#.TOP file

top_list = []
top_list.append(['#include "force_field.itp"'])
top_list.append([f'#include "{groupid}.itp" '])
top_list.append([''])
top_list.append(['[ system ]'])
top_list.append([f'{groupid}','',''])
top_list.append([''])
top_list.append(['[ molecules ]'])
top_list.append([groupid,1])


top_name=f'{structure}/topol.top'
with open(top_name, 'w') as file:
    for row in top_list:
        row_str = '\t'.join(map(str, row)) 
        file.write(row_str + '\n')
        
print("topol.top file printed")

###############################################################################
print(start)
#masses
masses = []
for i in range(len(unique_atoms)):
    atom_uff = unique_atoms[i][0]
    for x, b in enumerate(uff_table['atoms']):
        if atom_uff == b:
            atom_mass = uff_table['mass'][x]
            masses.append([atom_uff,atom_mass])
    
#Lennard Jones
atom_types = []
#Lennard Jones Potential in GROMACS UNITS
for n in masses:
    _atom=n[0]
    _mass=round(float(n[1]),3)
    for i, line in enumerate(uff_table['atoms']):
        if _atom == uff_table['atoms'][i]:
            _energy=round(float(uff_table['di'][i])*4.184,4)
            _radius=round((float(uff_table['xi'][i]))*0.89089871814*0.1,6)
            atom_types.append([_atom,'1',_mass,'0.000','A',_radius,_energy])
            
            
#BOND Order calculation
#assigning hybridisation to each atoms
atom_hybridisation=pd.DataFrame()
hyb_list = []
for i in range(len(label_list)):
    hyb_list.append(get_hybridisation(label_list[i]))

atom_hybridisation['atoms']=label_list
atom_hybridisation['hybridisation']=hyb_list

#now assigning bond orders for each pairs
BO_list=[]
for i in range(len(unique_bonds)):
    atom_A = unique_bonds[i][0]
    atom_B = unique_bonds[i][1]
    BO_list.append([atom_A,atom_B,get_BO(atom_A, atom_B)])

bond_types=[]

#DISTANCE: calculating the r_ij value 
    #r_ij = r_i + r_j + r_BO - r_EN
r_ij=[] 
r_ij_G = []

for i in range(len(unique_bonds)):
    atom_A = unique_bonds[i][0]
    atom_B = unique_bonds[i][1]
    _r = get_natural_bond(atom_A, atom_B)
    
    r_ij.append([atom_A, atom_B, round(_r,4)])
    r_ij_G.append([atom_A, atom_B, round(_r*0.1,4)])
                
    
#ENERGY: calculating the k_ij value
k_ij=[] #kcal/mol/A^2

for i in range(len(unique_bonds)):
    atom_A = unique_bonds[i][0]
    atom_B = unique_bonds[i][1]
    z_i = None
    z_j = None
    rab=[]

    for n, line in enumerate(uff_table['atoms']):
        if atom_A == uff_table['atoms'][n]:
            z_i = float(uff_table['zmm'][n])

    for a, ll in enumerate(uff_table['atoms']):
        if atom_B == uff_table['atoms'][a]:
            z_j = float(uff_table['zmm'][a]) 
    for b, lll in enumerate(r_ij):
        if atom_A in lll[0] and atom_B in lll[1]:
            rab=lll[2]
            
            k_ab=664.12*((z_i*z_j)/(rab))
            k_ij.append([atom_A, atom_B,k_ab])

#converting the UFF units to the RASPA units
k_ij_G = []
for i in range(len(k_ij)):
    kab=k_ij[i][2]*418.4    #kJ/mol/nm^2
    kab=round(kab,2)
    k_ij_G.append(kab)
    
for i, line in enumerate(k_ij):
    atom_A = line[0]
    atom_B = line[1]
    _energy = k_ij_G[i]
    _radius = r_ij_G[i][2]
    
    bond_types.append([atom_A, atom_B, '1', _radius, _energy])

#ANGLES

#obtaining the middle atom natural angle
theta=[] #in degrees
for i in range(len(unique_angles)):
    _atom=unique_angles[i][1]
    for n, row in enumerate(uff_table['atoms']):
        if _atom == uff_table['atoms'][n]:
            theta.append([_atom,uff_table['phi'][n]])

#calculate the force constant K_ijk 
K_ijk=[] #kcal/mol/rad*2

for i, line in enumerate(unique_angles):
    atom_A=line[0]  #i
    atom_B=line[1]  #j
    atom_C=line[2]  #k
    
    z_i = None
    z_k = None
    angle = None
    
    #calculating the natural bonds
    _r_ij=get_natural_bond(atom_A, atom_B)  #atom 1 and 2
    _r_jk=get_natural_bond(atom_B, atom_C)  #atom 2 and 3

    
    for n, row in enumerate(theta):
        if atom_B == theta[n][0]:
            angle=float(theta[n][1])
            theta_rad = math.radians(angle)
            
    for b, ll in enumerate(uff_table['atoms']): 
        if atom_A == uff_table['atoms'][b]:
            z_i = float(uff_table['zmm'][b])
            
    for d, llll in enumerate(uff_table['atoms']):
        if atom_C == uff_table['atoms'][d]:
            z_k = float(uff_table['zmm'][d])
            
            _r_ik=math.sqrt((_r_ij**2)*(_r_jk**2)-2*(_r_ij*_r_jk*math.cos(theta_rad)))  #atom 1 and 3
            
            
            #K_ijk equation
            _beta = 664.12/(_r_ij*_r_jk)
            
            _back = 3*_r_ij*_r_jk*(1-((math.cos(theta_rad))**2))-(((_r_ik)**2)*math.cos(theta_rad))
            
            _kijk= _beta*((z_i*z_k)/(_r_ik)**5)*(_r_ij*_r_jk)*_back
            
            K_ijk.append([atom_A, atom_B, atom_C,_kijk])
            
            
#final equation
angle_types = []

for a, line in enumerate(unique_angles):
    atom_A=line[0]
    atom_B=line[1]
    atom_C=line[2]

    #the terms needed to calculate the coefficients
    K = None
    angle = None #in degrees
    theta_rad = None #in radius, for calculation purposes
    deg_harm = None
    k_harm = None
    
    for n, row in enumerate(theta):
        if atom_B == row[0]:
            angle=float(row[1])
            theta_rad = math.radians(angle)
            
    for i, sub in enumerate(K_ijk):
        if atom_B == K_ijk[i][1]:
            K = round(4.184*K_ijk[i][3],2) #UFF (Kcal/mol) -> GROMACS (kJ/mol)
   
    if force_field_mixing == True:
        for b, l2 in enumerate(angle_ff):
            if [atom_A,atom_B,atom_C] in l2:
                deg_harm = l2[3]
                k_harm = l2[4]
                angle_types.append([atom_A, atom_B, atom_C, '1', deg_harm, k_harm])
                break
    else:
        continue
            
    if (len(atom_B) >=3) and (atom_B[2] in [1,2,4,6]):
        if atom_B[2] == 1:
            #linear
            k_harm = K
            deg_harm = 180   
            angle_types.append([atom_A, atom_B, atom_C, '1', deg_harm, k_harm])
        elif atom_B[2] == 2:
            #trigonal
            k_harm = 4*K/3
            deg_harm = 180*2/3
            angle_types.append([atom_A, atom_B, atom_C, '1', deg_harm, k_harm])
        elif (atom_B == 4) or (atom_B == 6):
            #square planar or octahedral
            k_harm = 4*K/3
            deg_harm = 180*2/3
            angle_types.append([atom_A, atom_B, atom_C, '1', deg_harm, k_harm])

    else:
        if angle == 180.0:
            k_harm = K
            deg_harm = 180   
            angle_types.append([atom_A, atom_B, atom_C, '1', deg_harm, k_harm])
        elif angle == 120.0:
            k_harm = 4*K/3
            deg_harm = math.pi*2/3
            angle_types.append([atom_A, atom_B, atom_C, '1', deg_harm, k_harm])

        else: 
            k_harm = K
            deg_harm = angle
            angle_types.append([atom_A, atom_B, atom_C, '1', deg_harm, k_harm])



#
#TORSION BEND

#geting the natural angle of the middle two

dihedral_types = []

for a, line in enumerate(unique_torsions):
    V_tot = None
    nval = None
    _phi = None
    
    k_gro = None
    deg_gro = None

    
    atom_A=line[0]
    atom_B=line[1]
    atom_C=line[2]
    atom_D=line[3]
    
    B_hyb = get_hybridisation(atom_B)
    C_hyb = get_hybridisation(atom_C)
    
    if force_field_mixing == True:
        for x, gamma in enumerate(torsion_ff):
            if ([atom_A,atom_B,atom_C,atom_D] in gamma) or ([atom_D,atom_C,atom_B,atom_A] in gamma):
                deg_gro = l2[4]
                k_gro = l2[5]
                nval = l2[6]
                
            elif (['*',atom_B,atom_C,'*'] in gamma) or (['*',atom_C,atom_B,'*'] in gamma):
                deg_gro = l2[4]
                k_gro = l2[5]
                nval = l2[6]
                dihedral_types.append([atom_A,atom_B, atom_C, atom_D, '1', deg_gro, k_gro, nval])
                break
    else:
        continue
    
    #checking if they have torsion
    if B_hyb in ['sp3', 'sp2'] and C_hyb in ['sp3', 'sp2']:
    
        #only the middle two atom matters
        if (B_hyb == 'sp3') and (C_hyb == 'sp3'): #sp3-sp3 pairs
            
            for elements in g6:    #checking for pair of group 6 atoms
                if elements in atom_B and elements in atom_C:   
                    if (atom_B[0] == 'O') or (atom_C[0] == 'O'): #O(sp3)-g6(sp3)
                        V_tot = math.sqrt(2.0*6.8) #kcal/mol
                        _phi = 0
                        nval = 2
                    else:   #g6(sp3)-g6(sp3)
                        V_tot = 2 #kcal/mol
                        _phi = 0
                        nval = 2
                
                else:       #for all of the other sp3-sp3 bonds
                    nval = 3
                    _phi = 180
                    vj = None
                    vk = None
                    for a, line in enumerate(uff_table['atoms']):
                        if atom_B == uff_table['atoms'][a]:
                            vj = float(uff_table['vsp3'][a])
                    for b, ll in enumerate(uff_table['atoms']):
                        if atom_C == uff_table['atoms'][b]:
                            vk = float(uff_table['vsp3'][b])
                            
                            V_tot=math.sqrt(vj*vk) #kcal/mol
        elif (B_hyb == 'sp2') and (C_hyb == 'sp2'): #sp2-sp2 pair
            vj = None
            vk = None
            nval = 2
            _phi = 180
            for c, l3 in enumerate(uff_table['atoms']):
                if atom_B == uff_table['atoms'][c]:
                    vj = float(uff_table['vsp2'][c])
            for d, l4 in enumerate(uff_table['atoms']):
                if atom_C == uff_table['atoms'][d]:
                    vk = float(uff_table['vsp2'][d])
                    
                    rbo_jk = get_BO(atom_B, atom_C)
                    V_tot = 5*math.sqrt(vj*vk)*(1+4.18*math.log(rbo_jk)) #kcal/mol
        
        elif (B_hyb =='sp2' and C_hyb == 'sp3') or (B_hyb == 'sp3' and C_hyb == 'sp2'): #sp2-sp3
        
            A_hyb = get_hybridisation(atom_A)
            D_hyb = get_hybridisation(atom_D)
            
            if (A_hyb == 'sp2' and B_hyb == 'sp2') or (C_hyb == 'sp2' and D_hyb == 'sp2'): #seeing if the sp2 has another sp2 connection
                V_tot = 2       #kcal/mol
                nval = 3
                _phi = 180
            
            else:
                nval = 6
                _phi = 0
                V_tot = 1       #kcal/mol
    
        k_gro = round(0.5*V_tot*4.184,4)
        deg_gro = nval*_phi-180
        
        dihedral_types.append([atom_A,atom_B, atom_C, atom_D, '1', deg_gro, k_gro, nval])
        
        
#INVERSION

inversion_types = []

if Inversion == True:
    inv_target = ['O_2','O_R', 'N_R', 'N_2']
    
    for a, line in enumerate(unique_inversions):
        _c0 = 0
        _c1 = 0
        _c2 = 0
        _k = 0  #kcal/mol
        k_gro = None #kj/mol
        A_gro = None
        
        atom_A=line[0]
        atom_B=line[1]
        atom_C=line[2]
        atom_D=line[3]
        
        if force_field_mixing == True:
            for b, l2 in enumerate(inversion_ff):
                if ([atom_A,atom_B,atom_C,atom_D] in l2) or ([atom_D,atom_C,atom_B,atom_A] in l2):
                    A_gro = l2[4]
                    k_gro = l2[5]
                    inversion_types.append([atom_A,atom_B, atom_C, atom_D, '4', k_gro, A_gro, '#inversion'])
                    break
        else:
            continue
        
        if (atom_A == 'C_2') or (atom_A == 'C_R'):
            _c0 = 1
            _c1 = -1
            _c2 = 0
            if 'O_2' in line:
                _k = 50 
            else: 
                _k = 6
            if _c2 == 0:
                break
            else:
                k_gro = 4*_c2*(1-(_c1/(4*_c2))**2)*k*4.184 #kj/mol
                A_gro = math.pi- math.acos(_c1/(4*c_c2))
                
        elif atom_A in inv_target:
            _c0 = 1
            _c1 = -1
            _c2 = 0
            _k = 6
            
            if _c2 == 0:
                break
            else:
                k_gro = 4*_c2*(1-(_c1/(4*_c2))**2)*k*4.184 #kj/mol
                A_gro = math.pi- math.acos(_c1/(4*c_c2))
            
            
        else:
            k_gro = 0 #kj/mol
            A_gro = 0
            
        
        inversion_types.append([atom_A,atom_B, atom_C, atom_D, '4', k_gro, A_gro, '#inversion'])

        

        
            


###############################################################################
#force_field.itp file

ff_list = []
ff_list.append(['[ defaults ]'])
ff_list.append(['; nbfunc','comb-rule','gen-pairs','fudgeLJ','fudgeQQ'])
ff_list.append(['1','2','No','0.5','0.5'])
ff_list.append(['[ atomtypes ]'])
ff_list.append([';name','bond_type','mass','charge','ptype','sigma','epsilon'])
ff_list += atom_types
ff_list.append([''])
ff_list.append(['[ bondtypes ]'])
ff_list += bond_types
ff_list.append([''])
ff_list.append(['[ angletypes ]'])
ff_list += angle_types
ff_list.append([''])
ff_list.append(['[ dihedraltypes ]'])
ff_list += dihedral_types
ff_list += inversion_types


ff_name=f'{structure}/force_field.itp'
with open(ff_name, 'w') as file:
    for row in ff_list:
        row_str = '\t'.join(map(str, row)) 
        file.write(row_str + '\n')
        
print("force_field.itp file printed")



#FULLCONNECTIVITY

atom_con = []
for i, row in enumerate(car_table['index']):
    atom_charge = charge_list[i][1]
    atom_con.append([row,car_table['uff names'][i],'1',groupid,car_table['uff names'][i],'1',atom_charge])

bond_con = []
for i, row in enumerate(bond_connectivity):
    atom_A = row[0]
    atom_B = row[1]
    
    index_A = None
    index_B = None
    for a, l1 in enumerate(car_table['labels']):
        if atom_A == l1:
            index_A = car_table['index'][a]
    for b, l2 in enumerate(car_table['labels']):
        if atom_B == l2:
            index_B = car_table['index'][b]
    
    bond_con.append([index_A, index_B, '1'])
    
angle_con = []
for i, row in enumerate(angle_connectivity):
    atom_A = row[0]
    atom_B = row[1]
    atom_C = row[2]
    
    index_A = None
    index_B = None
    index_C = None
    for a, l1 in enumerate(car_table['labels']):
        if atom_A == l1:
            index_A = car_table['index'][a]
    for b, l2 in enumerate(car_table['labels']):
        if atom_B == l2:
            index_B = car_table['index'][b]
    for c, l3 in enumerate(car_table['labels']):
        if atom_C == l3:
            index_C = car_table['index'][c]
    
    angle_con.append([index_A, index_B,index_C, '1'])
    
torsion_con = []
for i, row in enumerate(torsion_connectivity):
    atom_A = row[0]
    atom_B = row[1]
    atom_C = row[2]
    atom_D = row[3]
    
    index_A = None
    index_B = None
    index_C = None
    index_D = None
    for a, l1 in enumerate(car_table['labels']):
        if atom_A == l1:
            index_A = car_table['index'][a]
    for b, l2 in enumerate(car_table['labels']):
        if atom_B == l2:
            index_B = car_table['index'][b]
    for c, l3 in enumerate(car_table['labels']):
        if atom_C == l3:
            index_C = car_table['index'][c]
    for d, l3 in enumerate(car_table['labels']):
        if atom_D == l3:
            index_D = car_table['index'][d]
    
    torsion_con.append([index_A, index_B,index_C,index_D,'1'])
    
inversion_con = []
for i, row in enumerate(inversion_connectivity):
    atom_A = row[0]
    atom_B = row[1]
    atom_C = row[2]
    atom_D = row[3]
    
    index_A = None
    index_B = None
    index_C = None
    index_D = None
    for a, l1 in enumerate(car_table['labels']):
        if atom_A == l1:
            index_A = car_table['index'][a]
    for b, l2 in enumerate(car_table['labels']):
        if atom_B == l2:
            index_B = car_table['index'][b]
    for c, l3 in enumerate(car_table['labels']):
        if atom_C == l3:
            index_C = car_table['index'][c]
    for d, l3 in enumerate(car_table['labels']):
        if atom_D == l3:
            index_D = car_table['index'][d]
    
    inversion_con.append([index_A, index_B,index_C,index_D,'4'])
            
con_list = []
con_list.append(['[ moleculetype ]'])
con_list.append([';','','Name','nrexcl'])
con_list.append([f'{groupid}','','3'])
con_list.append([])
con_list.append(['[ atoms ]'])
con_list += atom_con
con_list.append([])
con_list.append(['[ bonds ]'])
con_list += bond_con
con_list.append([])
con_list.append(['[ angles ]'])
con_list += angle_con
con_list.append([])
con_list.append(['[ dihedrals ]'])
con_list += torsion_con
if Inversion == True:
    con_list.append(['[ dihedrals ]'])
    con_list += inversion_con

itpname = f'{structure}/{groupid}.itp'
with open(itpname, 'w') as file:
    for row in con_list:
        row_str = '\t'.join(map(str, row)) 
        file.write(row_str + '\n')

print(f"Data written to {itpname}")
            

   
    
