#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 11:23:13 2024
Generalised UFF force field conversion code to generate input suitable for RASPA
to get the connectivity for all types of framework
UFF_table obtained from CP2K DTFB section
"""
from functions import *
start = datetime.now()
print(start)

#input needed
structure= 'OHIHET'
#connectivity we want
bond_connectivity=True
angle_connectivity=True
torsion_connectivity=True
inv = 0
improper = 0


#####
inp = f'{structure}/input/file3.data' #the pdb file name, currently we're using the lmp file name from VMD as the first half of the code is not complete
pdb = f'{structure}/input/{structure}.pdb' #the pdb file name, currently we're using the lmp file name from VMD as the first half of the code is not complete



#Reading the UFF table data
_uff='data/uff_table'
with open(_uff,'r') as f:
    lines = f.readlines()
uff = []
for line in lines:
    words = [x.strip() for x in line.split(',')]
    uff.append(words)

uff_table = pd.DataFrame(uff, columns=['atoms', 'ri', 'phi', 'xi', 'di', 'psi', 'zmm', 'vsp3', 'vsp2', 'chi', 'nc'])
###############################################################################
#Reading the PDB file
#reading the data from writelammpsdata from VMD
with open(pdb,'r') as f:
    lines = f.readlines()
pdb_ = []
for line in lines:
    word = line.split()
    pdb_.append(word)
    
#The order of the atoms
atom_list=[]
target='ATOM'
for i,row in enumerate(pdb_):
    if target in row:
        atom_list.append([pdb_[i][1],pdb_[i][10],pdb_[i][2]])
        
unique_atoms=[]
for i, row in enumerate(atom_list):
    _atom = row[2]
    if _atom not in unique_atoms:
        unique_atoms.append(_atom)
print(f"there are {len(unique_atoms)} different atom types: {unique_atoms}")

##############################################################################
#reading the data from writelammpsdata from VMD
with open(inp,'r') as f:
    lines = f.readlines()
lmp = []
for line in lines:
    word = line.split()
    lmp.append(word)

target='atoms'
for i,row in enumerate(lmp):
    if target in row:
        atom_num=lmp[i][0]

target='bonds'
for i,row in enumerate(lmp):
    if target in row:
        bond_num=lmp[i][0]

target='angles'
for i,row in enumerate(lmp):
    if target in row:
        angle_num=lmp[i][0]
        
target='dihedrals'
for i,row in enumerate(lmp):
    if target in row:
        dihedral_num=lmp[i][0]
        
#the parameters
target = 'Masses' 
mass = []
for i,row in enumerate(lmp):
    if target in row:
        #line = row.index(target)
        mass=lmp[i+2:i+2+len(unique_atoms)]

print(f"there are {atom_num} atoms, {bond_num} bonds, {angle_num} angles and {dihedral_num} torsions observed by VMD")

###############################################################################
#getting separate list of bonds
    #writelammpsdata does not give the bonds for some reason
if bond_connectivity==True:
    target='Bonds'
    for i, row in enumerate(lmp):
        if target in row:
            _start=i+2
            _bondlist=lmp[_start:_start+int(bond_num)]
            bondlist=[[i[2], i[3]] for i in _bondlist]
    
    #assigning the lists of indexes to their respective atoms and remove the duplicates
    unique_bonds = []
    
    for i, row in enumerate(bondlist):
        _A = bondlist[i][0]
        _B = bondlist[i][1]
        atom_A = None
        atom_B = None
        
        for n, line in enumerate(atom_list):
            if _A in line:
                atom_A = line[2]
        for x, l in enumerate(atom_list):
            if _B in l:
                atom_B = l[2]
                newbond = [atom_A, atom_B]
    
                if newbond not in unique_bonds:
                    unique_bonds.append(newbond)
    print(f"there are {len(unique_bonds)} unique bonds")

if bond_connectivity==False:
    print('0 bond connectivity searched')

#getting the list of angles
if angle_connectivity==True:
    target = 'Angles'
    for i, row in enumerate(lmp):
        if target in row:
            _start=i+2
            _anglelist=lmp[_start:_start+int(angle_num)]
            anglelist=[[i[2], i[3],i[4]] for i in _anglelist]

    unique_angles=[]
    for i, row in enumerate(anglelist):
        _A = anglelist[i][0]
        _B = anglelist[i][1]
        _C = anglelist[i][2]
        atom_A = None
        atom_B = None
        atom_C = None
        
        for n, line in enumerate(atom_list):
            if _A in line:
                atom_A = line[2]
        for x, l in enumerate(atom_list):
            if _B in l:
                atom_B = l[2]       
        for y, sub in enumerate(atom_list):
            if _C in sub:
                atom_C = sub[2]
                newangle = [atom_A, atom_B, atom_C]
    
                if newangle not in unique_angles:
                    unique_angles.append(newangle)
    print(f"there are {len(unique_angles)} unique angles")
    
if angle_connectivity==False:
    print('0 angle connectivity searched')

#getting the list of dihedrals
if torsion_connectivity==True:
    target = 'Dihedrals'
    for i, row in enumerate(lmp):
        if target in row:
            _start=i+2
            _dihedrallist=lmp[_start:_start+int(dihedral_num)]
            dihedrallist=[[i[2], i[3],i[4],i[5]] for i in _dihedrallist]

    unique_torsions=[]
    for i, row in enumerate(dihedrallist):
        _A = dihedrallist[i][0]
        _B = dihedrallist[i][1]
        _C = dihedrallist[i][2]
        _D = dihedrallist[i][3]
        atom_A = None
        atom_B = None
        atom_C = None
        atom_D = None
        
        for n, line in enumerate(atom_list):
            if _A in line:
                atom_A = line[2]
        for x, l in enumerate(atom_list):
            if _B in l:
                atom_B = l[2]       
        for y, sub in enumerate(atom_list):
            if _C in sub:
                atom_C = sub[2]
        for dd, ll in enumerate(atom_list):
            if _D in ll:
                atom_D = ll[2]
                newtorsion = [atom_A, atom_B, atom_C, atom_D]
    
                if newtorsion not in unique_torsions:
                    unique_torsions.append(newtorsion)
    print(f"there are {len(unique_torsions)} unique torsions")
    
if torsion_connectivity==False:
    print('0 torsion connectivity searched')


#Making the pseudo atom files
pseudoatom=[]
for n in unique_atoms:
    _atom=n
    for x, ll in enumerate(mass):
        if _atom in ll:
            _mass=ll[1]
    for i, line in enumerate(uff_table['atoms']):
        if _atom == uff_table['atoms'][i]:
            pseudoatom.append([_atom,'yes',_atom,_atom,'0',_mass,uff_table['zmm'][i],
                               0.0,1.0,uff_table['ri'][i],uff_table['nc'][i],'absolute',0])
        
pseudo_atom=[]
pseudo_atom.append(['#number','of','pseudo','atoms'])
pseudo_atom.append(str(len(unique_atoms)))
pseudo_atom.append(['#type','print','as','chem','oxidation','mass','charge',
                    'polarization','B-factor','radii','connectivity','anisotropic-type',
                    'tinker-type'])
pseudo_atom += pseudoatom

paff_name=f'{structure}/output/pseudo_atoms.def'
with open(paff_name, 'w') as file:
    for row in pseudo_atom:
        row_str = '\t'.join(map(str, row)) 
        file.write(row_str + '\n')
        
print(f"Pseudo atoms Data written to {paff_name} in the output file")


#Lennard Jones Potential
    #in UFF units
LJ=[]
for n in unique_atoms:
    
    _atom=n
    for i, line in enumerate(uff_table['atoms']):
        if _atom == uff_table['atoms'][i]:
            _energy=round(float(uff_table['di'][i])*503.22271716452354,2)
            _radius=round((float(uff_table['xi'][i]))*0.89089871814,2)
            LJ.append([_atom,' ','LENNARD_JONES',_energy,_radius])
            
ff_mixing=[]
ff_mixing.append(['# general rule for shifted vs truncated'])
ff_mixing.append(['shifted'])
ff_mixing.append(['# general rule tailcorrections'])
ff_mixing.append(['no'])
ff_mixing.append(['# number of defined interactions'])
ff_mixing.append([str(len(unique_atoms))])
ff_mixing.append(['# type interaction, parameters.',
                  'IMPORTANT: define shortest matches first, so that more specific ones overwrites these'])
ff_mixing += LJ
ff_mixing.append(['# general mixing rule for Lennard-Jones'])
ff_mixing.append(['Lorentz-Berthelot'])

ffmr_name=f'{structure}/output/force_field_mixing_rules.def'
with open(ffmr_name, 'w') as file:
    for row in ff_mixing:
        row_str = '\t'.join(map(str, row)) 
        file.write(row_str + '\n')
        
print(f"Lennard Jones parameters written to {ffmr_name} in the output file")

##############################################################################
#BOND Order calculation
#assigning hybridisation to each atoms
atom_hybridisation=pd.DataFrame()
hyb_list = []
for i in range(len(unique_atoms)):
    hyb_list.append(get_hybridisation(unique_atoms[i]))

atom_hybridisation['atoms']=unique_atoms
atom_hybridisation['hybridisation']=hyb_list

#now assigning bond orders for each pairs
BO_list=[]
for i in range(len(unique_bonds)):
    atom_A = unique_bonds[i][0]
    atom_B = unique_bonds[i][1]
    BO_list.append([atom_A,atom_B,get_BO(atom_A, atom_B)])

            
##############################################################################
#BOND STRETCH

#DISTANCE: calculating the r_ij value 
    #r_ij = r_i + r_j + r_BO - r_EN
r_ij=[]

for i in range(len(unique_bonds)):
    atom_A = unique_bonds[i][0]
    atom_B = unique_bonds[i][1]
    
    r_ij.append([atom_A, atom_B, get_natural_bond(atom_A, atom_B)])
                
    
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
            
            k_ab=664.12*((z_i*z_j)/(rab)**3)
            k_ij.append([atom_A, atom_B,k_ab])

#converting the UFF units to the RASPA units
k_ij_R = []
for i in range(len(k_ij)):
    kab=k_ij[i][2]*503.22271716452354
    kab=round(kab,2)
    k_ij_R.append(kab)
    
raspa_bond=[]
for i in range(len(k_ij)):
    raspa_bond.append([k_ij[i][0],k_ij[i][1],'HARMONIC_BOND',k_ij_R[i],r_ij[i][2]])

#creating the bond stretch section in RASPA framework.def file
bond_stretch_R=[]
bond_stretch_R.append(['#bond stretch atom n1-n2, equilibrium distance, bondforce-constant'])
bond_stretch_R += raspa_bond


##############################################################################
#Angle Bend
#obtaining the middle atom natural angle
theta=[]
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
    _r_ik=get_natural_bond(atom_A, atom_C)  #atom 1 and 3
    
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
            
            
            #K_ijk equation
            _beta = 664.12/(_r_ij*_r_jk)
            
            _back = 3*_r_ij*_r_jk*(1-((math.cos(theta_rad))**2))-(((_r_ik)**2)*math.cos(theta_rad))
            
            _kijk= _beta*((z_i*z_k)/(_r_ik)**5)*(_r_ij*_r_jk)*_back
            
            K_ijk.append([atom_A, atom_B, atom_C,_kijk])
            
            
#final equation
angle_type1_list = []   #linear type: atom A, B, C, P0, P1, P2 (already in RASPA units)
angle_type2_list = []   #nonlinear type: atom A, B, C, K, C0, C1, C2
angle_list_R = []
#atom A, atom B, atom C, P0, P1, P2

for a, line in enumerate(unique_angles):
    atom_A=line[0]
    atom_B=line[1]
    atom_C=line[2]

    #the terms needed to calculate the coefficients
    K = None
    angle = None
    theta_rad = None
    _C0 = None
    _C1 = None
    _C2 = None
    _p0 = None
    _p1 = None
    _p2 = None
    
    for n, row in enumerate(theta):
        if atom_B == row[0]:
            angle=float(row[1])
            theta_rad = math.radians(angle)
            
    for i, sub in enumerate(K_ijk):
        if atom_B == K_ijk[i][1]:
            K = round(503.22271716452354*K_ijk[i][3],2) #UFF (Kcal/mol) -> RASPA (K)
            
    if (len(atom_B) >=3) and (atom_B[2] in [1,2,4,6]):
        if atom_B[2] == 1:
            _n = 1 #linear
        elif atom_B[2] == 2:
            _n = 3 #trigonal
        elif (atom_B == 4) or (atom_B == 6):
            _n = 4 #square planar or octahedral
            
            _p0 = K/(_n)**2 #in KELVIN
            _p1 = _n #dimensionless
            _p2 = 180.0 #degrees
            
            angle_type1_list.append([atom_A,atom_B,atom_C,_p0,_p1,_p2])
            angle_list_R.append([atom_A,atom_B,atom_C,'COSINE_BEND',round(_p0,2),round(_p1,2),round(_p2,2)])
    
    else:
        if angle == 180.0:
            _n = 1 #linear
            _p0 = K/(_n)**2 #in KELVIN
            _p1 = _n #dimensionless
            _p2 = 180.0 #degrees
            angle_type1_list.append([atom_A,atom_B,atom_C,'COSINE_BEND',_p1,_p2])
            angle_list_R.append([atom_A,atom_B,atom_C,'COSINE_BEND',round(_p0,2),round(_p1,2),round(_p2,2)])

        elif angle == 120.0:
            _n = 3
            _p0 = K/(_n)**2 #in KELVIN
            _p1 = _n #dimensionless
            _p2 = 180.0 #degrees
            angle_type1_list.append([atom_A,atom_B,atom_C,'COSINE_BEND',_p1,_p2])
            angle_list_R.append([atom_A,atom_B,atom_C,'COSINE_BEND',round(_p0,2),round(_p1,2),round(_p2,2)])

        else: 
            _c2 = 1/(4*(math.sin(theta_rad))**2)
            _c1 = -4*_c2*math.cos(theta_rad)
            _c0 = _c2*((2*(math.cos(theta_rad))**2)+1)
            
            angle_type2_list.append([atom_A,atom_B,atom_C,K,_c0,_c1,_c2])

#General nonlinear type conversion to raspa
for a, line in enumerate(angle_type2_list):
    atom_A = line[0]
    atom_B = line[1]
    atom_C = line[2]
    K = line[3]
    C_0 = line[4]
    C_1 = line[5]
    C_2 = line[6]
    _p0 = K*C_2 #K
    _p1 = 2 #dimensionless
    _p2 = 0 #angle
    angle_list_R.append([atom_A,atom_B,atom_C,'COSINE_BEND',round(_p0,2),round(_p1,2),round(_p2,2)])
    _b0 = K*C_1
    _b1 = 1 #dimensionless
    _b2 = 0 #degrees
    angle_list_R.append([atom_A,atom_B,atom_C,'COSINE_BEND',round(_b0,2),round(_b1,2),round(_b2,2)])
    _x0 = K*(C_0-C_2-C_1)
    _x1 = 0
    _x2 = -90 #degrees
    angle_list_R.append([atom_A,atom_B,atom_C,'COSINE_BEND',round(_x0,2),round(_x1,2),round(_x2,2)])
        
#######################################################
#torsion_E_phi
#geting the natural angle of the middle two

torsion_list_R = []

for a, line in enumerate(unique_torsions):
    V_tot = None
    nval = None
    _phi = None
    
    p_0 = None
    p_1 = None
    p_2 = None
    p_3 = None
    p_4 = None
    p_5 = None
    
    atom_A=line[0]
    atom_B=line[1]
    atom_C=line[2]
    atom_D=line[3]
    
    B_hyb = get_hybridisation(atom_B)
    C_hyb = get_hybridisation(atom_C)
    
    #checking if they have torsion
    if B_hyb in ['sp3', 'sp2'] and C_hyb in ['sp3', 'sp2']:
    
        #only the middle two atom matters
        if (B_hyb == 'sp3') and (C_hyb == 'sp3'): #sp3-sp3 pairs
            
            for elements in g6:    #checking for pair of group 6 atoms
                if elements in atom_B and elements in atom_C:   
                    if (atom_B[0] == 'O') or (atom_C[0] == 'O'): #O(sp3)-g6(sp3)
                        V_tot = math.sqrt(2.0*6.8)
                        _phi = 0
                        nval = 2
                    else:   #g6(sp3)-g6(sp3)
                        V_tot = 2
                        _phi = 0
                        nval = 2
                
                else:       #for all of the other sp3-sp3 bonds
                    nval = 3
                    _phi = math.pi
                    vj = None
                    vk = None
                    for a, line in enumerate(uff_table['atoms']):
                        if atom_B == uff_table['atoms'][a]:
                            vj = float(uff_table['vsp3'][a])
                    for b, ll in enumerate(uff_table['atoms']):
                        if atom_C == uff_table['atoms'][b]:
                            vk = float(uff_table['vsp3'][b])
                            
                            V_tot=math.sqrt(vj*vk)
        elif (B_hyb == 'sp2') and (C_hyb == 'sp2'): #sp2-sp2 pair
            vj = None
            vk = None
            nval = 2
            _phi = math.pi
            for c, l3 in enumerate(uff_table['atoms']):
                if atom_B == uff_table['atoms'][c]:
                    vj = float(uff_table['vsp2'][c])
            for d, l4 in enumerate(uff_table['atoms']):
                if atom_C == uff_table['atoms'][d]:
                    vk = float(uff_table['vsp2'][d])
                    
                    rbo_jk = get_BO(atom_B, atom_C)
                    V_tot = 5*math.sqrt(vj*vk)*(1+4.18*math.log(rbo_jk))
        
        elif (B_hyb =='sp2' and C_hyb == 'sp3') or (B_hyb == 'sp3' and C_hyb == 'sp2'): #sp2-sp3
        
            A_hyb = get_hybridisation(atom_A)
            D_hyb = get_hybridisation(atom_D)
            
            if (A_hyb == 'sp2' and B_hyb == 'sp2') or (C_hyb == 'sp2' and D_hyb == 'sp2'): #seeing if the sp2 has another sp2 connection
                V_tot = 2
                nval = 3
                _phi = math.pi
            
            else:
                nval = 6
                _phi = 0
                V_tot = 1
    
        if nval == 2:
            C_0 = 1
            K = 0.5*V_tot*503.22271716452354
            C_1 = -math.cos(_phi)
            C_2 = -math.cos(2*_phi)
            
            p_0 = round((K*(C_0-C_2)),2)
            p_1 = round((-K*C_1),2)
            p_2 = round((2*K*C_2),2)
            p_3 = 0
            p_4 = 0
            p_5 = 0
            
            torsion_list_R.append([atom_A,atom_B,atom_C,atom_D,'SIX_COSINE_DIHEDERAL',p_0,p_1,p_2,p_3,p_4,p_5])

        elif nval == 3:
            C_0 = 1
            K = 0.5*V_tot*503.22271716452354
            C_1 = -math.cos(_phi)
            C_2 = -math.cos(2*_phi)
            C_3 = -math.cos(3*_phi)
            
            p_0 = round((K*(C_0-C_2)),2)
            p_1 = round((-K*C_1),2)
            p_2 = round((2*K*C_2),2)
            p_3 = round((-4*K*C_3),2)
            p_4 = 0
            p_5 = 0
            
            torsion_list_R.append([atom_A,atom_B,atom_C,atom_D,'SIX_COSINE_DIHEDERAL',p_0,p_1,p_2,p_3,p_4,p_5])

            
        elif nval == 6:
            C_0 = 1
            K = 0.5*V_tot*503.22271716452354
            C_1 = -math.cos(_phi)
            C_2 = -math.cos(2*_phi)
            C_3 = -math.cos(3*_phi)
            C_4 = -math.cos(4*_phi)
            C_5 = -math.cos(5*_phi)
            C_6 = -math.cos(6*_phi)
            
            p_0 = round((K*(C_0-C_2+C_4-33*C_6)),2)
            p_1 = round((-K*(C_1-3*C_3+5*C_5)),2)
            p_2 = round((2*K*(2*C_2-8*C_4+18*C_6)),2)
            p_3 = round((-4*K*(4*C_3-20*C_5)),2)
            p_4 = round((K*(8*C_4-48*C_6)),2)
            p_5 = round((16*K*C_5),2)
            b_0 = round((32*K*C_6),2)
            b_1 = 6
            b_2 = 0
            
            torsion_list_R.append([atom_A,atom_B,atom_C,atom_D,'SIX_COSINE_DIHEDERAL',p_0,p_1,p_2,p_3,p_4,p_5])
            torsion_list_R.append([atom_A,atom_B,atom_C,atom_D,'CVFF_DIHEDRAL',b_0,b_1,b_2])
    
    else:
        p_0 = 0 
        p_1 = 0
        p_2 = 0
        p_3 = 0
        p_4 = 0
        p_5 = 0

##############################################################################
#creating the framework.def file
frame_def = []
frame_def.append(['#CoreShells bond  BondDipoles UreyBradley bend  inv  tors improper-torsion bond/bond bond/bend bend/bend stretch/torsion bend/torsion'])
frame_def.append([0,len(unique_bonds), 0, 0, 0, len(angle_list_R), inv, len(torsion_list_R), improper, 0, 0, 0, 0, 0])
frame_def += bond_stretch_R
frame_def.append(['#bond bending atom n1-n2-n3, equilibrium angle, bondforce-constant'])
frame_def += angle_list_R
frame_def.append(['#torsion atom n:Q1-n2-n3-n4, potential, coefficients'])
frame_def += torsion_list_R

framedef_name=f'{structure}/output/framework.def'
with open(framedef_name, 'w') as file:
    for row in frame_def:
        row_str = '\t'.join(map(str, row)) 
        file.write(row_str + '\n')
        
print(f"Framework connectivity and parameters written to {ffmr_name} in the output file")
print("--- %s time taken ---" % ((datetime.now() - start)))
                    
                
            
            
     

    
    
  
        


