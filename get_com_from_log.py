#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 16:16:53 2024

@author: samantha.bealldenn
"""
import numpy as np
import argparse
from openbabel import pybel
from get_xyz_from_log import log_to_xyz

def get_arguments():
    parser=argparse.ArgumentParser()
    #Adds essential information as arguments
    #gets filename from command line
    parser.add_argument("filename", nargs=1, action='store')
    #gets density functional from command line
    parser.add_argument("density_functional", choices=["b3lyp", "mn15", "b2plypd3", "b2plyp"], nargs=1, action="store")
    #gets field stength from command line
    parser.add_argument("field_strength", choices=["n4", "n3", "n2", "n1", 'nofield', 'p1', 'p2', 'p3', 'p4'], nargs=1, action="store")
    #gets solvent from command line
    parser.add_argument("solvent", choices=["acn", "dcm"], nargs=1, action='store')
    parser.add_argument("basis_set", choices=["+", "++", "pvdz", "pvtz"], nargs=1, action="store")
    
    group = parser.add_mutually_exclusive_group(required=False)
    #changes the regular minimization to a single point calculation
    group.add_argument("-s", "--spc", action="store_true")
    #adds the option to make the output file QST3, in which case we want to get the names of the other two files being used
    group.add_argument("-q", "--qst3", nargs=2, action='store')
    
    parser.add_argument("-f", "--freq", action="store_true")
    #adds the option to include frequency calculations, which factor entropy into the free energy calculations

    args=parser.parse_args()
    return args
    #args is the entire namespace

#density functional options
def get_density_functional(args, options):
    if args.density_functional.__contains__('b3lyp'):
        options['density_functional'] = "B3LYP"
        options['empirical_dispersion'] = "EmpiricalDispersion=GD3"
    elif args.density_functional.__contains__('mn15'):
        options['density_functional'] = "MN15"
        options['empirical_dispersion'] = ""
    elif args.density_functional.__contains__('b2plypd3'):
        options['density_functional'] = "B2PLYPD3"
        options['empirical_dispersion'] = ""
    elif args.density_functional.__contains__('b2plyp'):
        options['density_functional'] = "B2PLYP"
        options['empirical_dispersion'] = ""
    return options

    
#field strength options    
def get_field_strength(args, options):
    if args.field_strength == ['n4']:
        options['field_strength'] = "Field=Z-4"
    elif args.field_strength == ['n3']:
        options['field_strength'] = "Field=Z-3"
    elif args.field_strength == ['n2']:
        options['field_strength'] = "Field=Z-2"
    elif args.field_strength == ['n1']:
        options['field_strength'] = "Field=Z-1"
    elif args.field_strength == ['nofield']:
        options['field_strength'] = ""
    elif args.field_strength == ['p1']:
        options['field_strength'] = "Field=Z+1"
    elif args.field_strength == ['p2']:
        options['field_strength'] = "Field=Z+2"
    elif args.field_strength == ['p3']:
        options['field_strength'] = "Field=Z+3"
    elif args.field_strength == ['p4']:
        options['field_strength'] = "Field=Z+4"
    return options
    
#solvent options
def get_solvent(args, options):
    if args.solvent == ['acn']:
        options['solvent'] = "SCRF=(Solvent=Acetonitrile)"
    elif args.solvent == ['dcm']:
        options['solvent'] = "SCRF=(Solvent=Dichloromethane)"
    else:
        print("Unrecognized solvent")
    return options

def get_basis_set(args, options):
    if args.basis_set == ['+']:
        options['basis_set'] = "6-311+g(d,p)"
    elif args.basis_set == ['++']:
        options['basis_set'] = "6-311++g(d,p)"
    elif args.basis_set == ['pvdz']:
            options['basis_set'] = "AUG-cc-pVDZ"  
    elif args.basis_set == ['pvtz']:
            options['basis_set'] = "AUG-cc-pVTZ"
    return options

#turns the xyz file format into a z matrix format in a string
def get_coordinates(options, filename):
    with open(filename, 'r') as f:
        file_string=f.read()
    molecule = pybel.readstring("xyz", file_string)
    molecule_gzmat = molecule.write(format="gzmat")
    molecule_gzmat_split=molecule_gzmat.split("\n")
    del molecule_gzmat_split[0:5]
    del molecule_gzmat_split[-1]
    coordinates="\n".join(molecule_gzmat_split)
    return coordinates

#gets all the options for the name of the .com file
def get_dotcom_filename(args, options, filename_options):
    filename=args.filename[0]
    #gets product/reactant, ketone/aldehyde, enantiomer, and long/short arm parameters from the old file name
    if filename.__contains__("product"):
        filename_options['molecule_type'] = "_product"
    elif filename.__contains__("reactant"):
        filename_options['molecule_type'] = "_reactant"
    else:   
        filename_options['molecule_type'] = "_qst3"
   
    ## adds aldehyde and ketone to the name of the file
    ## in the case of QST3 files, adds aldehyde and ketone to the file    
    if filename.__contains__("ketone"):
        filename_options['pathway'] = "_ketone"
        options['product_type'] = "Ketone"
    elif filename.__contains__("aldehyde"):
        filename_options['pathway'] = "_aldehyde"
        options['product_type'] = "Aldehyde"
    else:
        filename_options['pathway'] = ""
        
    if filename.__contains__("_s_"):
        filename_options['enantiomer'] = "_s"
    elif filename.__contains__("_r_"):
        filename_options['enantiomer'] = "_r"
    else:
        filename_options['enantiomer'] = ""
        
    if filename.__contains__("_longarm_"):
        filename_options['reactant_type'] = "_longarm"
    elif filename.__contains__("_shortarm_"):
        filename_options['reactant_type'] = "_shortarm"
    else:
        filename_options['reactant_type'] = ""

#gets the density functional, solvent, and field strength options from the command line
    if args.density_functional == ['b3lyp']:
        filename_options['density_functional'] = '_b3lyp'
    elif args.density_functional == ['mn15']:
        filename_options['density_functional'] = '_mn15'
    elif args.density_functional == ['b2plyp']:
        filename_options['density_functional'] = '_b2plyp'
    elif args.density_functional == ['b2plypd3']:
        filename_options['density_functional'] = '_b2plypd3'     
    
    if args.solvent == ['acn']:
        filename_options['solvent'] = '_acn'
    elif args.solvent == ['dcm']:
        filename_options['solvent'] = '_dcm'
        
        
    if args.basis_set == ['+']:
        filename_options['basis_set'] = ""
    if args.basis_set == ['++']:
        filename_options['basis_set'] = "_++"
    if args.basis_set == ['pvdz']:
        filename_options['basis_set'] = "_pVDZ"  
    if args.basis_set == ['pvtz']:
        filename_options['basis_set'] = "_pVTZ"  
    
    if args.field_strength == ['n4']:
        filename_options['field_strength'] = "_neg4"
    elif args.field_strength == ['n3']:
        filename_options['field_strength'] = "_neg3"
    elif args.field_strength == ['n2']:
        filename_options['field_strength'] = "_neg2"
    elif args.field_strength == ['n1']:
        filename_options['field_strength'] = "_neg1"
    elif args.field_strength == ['nofield']:
        filename_options['field_strength'] = "_nofield"
    elif args.field_strength == ['p1']:
        filename_options['field_strength'] = "_pos1"
    elif args.field_strength == ['p2']:
        filename_options['field_strength'] = "_pos2"
    elif args.field_strength == ['p3']:
        filename_options['field_strength'] = "_pos3"
    elif args.field_strength == ['p4']:
        filename_options['field_strength'] = "_pos4"
    return filename_options

def get_dot_com(args, options, filename_options):
    #JUST reactant or product
    if args.spc == False and args.qst3 == None:
        if args.freq == False:
            syntax="%mem=24GB\n%NProcShared=32\n#n opt=Z-Matrix NoSymm {density_functional}/{basis_set} {empirical_dispersion} {field_strength} {solvent}\n\n {filename}\n\n{coordinates}\n".format_map(options)
            completed_filename=("cis-stilbene_oxide{molecule_type}{density_functional}{basis_set}{solvent}{pathway}{enantiomer}{reactant_type}{field_strength}.com").format_map(filename_options)
        elif args.freq == True:
            syntax="%mem=24GB\n%NProcShared=32\n#n opt=Z-Matrix NoSymm {density_functional}/{basis_set} {empirical_dispersion} {field_strength} {solvent} Freq\n\n {filename}\n\n{coordinates}\n".format_map(options)
            completed_filename=("cis-stilbene_oxide{molecule_type}{density_functional}{basis_set}{solvent}{pathway}{enantiomer}{reactant_type}{field_strength}_freq.com").format_map(filename_options)        
    elif args.spc == True:
        if args.freq == False:
            syntax="%mem=24GB\n%NProcShared=32\n#n NoSymm {density_functional}/{basis_set} {empirical_dispersion} {field_strength} {solvent} Polar\n\n {filename}\n\n{coordinates}\n".format_map(options)
            completed_filename=("cis-stilbene_oxide{molecule_type}{density_functional}{basis_set}{solvent}{pathway}{enantiomer}{reactant_type}{field_strength}_spc_polar.com").format_map(filename_options)
        elif args.freq == True:
            syntax="%mem=24GB\n%NProcShared=32\n#n NoSymm {density_functional}/{basis_set} {empirical_dispersion} {field_strength} {solvent} Freq\n\n {filename}\n\n{coordinates}\n".format_map(options)
            completed_filename=("cis-stilbene_oxide{molecule_type}{density_functional}{basis_set}{solvent}{pathway}{enantiomer}{reactant_type}{field_strength}_spc_freq.com").format_map(filename_options)
    elif args.qst3 != None:
        if args.freq == False:
            syntax="%mem=24GB\n%NProcShared=32\n#n opt=(Z-Matrix,QST3,calcfc) NoSymm {density_functional}/{basis_set} {empirical_dispersion} {field_strength} {solvent}\n\nStarting Material\n\n{starting_material}\n{product_type} Product\n\n{product}\nSaddle Point Guess\n\n{saddle_point}\n".format_map(options)
            completed_filename=("cis-stilbene_oxide{molecule_type}{density_functional}{basis_set}{solvent}{pathway}{enantiomer}{reactant_type}{field_strength}.com").format_map(filename_options)
        elif args.freq == True:
            syntax="%mem=24GB\n%NProcShared=32\n#n opt=(Z-Matrix,QST3,calcfc) NoSymm {density_functional}/{basis_set} {empirical_dispersion} {field_strength} {solvent} Freq\n\nStarting Material\n\n{starting_material}\n{product_type} Product\n\n{product}\nSaddle Point Guess\n\n{saddle_point}\n".format_map(options)
            completed_filename=("cis-stilbene_oxide{molecule_type}{density_functional}{basis_set}{solvent}{pathway}{enantiomer}{reactant_type}{field_strength}_freq.com").format_map(filename_options)
    print("Prepared .com file name:", completed_filename)     
    with open(completed_filename, 'w') as f:
        f.write(syntax)

def read_xyz(filename):
    """ Read in an xyz file.  Returns a list of the atom column and 
    numpy array containing the atomic coordinates
    """
    f = open(filename,'r')

    n_atoms = int(f.readline())
    f.readline()
    atoms_col = []
    coords = [[],[],[]]

    for line in f:
        words = line.split()
        atoms_col.append(words[0])
        coords[0].append(float(words[1]))
        coords[1].append(float(words[2]))
        coords[2].append(float(words[3]))
    return atoms_col, np.array(coords).T

def write_xyz(filename,atom_col,coords,write_file=True):
    """ Write an xyz file, given a list containing the atom names and a
    numpy array containing the atomic coordinates
    """
    coords_to_write = coords.T
    natoms = len(atom_col)
    file_str = ""
    file_str += str(natoms) + '\n'
    file_str += 'generated by a python script \n'
    if write_file:
        f = open(filename,'w')
        f.write(file_str)

    for i in range(natoms):
        atom = atom_col[i]
        x = str(coords_to_write[0][i])
        y = str(coords_to_write[1][i])
        z = str(coords_to_write[2][i])
        sep = '     '
        line = atom + sep + x + sep + y + sep + z + '\n'
        file_str += line
        if write_file:
            f.write(line)
    return file_str


def normalize(v):
    """Normalize a vector."""
    norm = np.linalg.norm(v)
    if norm == 0: 
        return v
    return v / norm

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation 
    about the given axis by theta radians.
    """
    axis = normalize(axis)
    a = np.cos(theta / 2.0)
    b, c, d = -axis * np.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

def align_bond_to_z_axis(coords):
    bond_vector = coords[1] - coords[0]     # Calculate bond vector
    z_axis = np.array([0, 0, 1])            # Z-axis vector
    
    axis = np.cross(bond_vector, z_axis)    # Axis of rotation 
    angle = np.arccos(np.dot(normalize(bond_vector), z_axis))
    
    R = rotation_matrix(axis, angle)        # Rotation matrix
    
    # Apply rotation to all coordinates
    # transpose R for right-hand multiplication
    rotated_coords = np.dot(coords, R.T)  
    
    return rotated_coords

def rotation_matrix_z(theta):
    """Return the rotation matrix for a counterclockwise rotation around 
    the z-axis by theta radians."""
    return np.array([[np.cos(theta), -np.sin(theta), 0],
                     [np.sin(theta),  np.cos(theta), 0],
                     [0,              0,             1]])

def align_bond_to_x_axis(coords):
    # Calculate the bond vector between the first and third atoms
    bond_vector = coords[2] - coords[0]
    
    # Angle between the projected bond vector and the x-axis
    angle = np.arctan2(bond_vector[1], bond_vector[0])
    
    #  Since we want the bond vector to align with the x-axis,
    #  we rotate by -angle
    R = rotation_matrix_z(-angle)
    
    # Apply rotation to all coordinates
    # transpose R for right-hand multiplication
    rotated_coords = np.dot(coords, R.T)  
    
    return rotated_coords

def round_and_clean_up(coords, tol=1e-15, decimals = 6):
    result = np.zeros_like(coords)
    for i, atom_vector in enumerate(coords):
        for j, component in enumerate(atom_vector):
            if np.abs(component) <= tol:
                result[i,j] = 0. 
            else:
                result[i,j] = np.round(component, decimals)
    return result



def align_xyz(filename):
    atom_col, coords = read_xyz(filename)
    coords -= coords[0]
    rotated_coords_z = align_bond_to_z_axis(coords)
    rotated_coords = align_bond_to_x_axis(rotated_coords_z)
    rotated_coords = round_and_clean_up(rotated_coords)
    
    name=filename.split('.')[0] 
    aligned_filename = name + '_aligned.xyz'

    write_xyz(aligned_filename, atom_col, rotated_coords)
    return aligned_filename

def aligned_to_com(aligned_filename, args):
    options = dict()
    filename_options = dict()
    get_density_functional(args, options)
    get_field_strength(args, options)
    get_solvent(args, options)
    get_basis_set(args, options)
    
    if args.qst3 == None:
        filename = aligned_filename
        options['filename'] = filename
        structure=get_coordinates(options, filename)
        options['coordinates'] = structure
        get_dotcom_filename(args, options, filename_options)
        get_dot_com(args, options, filename_options)
        
    elif args.qst3 != None:
        saddle_point=aligned_filename
        reactant=args.qst3[0]
        product=args.qst3[1]
        reactant_structure=get_coordinates(options, reactant)
        options['starting_material'] = reactant_structure
        product_structure=get_coordinates(options, product)
        options['product']=product_structure
        saddle_point_structure=get_coordinates(options, saddle_point)
        options['saddle_point'] = saddle_point_structure
        get_dotcom_filename(args, options, filename_options)
        get_dot_com(args, options, filename_options)

def main():
    args=get_arguments()
    filename=args.filename[0]
    #notation for getting one specific argument is args.argument/option
    if filename.__contains__(".log"):
        xyz_filename = log_to_xyz(args)
        aligned_filename=align_xyz(xyz_filename)
        aligned_to_com(aligned_filename, args)
    elif filename.__contains__(".xyz"):
        aligned_filename=align_xyz(filename)
        aligned_to_com(aligned_filename, args)
    else:
        print("invalid file type")
    
if __name__ == "__main__":
    main()
