#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 16:16:53 2024

@author: samantha.bealldenn
"""
import re
import numpy as np
import argparse
from openbabel import pybel

def get_arguments():
    parser=argparse.ArgumentParser()
    #Adds essential information as arguments
    #gets filename from command line
    parser.add_argument("filename", nargs=1, action='store')
    #gets density functional from command line
    parser.add_argument("density_functional", choices=["b3lyp", "mn15"], nargs=1, action="store")
    #gets field stength from command line
    parser.add_argument("field_strength", choices=["n4", "n3", "n2", "n1", 'nofield', 'p1', 'p2', 'p3', 'p4'], nargs=1, action="store")
    #gets solvent from command line
    parser.add_argument("solvent", choices=["acn", "dcm"], nargs=1, action='store')
    
    group = parser.add_mutually_exclusive_group(required=False)
    #changes the regular minimization to a single point calculation
    group.add_argument("-s", "--spc", action="store_true")
    #adds the option to make the output file QST3, in which case we want to get the names of the other two files being used
    group.add_argument("-q", "--qst3", nargs=2, action='store')

    args=parser.parse_args()
    return args
    #args is the entire namespace

#density functional options
def get_density_functional(args, options):
    if args.density_functional.__contains__('b3lyp'):
        options['density_functional'] = "B3LYP/6-311+g(d,p) EmpiricalDispersion=GD3"
    elif args.density_functional.__contains__('mn15'):
        options['density_functional'] = "MN15/6-311+g(d,p)"
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

def get_termination_status(args):
    filename=args.filename[0]
    with open(filename) as f:
        for line in f:
            if re.search(" Normal termination of Gaussian 16 at", line):
                termination_status="normal"
                return(termination_status)
            elif re.search(" Error termination request processed by link 9999", line):
                termination_status="error_9999"
                return(termination_status)
            elif re.search(" Error termination request processed by link l103", line):
                print(line)

def get_xyz_filename(args, termination_status):
    log_filename=args.filename[0]
    name=log_filename.split('.')[0] 
    if termination_status == "normal":
        xyz_filename = name + '_structure.xyz'
    elif termination_status == "error_9999":
        xyz_filename = name + '_9999structure.xyz'
    return xyz_filename

def get_density_function(args):
    density_functional=args.density_functional[0]
    if density_functional.__contains__("mn15"):
       string_to_match="RMN15"
    elif density_functional.__contains__("b3lyp"):
       string_to_match="RB3LYP"
    else:
       print("enter a valid density function")
    return string_to_match

#The following two functions are only for normal termination

def get_energy_of_last_structure(args, string_to_match):
    filename=args.filename[0]
    matches = []
    with open(filename) as f:
    	for line in f:
    		if re.search(string_to_match, line):
    			matches.append(line)			
#match the density functional I found to the lines that contain the structure energies 
    just_energies=()
    for line in matches:
            each_line=line.split()
            if len(each_line) > 4:
#testing if the line contains energies, not just the matching string
               each_energy=float(each_line[4])
#Grab only the energy from the lines that contain the density functional match
               just_energies= just_energies + (each_energy,)
#add all the energies as strings into a tuple
            else:
               continue     
    lowest_energy=just_energies[-1]
#grab the last energy
    kilojoules_energy=(lowest_energy*2625.5)
#convert the energy to kJ/mol
    print(kilojoules_energy, "kJ/mol")

def get_last_coordinates(args):
    filename=args.filename[0]
# gets the last set of coordinates
# includes numbers and spaces that need to be taken out
# is a tuple, each line/atom is a string in it
    with open(filename, 'r') as f:
        file_string=f.read()      
    split_file=file_string.split("---------------------------------------------------------------------")
    last_coordinates=split_file[-2]
    return last_coordinates

#The following two functions are for 9999 errors only

def get_index_of_lowest_energy(args, string_to_match):
    filename=args.filename[0]
    matches = []
    with open(filename) as f:
    	for line in f:
    		if re.search(string_to_match, line):
    			matches.append(line)
#match the density functional I found to the lines that contain the structure energies 
    just_energies=()
    for line in matches:
            each_line=line.split()
            if len(each_line) > 4:
#testing if the line contains energies, not just the matching string
               each_energy=float(each_line[4])
#grab only the energy
               just_energies= just_energies + (each_energy,)
#add all energies as string to a tuple
            else:
               continue 
    lowest_energy=min(just_energies)
#determine which energy is lowest
    lowest_energy_index=just_energies.index(lowest_energy)
#determine the index of the lowest energy
    lowest_energy_kJ=(lowest_energy*2625.5)
#convert the lowest energy to kilojoules per mole
    print(lowest_energy_kJ, "kJ/mol")
    return lowest_energy_index

def get_coordinates_of_lowest_energy(args, lowest_energy_index):
    filename=args.filename[0]
    sets_of_coords=()
    with open(filename, 'r') as f:
        file_string=f.read()
    split_file_by_paragraph=file_string.split("---------------------------------------------------------------------")
    for paragraph in split_file_by_paragraph:
        paragraph_by_line=paragraph.split("\n")
        if len(paragraph_by_line) == (33):
#ensures that the group of coordinates is the correct length (# of atoms +2), which means that it is the correct output format
            coordinates="\n".join(paragraph_by_line)
            sets_of_coords=sets_of_coords+(coordinates,)
#takes all the sets of coordinates that are in the correct format
    if (lowest_energy_index*3-2) > len(sets_of_coords):
        lowest_energy_coordinates= sets_of_coords[lowest_energy_index]
# checks to see if the sets_of_coords grabbed only the correctly formatted sets of coordinates
# sometimes it grabs extra coordinates that are not correct, which the formula index*3-2 corrects for
    elif (lowest_energy_index*3-2) < len(sets_of_coords):
        lowest_energy_coordinates=sets_of_coords[(lowest_energy_index*3-2)]
#weird formula that must be used as there are 3 sets of coordinates in the correct format that will be added per iteration done by Gaussian
#we want the first one in the set of 3, so we multiply the index by 3 to get to the correctly group of coordinate sets
#we subtract by 2 because the 1st set of coords in that group of 3 is in the format we want
    return lowest_energy_coordinates

#The following two functions are for both 9999 errors and normal termination 

def turn_coordinates_to_file(coordinates, xyz_filename):
#removes extraneous info from each line
    coordinate_lines=coordinates.split("\n") 
    del coordinate_lines[0]
    coordinate_lines.pop()
# first and last entry in the coordinates chunk are unneccessary
    clean_file=()
    for line in coordinate_lines:
       each_line=line.split()
# take out blank spaces & index of each atom in the structure
       del each_line[0]
       del each_line[1]
#change atomic numbers to their letters
       if each_line[0] == "1":
           each_line[0]= "H"
       if each_line[0] == "5":
           each_line[0]= "B"
       if each_line[0] == "6":
           each_line[0]= "C"
       if each_line[0] == "8":
           each_line[0]= "O"
       if each_line[0] == "9":
           each_line[0]= "F"
#put the atoms back together, adding each atom line as a string to a tuple
       each_line=" ".join(each_line)
       clean_file= clean_file + (each_line,)
#put the # of atoms at the top of the file & add that to the cleaned up coordinates
    vmd_syntax = (str(len(clean_file)), "",)
    vmd_file=(vmd_syntax + clean_file)
    complete_file="\n".join(vmd_file)
#joins all the tuples to be one string
    file_xyzstring=str(complete_file)
#writes file as .xyz
    with open(xyz_filename, 'w') as f:
        f.write(file_xyzstring)

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
        filename_options['molecule_type'] = "product"
    elif filename.__contains__("reactant"):
        filename_options['molecule_type'] = "reactant"
    else:   
        filename_options['molecule_type'] = "qst3"
   
    ## adds aldehyde and ketone to the name of the file
    ## in the case of QST3 files, adds aldehyde and ketone to the file    
    if filename.__contains__("ketone"):
        filename_options['pathway'] = "ketone"
        options['product_type'] = "Ketone"
    elif filename.__contains__("aldehyde"):
        filename_options['pathway'] = "aldehyde"
        options['product_type'] = "Aldehyde"
    else:
        filename_options['pathway'] = ""
        
    if filename.__contains__("_s_"):
        filename_options['enantiomer'] = "s"
    elif filename.__contains__("_r_"):
        filename_options['enantiomer'] = "r"
    else:
        filename_options['enantiomer'] = ""
        
    if filename.__contains__("_longarm_"):
        filename_options['reactant_type'] = "longarm"
    elif filename.__contains__("_shortarm_"):
        filename_options['reactant_type'] = "shortarm"
    else:
        filename_options['reactant_type'] = ""

#gets the density functional, solvent, and field strength options from the command line
    if args.density_functional == ['b3lyp']:
        filename_options['density_functional'] = 'b3lyp'
    elif args.density_functional == ['mn15']:
        filename_options['density_functional'] = 'mn15'
    
    if args.solvent == ['acn']:
        filename_options['solvent'] = 'acn'
    elif args.solvent == ['dcm']:
        filename_options['solvent'] = 'dcm'
    
    if args.field_strength == ['n4']:
        filename_options['field_strength'] = "neg4"
    elif args.field_strength == ['n3']:
        filename_options['field_strength'] = "neg3"
    elif args.field_strength == ['n2']:
        filename_options['field_strength'] = "neg2"
    elif args.field_strength == ['n1']:
        filename_options['field_strength'] = "neg1"
    elif args.field_strength == ['nofield']:
        filename_options['field_strength'] = "nofield"
    elif args.field_strength == ['p1']:
        filename_options['field_strength'] = "pos1"
    elif args.field_strength == ['p2']:
        filename_options['field_strength'] = "pos2"
    elif args.field_strength == ['p3']:
        filename_options['field_strength'] = "pos3"
    elif args.field_strength == ['p4']:
        filename_options['field_strength'] = "pos4"
    return filename_options

def get_dot_com(args, options, filename_options):
    #JUST reactant or product
    if args.spc == False and args.qst3 == None:
        syntax="%mem=24GB\n%NProcShared=32\n#n opt=Z-Matrix NoSymm {density_functional} {field_strength} {solvent}\n\n {filename}\n\n{coordinates}\n".format_map(options)
        completed_filename=("cis-stilbene_oxide_{molecule_type}_{density_functional}_{solvent}_{pathway}_{enantiomer}_{reactant_type}_{field_strength}.com").format_map(filename_options)
    elif args.spc == True:
        syntax="%mem=24GB\n%NProcShared=32\n#n NoSymm {density_functional} {field_strength} {solvent} Polar\n\n {filename}\n\n{coordinates}\n".format_map(options)
        completed_filename=("cis-stilbene_oxide_{molecule_type}_{density_functional}_{solvent}_{pathway}_{enantiomer}_{reactant_type}_{field_strength}_spc.com").format_map(filename_options)
    elif args.qst3 != None:
        syntax="%mem=24GB\n%NProcShared=32\n#n opt=(Z-Matrix,QST3) NoSymm {density_functional} {field_strength} {solvent}\n\nStarting Material\n\n{starting_material}\n{product_type} Product\n\n{product}\nSaddle Point Guess\n\n{saddle_point}\n".format_map(options)
        completed_filename=("cis-stilbene_oxide_qst3_{density_functional}_{solvent}_{pathway}_{enantiomer}_{reactant_type}_{field_strength}.com").format_map(filename_options)
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

def log_to_xyz(args):
    termination_status=get_termination_status(args)
    xyz_filename=get_xyz_filename(args, termination_status)
    string_to_match=get_density_function(args)
    if termination_status == "normal":
        get_energy_of_last_structure(args, string_to_match)
        last_coordinates=get_last_coordinates(args)
        turn_coordinates_to_file(last_coordinates, xyz_filename)
    elif termination_status == "error_9999":
        lowest_energy_index=get_index_of_lowest_energy(args, string_to_match)
        lowest_energy_coordinates=get_coordinates_of_lowest_energy(args, lowest_energy_index)
        turn_coordinates_to_file(lowest_energy_coordinates, xyz_filename)
    else:
        print("Unknown error or file is still running")
    return xyz_filename

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
    print(filename)
    #notation for getting one specific argument is args.argument/option
    if filename.__contains__(".log"):
        xyz_filename = log_to_xyz(args)
        aligned_filename=align_xyz(xyz_filename)
        aligned_to_com(aligned_filename, args)
    elif filename.__contains__(".xyz"):
        aligned_filename=align_xyz(filename)
        aligned_to_com(aligned_filename, args)
    else:
        print("invalid")
    
if __name__ == "__main__":
    main()