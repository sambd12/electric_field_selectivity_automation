#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 16:57:36 2024

@author: samantha.bealldenn
"""

import re
import argparse
from decomposing_energy import decompose_energy, get_z_dipole, get_z_polar
from get_internal_energy import get_energy_of_last_structure
from openbabel import pybel


def get_arguments():
    parser=argparse.ArgumentParser()
    #Adds essential information as arguments
    #gets filename from command line
    parser.add_argument("filename", nargs=1, action='store')
    #gets density functional from command line
    args=parser.parse_args()
    return args

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
            elif re.search("l103", line):
                termination_status="error_l103"
                return(termination_status)

def get_xyz_filename(args, termination_status):
    log_filename=args.filename[0]
    name=log_filename.split('.')[0] 
    if termination_status == "normal":
        xyz_filename = name + '_structure.xyz'
    elif termination_status == "error_9999":
        xyz_filename = name + '_9999structure.xyz'
    elif termination_status == "error_l103":
        xyz_filename = name + '_103structure.xyz'
    return xyz_filename

def get_density_function(args):
    filename=args.filename[0]
    if filename.__contains__("mn15"):
       hybrid_status = False
    elif filename.__contains__("b3lyp"):
       hybrid_status = False
    elif filename.__contains__("b2plyp"):
        hybrid_status = True
    else:
       print("Filename must contain a valid density functional")
    return hybrid_status

def get_freq_status(args):
    filename=args.filename[0]
    if filename.__contains__("freq") or filename.__contains__("_hr"):
        freq_status = True
    else:
        freq_status = False
    return freq_status

#The following two functions are only for normal termination
        
        
def get_molecule_length(args):
    filename=args.filename[0]
    with open (filename, 'r') as f:
        file_string=f.read()
    split_file_by_line=file_string.split("\n")
    for line in split_file_by_line:
        if re.search(" Symbolic Z-matrix:", line):
            beginning_of_molecule_index = split_file_by_line.index(line)
            break
    for line in split_file_by_line:
        if re.search(" Add virtual bond connecting atoms", line):
            end_of_molecule_index = split_file_by_line.index(line) - 1 
            break
        elif re.search(" The following ModRedundant input section has been read:", line):
            end_of_molecule_index = split_file_by_line.index(line) - 1 
            break
        elif re.search("Variables:", line):
            end_of_molecule_index = split_file_by_line.index(line)
            break
    ## this is the line where the free energies start to be listed, so we want to document this index
    input_molecule= split_file_by_line[beginning_of_molecule_index:end_of_molecule_index]
    
    molecule_length=len(input_molecule)

    ## grab all free energies that are printed
    return molecule_length

def get_last_coordinates(args, molecule_length):
    filename=args.filename[0]
    sets_of_coords = []
# gets the last set of coordinates
# includes numbers and spaces that need to be taken out
    with open(filename, 'r') as f:
        file_string=f.read()      
    split_file_by_paragraph=file_string.split("---------------------------------------------------------------------")
## each set of coordinates is surrounded on either side by -------...
## splits the entire file by these and looks at what is in between each
    for paragraph in split_file_by_paragraph:
          paragraph_by_line=paragraph.split("\n")
## splits each chunk into lines
          if len(paragraph_by_line) == molecule_length and "-" not in paragraph_by_line[0]:
##ensures that the group of coordinates is the correct length (# of atoms +2)
## and that we grabbed a set of coordinates, not another analysis that happened to be that long
## if we have extra - separating our sections, they are not the correct format of coordinates
              coordinates="\n".join(paragraph_by_line)
##rejoins the lines
              sets_of_coords.append(coordinates)
##adds the set of coordinates to the possible matches
    last_coordinates=sets_of_coords[-1]
##finds the last match to get the last coordinates
    return last_coordinates

#The following two functions are for 9999 errors only

def get_index_of_lowest_energy(args, hybrid_status):
#determine which energy is lowest
#determine the index of the lowest energy
#convert the lowest energy to kilojoules per mole
    filename=args.filename[0]
    matches = []
    if hybrid_status == False:
        string_to_match = "SCF Done:"
        with open(filename) as f:
        	for line in f:
        		if re.search(string_to_match, line):
        			matches.append(line)       
        if len(matches) <= 5:
            print("Warning! Gaussian completed less than 5 cycles before erroring!")
        just_energies=()
        for line in matches:
            each_line=line.split()
            each_energy=float(each_line[4])
            just_energies= just_energies + (each_energy,)
        lowest_energy=min(just_energies)
        lowest_energy_index=just_energies.index(lowest_energy)
        lowest_energy_kJ=(lowest_energy*2625.5)
        print("Internal Energy:", lowest_energy_kJ, "kJ/mol")
        return lowest_energy_index
    if hybrid_status == True:
        if filename.__contains__("b2plyp"):
            string_to_match = "B2PLYP"
        if filename.__contains__("b2plypd3"):
            string_to_match = "B2PLYPD3"
        with open(filename) as f:
            for line in f:
            	if re.search(string_to_match, line) and re.search("E2", line):
            			matches.append(line)
        if len(matches) <= 5:
            print("Warning! Gaussian completed less than 5 cycles before erroring!")
        just_energies=()
        ## first line comes from the input of the file, need to delete that one
        for line in matches:
                each_line=line.split()
                hartrees_scientific_energy=each_line[-1]
                hartrees_scientific_split=hartrees_scientific_energy.split("D+0")
            # divide up the scientific notation into the value and the scalar
                hartrees_decimal=float(hartrees_scientific_split[0])
                hartrees_scalar=float(hartrees_scientific_split[1])
                hartrees=hartrees_decimal*10**(hartrees_scalar)
            # multiply and add to the list of energies
                just_energies= just_energies + (hartrees,)
        lowest_energy=min(just_energies)
            # get the lowest of the energies
        lowest_energy_index=just_energies.index(lowest_energy)
        lowest_energy_kJ=(lowest_energy*2625.5)
        print("Internal Energy:", lowest_energy_kJ, "kJ/mol")
        return lowest_energy_index

def get_coordinates_of_lowest_energy(args, lowest_energy_index, molecule_length):
    filename=args.filename[0]
    sets_of_coords=()
    with open(filename, 'r') as f:
        file_string=f.read()
    split_file_by_paragraph=file_string.split("---------------------------------------------------------------------")
    for paragraph in split_file_by_paragraph:
        paragraph_by_line=paragraph.split("\n")
        if len(paragraph_by_line) == molecule_length and "-" not in paragraph_by_line[0]:
#ensures that the group of coordinates is the correct length (# of atoms +2), which means that it is the correct output format
## and that we grabbed a set of coordinates, not another analysis that happened to be that long
## if we have extra - separating our sections, they are not the correct format of coordinates
            coordinates="\n".join(paragraph_by_line)
            sets_of_coords=sets_of_coords+(coordinates,)
    lowest_energy_coordinates= sets_of_coords[lowest_energy_index]    
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
       elif each_line[0] == "5":
           each_line[0]= "B"
       elif each_line[0] == "6":
           each_line[0]= "C"
       elif each_line[0] == "7":
            each_line[0] = "N"
       elif each_line[0] == "8":
           each_line[0]= "O"
       elif each_line[0] == "9":
           each_line[0]= "F"
       elif each_line[0] == "13":
           each_line[0] = "Al"
       elif each_line[0] == "15":
            each_line[0] = "P"
       elif each_line[0] == "17":
           each_line[0] = "Cl"
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
    print("Structure coordinates:", xyz_filename)
    with open(xyz_filename, 'w') as f:
        f.write(file_xyzstring)
    
    
def change_file_format(args, zmatrix):
    filename=args.filename[0]
    name=args.filename[0].split('.')[0] 
    xyz_filename = name + '_structurefromcom.xyz'
    
    zmatrix_syntax=['!Put Keywords Here, check Charge and Multiplicity.', "#", "", "Comment goes here", "" ]
    
    with open (filename, 'r') as f:
        file_string=f.read()
    split_file_by_line = file_string.split("\n")
    del split_file_by_line[0:7]
    del split_file_by_line[-1]
    
    
    if "qst3" not in filename and zmatrix == False:
        if zmatrix == False:
             del split_file_by_line[-1]
             atoms=[]
             for line in split_file_by_line:
                 each_line=line.split()
                 if len(each_line) == 4:
                     atoms.append(line)
             num_atoms = (len(atoms))
             xyz_syntax = [str(num_atoms), '']
             atoms = xyz_syntax + atoms
             
        elif zmatrix == True:
             split_file_by_line = zmatrix_syntax + split_file_by_line
             file_string="\n".join(split_file_by_line)
             molecule = pybel.readstring("gzmat", file_string)
             molecule_xyz = molecule.write(format="xyz")
             atoms = molecule_xyz.split("\n")
             del atoms[-1]
             
        molecule_xyz = "\n".join(atoms)
        file_xyzstring=str(molecule_xyz)
        print("Structure coordinates:", xyz_filename)
        with open(xyz_filename, 'w') as f:
            f.write(file_xyzstring)
         
    elif filename.__contains__("qst3"):
        if zmatrix == False:
            del split_file_by_line[-1]
            atoms=[]
            for line in split_file_by_line:
                each_line=line.split()
                if len(each_line) == 4:
                    atoms.append(line)
            each_atom = int(len(atoms)/3)
            xyz_syntax = [str(each_atom), '']
            reactant = atoms[:each_atom]
            product = atoms[each_atom:each_atom*2]
            ts = atoms[each_atom*2:]
         
        elif zmatrix == True:
            for line in split_file_by_line:
                if line.__contains__("Product"):
                    product_start = split_file_by_line.index(line)
                if line.__contains__("Saddle Point Guess"):
                    TS_start = split_file_by_line.index(line)
                    
            reactant=['0 1'] + split_file_by_line[:product_start-1]
            product=split_file_by_line[product_start+2:TS_start-1]
            ts=split_file_by_line[TS_start+2:-1]
        
        #put molecules into a list
        molecule_list = [reactant, product, ts]
        name = filename.split('.')[0]
        name = name.split('_')
        
        #change molecule names
        
        reactant_name = name[:2] + ["reactant"] + name[3:6] + name[7:]
        reactant_name = "_".join(reactant_name) + '_structurefromcom.xyz'
        
        product_name = name[:2] + ["product"] + name[3:8] + name[9:]
        product_name = "_".join(product_name) + '_structurefromcom.xyz'

        qst3_name = "_".join(name) + '_structurefromcom.xyz'

        molecule_names = [reactant_name, product_name, qst3_name]
        for i, mol in enumerate(molecule_list):
            if zmatrix == False:
                mol = xyz_syntax + mol
                molecule_xyz="\n".join(mol)
            elif zmatrix == True:
                mol = zmatrix_syntax + mol
                mol_string="\n".join(mol)
                molecule = pybel.readstring("gzmat", mol_string)
                molecule_xyz = molecule.write(format="xyz")
                molecule_xyz = molecule_xyz.split("\n")
                del molecule_xyz[-1]
                molecule_xyz = "\n".join(molecule_xyz)
            file_xyzstring=str(molecule_xyz)
            print("Structure coordinates:", molecule_names[i])
            with open(molecule_names[i], 'w') as f:
                f.write(file_xyzstring)
        

def log_to_xyz(args):
    filename=args.filename[0]
    termination_status=get_termination_status(args)
    xyz_filename=get_xyz_filename(args, termination_status)
    hybrid_status=get_density_function(args)
    freq_status=get_freq_status(args)
    molecule_length=get_molecule_length(args)
    
    if termination_status == "normal":
        get_energy_of_last_structure(filename, hybrid_status)
        last_coordinates=get_last_coordinates(args, molecule_length)
        turn_coordinates_to_file(last_coordinates, xyz_filename)
    elif termination_status == "error_9999" or termination_status == "error_l103":
        lowest_energy_index=get_index_of_lowest_energy(args, hybrid_status)
        lowest_energy_coordinates=get_coordinates_of_lowest_energy(args, lowest_energy_index, molecule_length)
        turn_coordinates_to_file(lowest_energy_coordinates, xyz_filename)
    else:
        print("Unknown error or file is still running")
        
    if termination_status == "normal" and freq_status == True:
          decompose_energy(args)
          
    return xyz_filename
     

def com_to_xyz(args):
    filename=args.filename[0]
    
    with open (filename, 'r') as f:
        file_string=f.read()
    split_file_by_line = file_string.split("\n")
    
    zmatrix = False
    for line in split_file_by_line:
        if re.search("Variables:", line):
            zmatrix = True
            break
        
    change_file_format(args, zmatrix)
        
    
def main():
    args=get_arguments()
    filename=args.filename[0]
    #notation for getting one specific argument is args.argument/option
    #determine how to proceed, get the appropriate name, and the density functional
    if filename.__contains__(".log"):
        log_to_xyz(args)
    elif filename.__contains__(".com"):
        com_to_xyz(args)
    else:
        print("invalid filetype. Please submit a .log or .com")
          
if __name__ == "__main__":
    main()
    
