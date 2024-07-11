#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 16:57:36 2024

@author: samantha.bealldenn
"""

import re
import argparse
from decomposing_energy import decompose_energy

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

def get_energy_of_last_structure(args, hybrid_status):
    #Find the lines that printed energies 
    #Grab only the energy
    #add all the energies as strings into a tuple  
    #grab the last energy
    #convert the energy to kJ/mol
    filename=args.filename[0]
    matches = []
    if hybrid_status == False:
        string_to_match = "SCF Done:"
        with open(filename) as f:
        	for line in f:
        		if re.search(string_to_match, line):
        			matches.append(line)			
        just_energies=()
        for line in matches:
            each_line=line.split()
            each_energy=float(each_line[4])
            just_energies= just_energies + (each_energy,)
        lowest_energy=just_energies[-1]
        kilojoules_energy=(lowest_energy*2625.5)
        print("\nInternal Energy:", kilojoules_energy, "kJ/mol")
    if hybrid_status == True:
        if filename.__contains__("b2plyp"):
            string_to_match = "B2PLYP"
            string_unmatched = "RB2PLYP"
        if filename.__contains__("b2plypd3"):
            string_to_match = "B2PLYPD3"
            string_unmatched = "RB2PLYPD3"
        with open(filename) as f:
            for line in f:
            	if re.search(string_to_match, line) and string_unmatched not in line and "#n" not in line:
            			matches.append(line)
            just_energies=()
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
        lowest_energy=just_energies[-1]
        kilojoules_energy=(lowest_energy*2625.5)
        print("Internal Energy:", kilojoules_energy, "kJ/mol")
        
        
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
        if re.search(" The following ModRedundant input section has been read:", line):
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
            string_unmatched = "RB2PLYP"
        if filename.__contains__("b2plypd3"):
            string_to_match = "B2PLYPD3"
            string_unmatched = "RB2PLYPD3"
        with open(filename) as f:
            for line in f:
            	if re.search(string_to_match, line) and string_unmatched not in line and "#n" not in line:
            			matches.append(line)
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
       if each_line[0] == "5":
           each_line[0]= "B"
       if each_line[0] == "6":
           each_line[0]= "C"
       if each_line[0] == "8":
           each_line[0]= "O"
       if each_line[0] == "9":
           each_line[0]= "F"
       if each_line[0] == "13":
           each_line[0] = "Al"
       if each_line[0] == "17":
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
        
def get_low_frequencies(args):
    filename=args.filename[0]
    string_to_match="Low frequencies"
    matches = []
    with open(filename) as f:
    	for line in f:
    		if re.search(string_to_match, line):
    			matches.append(line)
    print(matches[0],matches[1])
        
        
def get_z_dipole(args):
    filename=args.filename[0]
    with open (filename, 'r') as f:
        file_string=f.read()
    split_file_by_line=file_string.split("\n")
    electric_dipole_line_index=split_file_by_line.index(" Electric dipole moment (dipole orientation):")
# find the chunk of dipoles in dipole orientation
    z_dipole_line_index=electric_dipole_line_index+6
    z_dipole_line=split_file_by_line[z_dipole_line_index]
# find the line containing only the z dipole
    z_dipole_line_split=z_dipole_line.split(" ")
# line contains z and values with different unit systems, divide that line so each value is in its own index
    z_dipole_scientific=z_dipole_line_split[13]
# grab the dipole moment in atomic units
    z_dipole_scientific_split=z_dipole_scientific.split("D+0")
# divide up the scientific notation into the value and the scalar
    z_dipole_decimal=float(z_dipole_scientific_split[0])
    z_dipole_scalar=float(z_dipole_scientific_split[1])
    z_dipole=z_dipole_decimal*10**(z_dipole_scalar)
# get the dipole moment out of scientific notation by multiplying the value by 10^ scalar
    print("The z dipole moment is", z_dipole, 'atomic units.')
    
def get_z_polar(args):
    filename=args.filename[0]
    with open (filename, 'r') as f:
        file_string=f.read()
    split_file_by_line=file_string.split("\n")
    electric_polar_line_index=split_file_by_line.index(" Dipole polarizability, Alpha (dipole orientation).")
# find the chunk of polarizabilities in dipole orientation
    zz_polarizability_line_index=electric_polar_line_index+11
    zz_polarizability_line=split_file_by_line[zz_polarizability_line_index]
# find the line containing only the zz polarizability in dipole orientation
    zz_polarizability_line_split=zz_polarizability_line.split(" ")
# line contains zz and values with different unit systems, divided it so each is in its own index
    zz_polar_scientific=zz_polarizability_line_split[12]
#grab the zz polarizability in atomic units
    zz_polar_scientific_split=zz_polar_scientific.split("D+0")
#divide up the scientific notation its the value and the scalar
    zz_polar_decimal=float(zz_polar_scientific_split[0])
    zz_polar_scalar=float(zz_polar_scientific_split[1])
    zz_polar=zz_polar_decimal*10**(zz_polar_scalar)
#get the dpiole moment out of scientific notation by multiplying the value by 10^ scalar
    print("The zz polarizability is", zz_polar, 'atomic units.')

def log_to_xyz(args):
    termination_status=get_termination_status(args)
    xyz_filename=get_xyz_filename(args, termination_status)
    hybrid_status=get_density_function(args)
    freq_status=get_freq_status(args)
    molecule_length=get_molecule_length(args)
    if termination_status == "normal":
        get_energy_of_last_structure(args, hybrid_status)
        last_coordinates=get_last_coordinates(args, molecule_length)
        turn_coordinates_to_file(last_coordinates, xyz_filename)
    elif termination_status == "error_9999" or termination_status == "error_l103":
        lowest_energy_index=get_index_of_lowest_energy(args, hybrid_status)
        lowest_energy_coordinates=get_coordinates_of_lowest_energy(args, lowest_energy_index, molecule_length)
        turn_coordinates_to_file(lowest_energy_coordinates, xyz_filename)
    else:
        print("Unknown error or file is still running")
    if termination_status == "normal" and freq_status == True:
          get_z_dipole(args)
          get_z_polar(args)
          decompose_energy(args)
          get_low_frequencies(args)
    return xyz_filename

def main():
    args=get_arguments()
    #notation for getting one specific argument is args.argument/option
    #determine how to proceed, get the appropriate name, and the density functional
    log_to_xyz(args)
          
if __name__ == "__main__":
    main()
    
