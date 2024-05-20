#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 17:42:04 2024

@author: samantha.bealldenn
"""

import re
import argparse

def get_arguments():
    parser=argparse.ArgumentParser()
    #Adds essential information as arguments
    #gets filename from command line
    parser.add_argument("filename", nargs=1, action='store')
    #gets density functional from command line
    parser.add_argument("density_functional", choices=["b3lyp", "mn15"], nargs=1, action="store")
    #gets field stength from command line
    parser.add_argument("-p", "--polar", action="store_true")
    args=parser.parse_args()
    return args

def get_xyz_filename(args):
    log_filename=args.filename[0]
    name=log_filename.split('.')[0] 
    xyz_filename = name + '_structure.xyz'
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
    if args.polar == True:
        last_coordinates=split_file[4]
    elif args.polar == False:
        last_coordinates=split_file[-4]
    return last_coordinates


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
    with open(xyz_filename, 'w') as f:
        f.write(file_xyzstring)

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
    return z_dipole
    
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
    return zz_polar


def main():
    args=get_arguments()
    #notation for getting one specific argument is args.argument/option
    # get the name, and the density functional
    
    xyz_filename=get_xyz_filename(args)
    string_to_match=get_density_function(args)
    get_energy_of_last_structure(args, string_to_match)
    last_coordinates=get_last_coordinates(args)
    turn_coordinates_to_file(last_coordinates, xyz_filename)
    if args.polar == True:
        z_dipole=get_z_dipole(args)
        zz_polarizability=get_z_polar(args)
        print("The z dipole moment is", z_dipole, 'atomic units.')
        print("The zz polarizability is", zz_polarizability, 'atomic units.')
        
        
if __name__ == "__main__":
    main()



