#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 16:57:36 2024

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
            elif re.search(" Error termination request processed by link l103", line):
                print(line)

def get_filename(args, termination_status):
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


def main():
    args=get_arguments()
    #notation for getting one specific argument is args.argument/option
    #determine how to proceed, get the appropriate name, and the density functional
    
    termination_status=get_termination_status(args)
    xyz_filename=get_filename(args, termination_status)
    string_to_match=get_density_function(args)
    
    #Functions depending on how the script proceeeds
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
        
if __name__ == "__main__":
    main()
    