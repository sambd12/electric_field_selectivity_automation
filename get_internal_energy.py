#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 15 11:28:59 2025

@author: samantha.bealldenn
"""


import re

def get_energy_of_last_structure(filename, hybrid_status):
    #Find the lines that printed energies 
    #Grab only the energy
    #add all the energies as strings into a tuple  
    #grab the last energy
    #convert the energy to kJ/mol
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
    if hybrid_status == True:
        if filename.__contains__("b2plyp"):
            string_to_match = "B2PLYP"
        if filename.__contains__("b2plypd3"):
            string_to_match = "B2PLYPD3"
        with open(filename) as f:
            for line in f:
            	if re.search(string_to_match, line) and re.search("E2", line):
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
    return kilojoules_energy