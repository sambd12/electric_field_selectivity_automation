#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 15:57:24 2024

@author: samantha.bealldenn
"""

import re
import argparse
import sys
import pandas as pd

sys.getdefaultencoding()

def get_arguments():
    parser=argparse.ArgumentParser()
    #Adds essential information as arguments
    #gets filename from command line
    parser.add_argument("filename", nargs="+", action='store')
    #gets density functional from command line
    parser.add_argument("-s", "--spreadsheet", nargs=1, action='store')
    args=parser.parse_args()
    return args

def get_energy_breakdown(filename):
    with open (filename, 'r') as f:
        file_string=f.read()
    split_file_by_line=file_string.split("\n")
    for line in split_file_by_line:
        if re.search("                      Q", line):
            end_of_energy_corrections = split_file_by_line.index(line)
    ## This is the line where the energies change to Q, so we want to document the index to be able to collect the entire group of energies
    for line in split_file_by_line:
        if re.search(" Zero-point correction=", line):
            beginning_of_energy_corrections = split_file_by_line.index(line)
    ## this is the line where the free energies start to be listed, so we want to document this index
    energy_breakdown_lines= split_file_by_line[beginning_of_energy_corrections:end_of_energy_corrections]
    ## grab all free energies that are printed
    return energy_breakdown_lines

def get_free_energies(energy_breakdown_list):
    energy_breakdown = "\n".join(energy_breakdown_list)
    free_energies = energy_breakdown.split("                     E (Thermal)")[0]
    free_energies_list = free_energies.split("\n")
    for line in free_energies_list:
        if re.search("Sum of electronic and thermal Free Energies=", line):
            free_energy_line=line.split()
            free_energy=float(free_energy_line[-1])
            free_energy_kJ=(free_energy*2625.5)
        elif re.search("Sum of electronic and thermal Enthalpies=", line):
            enthalpy_line=line.split()
            enthalpy=float(enthalpy_line[-1])
            enthalpy_kJ=(enthalpy*2625.5)
        elif re.search("Sum of electronic and thermal Energies=", line):
            electronic_thermal_line=line.split()
            electronic_thermal_energy=float(electronic_thermal_line[-1])
            electronic_thermal_kJ=(electronic_thermal_energy*2625.5)
        elif re.search("Sum of electronic and zero-point Energies=", line):
             electronic_zeropt_line=line.split()
             electronic_zeropt_energy=float(electronic_zeropt_line[-1])
             electronic_zeropt_kJ=(electronic_zeropt_energy*2625.5)
        elif re.search("Zero-point correction=", line):
             zeropt_line=line.split()
             zeropt_energy=float(zeropt_line[-2])
             zero_pt_kJ=(zeropt_energy*2625.5)
        
    minus_temp_delta_s = free_energy_kJ - enthalpy_kJ
    thermal_energy = electronic_thermal_kJ - electronic_zeropt_kJ
    electronic_energy = electronic_zeropt_kJ - zero_pt_kJ
    
    print("Zero-point correction:", zero_pt_kJ, "kJ/mol")
    print("Electronic Energy:", electronic_energy, "kJ/mol")
    print("Thermal Energy:", thermal_energy, "kJ/mol")
    print("-T\u0394S=", minus_temp_delta_s, "kJ/mol")
    print("Free Energy:", free_energy_kJ, "kJ/mol")
    all_free_energies = [zero_pt_kJ, electronic_energy, thermal_energy, minus_temp_delta_s, free_energy_kJ]
    return all_free_energies
    
    
def get_internal_rotation_corrections(energy_breakdown_list, free_energy_kJ):
    energy_breakdown = "\n".join(energy_breakdown_list)
    internal_rotation_correction = energy_breakdown.split("                     E (Thermal)")[1]
    internal_rotation_correction_list = internal_rotation_correction.split("\n")
    for line in internal_rotation_correction_list:
        if re.search("internal rot", line):
            total_corrected_line=line.split()
            total_corrected=float(total_corrected_line[2])
        elif re.search("Total", line):
            total_thermal_line=line.split()
            total_thermal=float(total_thermal_line[1])
    ## gets only the correction for rotation
    internal_rot_correction = (total_corrected - total_thermal)
    ## converts to kJ
    internal_rot_correction_kJ = (internal_rot_correction*4.184) 
    ## adds rotational correction to the free energy
    total_free_energy= internal_rot_correction_kJ + free_energy_kJ
    print("Free Energy corrected by internal rotation:", total_free_energy, "kJ/mol")
    print("Correction by internal rotation:", internal_rot_correction_kJ, "kJ/mol")
    internal_rotation_corrections = [total_free_energy, internal_rot_correction_kJ]
    return internal_rotation_corrections

def parse_filename_for_info(filename):
    filename_info = []
    filename_split = filename.split("_")
    density_functional = filename_split[3]
    basis_set = filename_split[4]
    solvent = filename_split[5]
    reaction_pathway = filename_split[6]
    if len(filename_split) == 11:
        field_strength = filename_split[-3]
    elif len(filename_split) == 10:
        field_strength = filename_split[-2]
        
    if filename.__contains__('qst3') or filename.__contains__('ts'):
        structure_type = 'transition state'
        if len(filename_split) == 12:
            field_strength = filename_split[-3]
        elif len(filename_split) == 11:
            field_strength = filename_split[-2]
    elif filename.__contains__('reactant'):
        structure_type = 'reactant'
        reaction_pathway = "N/A"
    elif filename.__contains__('product'):
        structure_type = "product"
        reactant_conformer = "N/A"
    elif filename.__contains__('epoxide'):
        structure_type ='epoxide'
        reactant_conformer = "N/A"

    if filename.__contains__("longarm"):
        reactant_conformer = "long arm"
    elif filename.__contains__("shortarm"):
        reactant_conformer = "short arm"
    filename_info = [density_functional, basis_set, solvent, reactant_conformer, structure_type, reaction_pathway, field_strength]
    return filename_info

def write_energies_to_csv(tuple_of_energies, args):
    csv_filename = args.spreadsheet[0]
    df = pd.DataFrame(tuple_of_energies, 
                 columns=['Density Functional', 'Basis Set', 'Solvent', 'Reactant Conformer', 'Structure Type', 'Reaction Pathway', 'Field Strength', 'Zero-point Correction', 'Electronic Energy', 'Thermal Energy', 'minusT Delta S', 'Free Energy', 'Free Energy Corrected by Int. Rot.', 'Correction by Int. Rot.' ])
    df.to_csv(csv_filename, index=False)
    
def decompose_energy(args):
    tuple_of_energies = []
    ## each file and its energies will go into this list 
    for f in range(len(args.filename)):
        filename = args.filename[f]
        print("\nEnergies for", filename)
        energy_breakdown_list=get_energy_breakdown(filename)
        if filename.__contains__("_hr"):
            all_free_energies=get_free_energies(energy_breakdown_list)
            free_energy_kJ = all_free_energies[-1]
            internal_rotation_corrections=get_internal_rotation_corrections(energy_breakdown_list, free_energy_kJ)
        elif filename.__contains__("_freq"):
            all_free_energies=get_free_energies(energy_breakdown_list)
            internal_rotation_corrections=["N/A", "N/A"]
        filename_info=parse_filename_for_info(filename)
        all_energies_by_filename = filename_info + all_free_energies + internal_rotation_corrections
        tuple_of_energies.append(all_energies_by_filename)
    return tuple_of_energies
        
def main():
    args=get_arguments()
    tuple_of_energies=decompose_energy(args)
    if args.spreadsheet != None:
        write_energies_to_csv(tuple_of_energies, args)

if __name__ == "__main__":
    main()
    