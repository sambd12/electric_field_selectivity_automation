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
from corrected_free_energy import get_entropy_corrected_G
from get_tunneling_info import get_tunneling_information
from get_internal_energy import get_energy_of_last_structure

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

def get_z_dipole(filename):
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
    return z_dipole
    
def get_z_polar(filename):
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
    return zz_polar

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
    
    print("Free Energy:", free_energy_kJ, "kJ/mol")
    print("Zero-point correction:", zero_pt_kJ, "kJ/mol")
    print("Electronic Energy:", electronic_energy, "kJ/mol")
    print("Thermal Energy:", thermal_energy, "kJ/mol")
    print("-T\u0394S=", minus_temp_delta_s, "kJ/mol")
    
    free_E_breakdown = [zero_pt_kJ, electronic_energy, thermal_energy, minus_temp_delta_s]
    return free_energy_kJ, free_E_breakdown
    
    
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

def get_low_frequencies(filename):
    string_to_match="Low frequencies"
    matches = []
    with open(filename) as f:
    	for line in f:
    		if re.search(string_to_match, line) and "Full mass-weighted force constant matrix" not in line:
    			matches.append(line)
    print(matches[0],matches[1])
    split_frequencies=matches[0].split()
    first_frequency=abs(float(split_frequencies[3]))
    return first_frequency

def parse_filename_for_info(filename):
    filename_info = []
    filename_split = filename.split("_")
    
    structure_type = filename_split[2]
    density_functional = filename_split[3]
    basis_set = filename_split[4]
    solvent = filename_split[5]
    reaction_pathway = filename_split[6]

    if len(filename_split) <= 10:
        field_strength = filename_split[8]
    else:
        field_strength = filename_split[9]
        

    if filename.__contains__("longarm1"):
        reactant_conformer = "long arm1"
    elif filename.__contains__("longarm2"):
        reactant_conformer = "long arm2"
    elif filename.__contains__("longarm"):
        reactant_conformer = "long arm"
    elif filename.__contains__("shortarm1"):
        reactant_conformer = "short arm1"
    elif filename.__contains__("shortarm2"):
        reactant_conformer = "short arm2"
    elif filename.__contains__("shortarm"):
        reactant_conformer = "short arm"
        
    if field_strength == "nofield":
        field_integer = 0
    elif field_strength.__contains__("neg"):
        field_direction = "-"
        field_strength = field_strength.split("g")[1]
        field_integer = int(field_direction + field_strength)
    elif field_strength.__contains__("pos"):
        field_strength = field_strength.split("s")[1]
        field_integer = int(field_strength)
        
    filename_info = [density_functional, basis_set, solvent, reactant_conformer, structure_type, reaction_pathway, field_integer]
    return filename_info

def write_energies_to_csv(info_by_file, args):
    csv_filename = args.spreadsheet[0]
    filename = args.filename[0]
    
    files_by_field_strength = sorted(info_by_file, reverse=True)
    sorted_files = [info_by_file[key] for key in files_by_field_strength]
    
    if filename.__contains__("_hr"):
        df = pd.DataFrame(sorted_files, 
                 columns=['Density Functional', 'Basis Set', 'Solvent', 'Reactant Conformer', 'Structure Type', 'Reaction Pathway', 'Field Strength','Internal Energy','Free Energy','Free Energy Corrected by Int. Rot.', 'Correction by Int. Rot.','Quasi-Rho G','Quasi-Harmonic G','Zero-point correction','Electronic Energy','Thermal Energy','minusT Delta S', 'z dipole','zz polar', 'Standard Wigner Tunneling Coeff', "Truncated Wigner Tunneling Coeff" ])
    elif filename.__contains__("_freq"):
        df = pd.DataFrame(sorted_files, 
                 columns=['Density Functional', 'Basis Set', 'Solvent', 'Reactant Conformer', 'Structure Type', 'Reaction Pathway', 'Field Strength','Internal Energy','Free Energy','Quasi-Rho G','Quasi-Harmonic G','Zero-point correction','Electronic Energy','Thermal Energy','minusT Delta S','z dipole','zz polar','Standard Wigner Tunneling Coeff', "Truncated Wigner Tunneling Coeff"])
    df.to_csv(csv_filename, index=False)
    
def decompose_energy(args):
    info_by_file = dict()
    ## each file and its energies will go into this list 
    for f in range(len(args.filename)):
        filename = args.filename[f]
        print("\nEnergies for", filename)
    
        energy_breakdown_list=get_energy_breakdown(filename)

        filename_info=parse_filename_for_info(filename)
        field_strength=filename_info[-1]
        
        if filename.__contains__("mn15"):
           hybrid_status = False
        elif filename.__contains__("b3lyp"):
           hybrid_status = False
        elif filename.__contains__("b2plyp"):
            hybrid_status = True
            
        internal_energy = get_energy_of_last_structure(filename, hybrid_status)
        free_energy_kJ, free_E_breakdown = get_free_energies(energy_breakdown_list)
        
        energies = [internal_energy, free_energy_kJ]
        
        
        entropy_corrected_G = get_entropy_corrected_G(filename, temperature=None, w0=100.)
        
        dipole = get_z_dipole(filename)
        polarizability = get_z_polar(filename)
        properties = [dipole, polarizability]
        
        if filename.__contains__("_hr"):
            
            internal_rotation_corrections = get_internal_rotation_corrections(energy_breakdown_list, free_energy_kJ)
            all_energies_by_filename = filename_info + energies + internal_rotation_corrections +  entropy_corrected_G + free_E_breakdown
        
        elif filename.__contains__("_freq"):
            
            all_energies_by_filename = filename_info + energies + entropy_corrected_G + free_E_breakdown
        
        first_freq=get_low_frequencies(filename)
        tunneling_info=get_tunneling_information(first_freq)
        wigner_coeffs=tunneling_info[1:]
        all_info_by_filename = all_energies_by_filename + properties + wigner_coeffs
        info_by_file[field_strength] = all_info_by_filename
 
    return info_by_file
        
def main():
    args=get_arguments()
    dictionary_of_info=decompose_energy(args)
    if args.spreadsheet != None:
        write_energies_to_csv(dictionary_of_info, args)

if __name__ == "__main__":
    main()
    