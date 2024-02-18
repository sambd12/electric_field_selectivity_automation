#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 17:42:04 2024

@author: samantha.bealldenn
"""

import sys

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
    return zz_polar

filename=sys.argv[1]
z_dipole=get_z_dipole(filename)
zz_polarizability=get_z_polar(filename)
print("The z dipole moment is", z_dipole, 'atomic units.')
print("The zz polarizability is", zz_polarizability, 'atomic units.')



