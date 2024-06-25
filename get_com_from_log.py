#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 16:16:53 2024

@author: samantha.bealldenn
"""

import argparse
from openbabel import pybel
from get_xyz_from_log import log_to_xyz
from align import align_xyz

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
    parser.add_argument("solvent", choices=["acn", "dcm", "nosolv"], nargs=1, action='store')
    parser.add_argument("basis_set", choices=["+", "++", "pvdz", "pvtz", "pvqz"], nargs=1, action="store")
    
    group = parser.add_mutually_exclusive_group(required=False)
    #changes the regular minimization to a single point calculation
    group.add_argument("-s", "--spc", action="store_true")
    #adds the option to make the output file QST3, in which case we want to get the names of the other two files being used
    group.add_argument("-q", "--qst3", nargs=2, action='store')
    
    group2=parser.add_mutually_exclusive_group(required=False)
    #adds the option to include frequency calculations, which factor entropy into the free energy calculations.
    # can't have just frequency and hindered rotor
    group2.add_argument("-f", "--freq", action="store_true")
    group2.add_argument("-h", "--hindered_rotor", action='store_true')
    
    parser.add_argument("-m", "--memory", nargs=2, action='store')

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
    elif args.solvent == ['nosolv']:
        options['solvent'] = ""
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
    elif args.basis_set == ['pvqz']:
            options['basis_set'] = "AUG-cc-pVQZ"
    return options

def get_frequency(args, options):
    if args.freq == True:
        options['frequency'] = "Freq"
        options['hindered_rotor'] = ""
    elif args.freq == False and args.hinderedrotor == False:
        options['frequency'] = ""
        options['hindered_rotor'] = ""
    elif args.hindered_rotor == True:
        options['frequency'] = "Freq=(HinderedRotor,ReadHinderedRotor)"
        options['hindered_rotor'] = "\n1.0\n1 2 3 1 1\n"
    return options

def get_memory_and_processors(args,options):
    if args.memory == None:
        options['memory'] = "24"
        options['processors'] = "32"
    if args.memory != None:
        options['memory'] = args.memory[0]
        options['processors'] = args.memory[1]
        
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
        options['product_type'] = "Other"
        
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
    elif args.solvent ==['nosolv']:
        filename_options['solvent'] = "_nosolv"
        
        
    if args.basis_set == ['+']:
        filename_options['basis_set'] = ""
    if args.basis_set == ['++']:
        filename_options['basis_set'] = "_++"
    if args.basis_set == ['pvdz']:
        filename_options['basis_set'] = "_pVDZ"  
    if args.basis_set == ['pvtz']:
        filename_options['basis_set'] = "_pVTZ"  
    if args.basis_set == ['pvqz']:
        filename_options['basis_set'] = "_pVQZ" 
    
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
        
    if args.freq == True:
        filename_options['frequency'] = "_freq"
    elif args.hindered_rotor == True:
        filename_options['frequency'] = "_hr"    
    elif args.freq == False and args.hindered_rotor == False:
        filename_options['frequency'] = ""
    return filename_options

def get_dot_com(args, options, filename_options):
    #JUST reactant or product
    if args.spc == False and args.qst3 == None:
        syntax="%mem={memory}GB\n%NProcShared={processors}\n%chk=min.chk\n#n opt=Z-Matrix NoSymm {density_functional}/{basis_set} {empirical_dispersion} {field_strength} {solvent} {frequency}\n\n {filename}\n\n{coordinates}\n{hindered_rotor}".format_map(options)
        completed_filename=("cis-stilbene_oxide{molecule_type}{density_functional}{basis_set}{solvent}{pathway}{enantiomer}{reactant_type}{field_strength}{frequency}.com").format_map(filename_options)        
    elif args.spc == True:
        if args.freq == False and args.hindered_rotor == False:
            syntax="%mem={memory}GB\n%NProcShared={processors}\n%chk=min.chk\n#n NoSymm {density_functional}/{basis_set} {empirical_dispersion} {field_strength} {solvent} Polar\n\n {filename}\n\n{coordinates}\n".format_map(options)
            completed_filename=("cis-stilbene_oxide{molecule_type}{density_functional}{basis_set}{solvent}{pathway}{enantiomer}{reactant_type}{field_strength}_spc_polar.com").format_map(filename_options)
        elif args.freq == True or args.hindered_rotor == True:
            syntax="%mem={memory}GB\n%NProcShared={processors}\n%chk=min.chk\n#n NoSymm {density_functional}/{basis_set} {empirical_dispersion} {field_strength} {solvent} {frequency}\n\n {filename}\n\n{coordinates}\n{hindered_rotor}".format_map(options)
            completed_filename=("cis-stilbene_oxide{molecule_type}{density_functional}{basis_set}{solvent}{pathway}{enantiomer}{reactant_type}{field_strength}_spc_{frequency}.com").format_map(filename_options)
    elif args.qst3 != None:
        syntax="%mem={memory}GB\n%NProcShared={processors}\n%chk=min.chk\n#n opt=(Z-Matrix,QST3,calcfc) NoSymm {density_functional}/{basis_set} {empirical_dispersion} {field_strength} {solvent} {frequency}\n\nStarting Material\n\n{starting_material}\n{product_type} Product\n\n{product}\nSaddle Point Guess\n\n{saddle_point}\n{hindered_rotor}".format_map(options)
        completed_filename=("cis-stilbene_oxide{molecule_type}{density_functional}{basis_set}{solvent}{pathway}{enantiomer}{reactant_type}{field_strength}{frequency}.com").format_map(filename_options)
    print("Prepared .com file name:", completed_filename)     
    with open(completed_filename, 'w') as f:
        f.write(syntax)

def aligned_to_com(aligned_filename, args):
    options = dict()
    filename_options = dict()
    get_density_functional(args, options)
    get_field_strength(args, options)
    get_solvent(args, options)
    get_basis_set(args, options)
    get_frequency(args, options)
    get_memory_and_processors(args,options)
    
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

def determine_filetype(args):
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

def main():
    args=get_arguments()
    determine_filetype(args)
 
    
if __name__ == "__main__":
    main()
