#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 14:36:02 2024

@author: samantha.bealldenn
"""

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
        
def main():
    args=get_arguments()
    #notation for getting one specific argument is args.argument/option
    options = dict()
    filename_options = dict()
    get_density_functional(args, options)
    get_field_strength(args, options)
    get_solvent(args, options)
    
    if args.qst3 == None:
        filename=args.filename[0]
        options['filename']=filename
        structure=get_coordinates(options, filename)
        options['coordinates'] = structure
        get_dotcom_filename(args, options, filename_options)
        get_dot_com(args, options, filename_options)
    
    elif args.qst3 != None:
        saddle_point=args.filename[0]
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

if __name__ == "__main__":
    main()
