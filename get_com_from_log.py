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
    # field strength should be nofield, nf, negx, nx, px, or posx, with x being an integer
    parser.add_argument("field_strength", nargs=1, action="store")
    #gets solvent from command line
    parser.add_argument("solvent", choices=["acn", "dcm", "benz", "tcm", "nosolv"], nargs=1, action='store')
    parser.add_argument("basis_set", choices=["+", "++", "pvdz", "pvtz", "pvqz"], nargs=1, action="store")
    ## gives chance to change the memory and number of processors
    parser.add_argument("-m", "--memory", nargs=2, action='store')
    ## gives chance to add a comment to the end of the filename from the command line
    parser.add_argument("-c", "--comment", nargs=1, action='store')
    
    group = parser.add_mutually_exclusive_group(required=False)
    #changes the regular minimization to a single point calculation
    group.add_argument("-s", "--spc", action="store_true")
    #adds the option to make the output file QST3, in which case we want to get the names of the other two files being used
    group.add_argument("-q", "--qst3", nargs=2, action='store')
    # adds the option to do a transition state calculation, which does not require a reactant and product
    group.add_argument("-ts", "--transition_state", action="store_true")
    
    #adds the option to include frequency calculations, which factor entropy into the free energy calculations.
    # can't have just frequency and hindered rotor
    group2 = parser.add_mutually_exclusive_group(required=False)
    
    group2.add_argument("-nf", "--nofreq", action="store_false")
    # make frequency the default, so you can do -nf to take away frequency
    group2.add_argument("-hr", "--hindered_rotor", action="store_true")

    args=parser.parse_args()
    return args
    #args is the entire namespace

#density functional options
def get_density_functional(args, options, filename_options):
    if args.density_functional.__contains__('b3lyp'):
        options['density_functional'] = "B3LYP"
        options['empirical_dispersion'] = "EmpiricalDispersion=GD3"
        filename_options['density_functional'] = '_b3lyp'
    elif args.density_functional.__contains__('mn15'):
        options['density_functional'] = "MN15"
        options['empirical_dispersion'] = ""
        filename_options['density_functional'] = '_mn15'
    elif args.density_functional.__contains__('b2plypd3'):
        options['density_functional'] = "B2PLYPD3"
        options['empirical_dispersion'] = ""
        filename_options['density_functional'] = '_b2plypd3'
    elif args.density_functional.__contains__('b2plyp'):
        options['density_functional'] = "B2PLYP"
        options['empirical_dispersion'] = ""
        filename_options['density_functional'] = '_b2plyp'         
    return(options, filename_options)

    
#field strength options    
def get_field_strength(args, options, filename_options):
    electric_field = args.field_strength[0]
    ## EVALUATION OF INPUT MUST BE IN THIS ORDER
    ## "n" will be contained if the input is nf/nofield/neg/n, so we check for n after we check for nofield/nf
    ## instead of matching our input to a list of options, we will take the field and split it into the direction and strength
    ## then recombine it with the proper syntax to be in the file name or in the .com file.
    if electric_field.__contains__("nofield") or electric_field.__contains__("nf"):
        filename_options['field_strength'] = "_nofield"
        options['field_strength'] = ""
    elif electric_field.__contains__("neg"):
        field_direction = "Field=Z-"
        field_direction_for_filename = "_neg"
        field_strength = electric_field.split("g")[1]
        filename_options['field_strength'] = field_direction_for_filename + field_strength
        options['field_strength'] = field_direction + field_strength
    elif electric_field.__contains__("n"):
        field_direction = "Field=Z-"
        field_direction_for_filename = "_neg"
        field_strength = electric_field.split("n")[1]
        filename_options['field_strength'] = field_direction_for_filename + field_strength
        options['field_strength'] = field_direction + field_strength
    elif electric_field.__contains__("pos"):
        field_direction = "Field=Z+"
        field_direction_for_filename = "_pos"
        field_strength = electric_field.split("s")[1]
        filename_options['field_strength'] = field_direction_for_filename + field_strength
        options['field_strength'] = field_direction + field_strength
    elif electric_field.__contains__("p"):
        field_direction = "Field=Z+"
        field_direction_for_filename = "_pos"
        field_strength = electric_field.split("p")[1]
        filename_options['field_strength'] = field_direction_for_filename + field_strength
        options['field_strength'] = field_direction + field_strength
    return(options, filename_options)
    
#solvent options
def get_solvent(args, options, filename_options):
    if args.solvent == ['acn']:
        options['solvent'] = "SCRF=(Solvent=Acetonitrile)"
        filename_options['solvent'] = '_acn'
    elif args.solvent == ['dcm']:
        options['solvent'] = "SCRF=(Solvent=Dichloromethane)"
        filename_options['solvent'] = '_dcm'
    elif args.solvent == ['nosolv']:
        options['solvent'] = ""
        filename_options['solvent'] = "_nosolv"
    elif args.solvent == ['benz']:
        options['solvent'] = "SCRF=(Solvent=Benzene)"
        filename_options['solvent'] = "_benzene"
    elif args.solvent == ['tcm']:
        options['solvent'] = "SCRF=(Solvent=Chloroform)"
        filename_options['solvent'] = "_tcm"
    else:
        print("Unrecognized solvent")        
    return(options, filename_options)

def get_basis_set(args, options, filename_options):
    if args.basis_set == ['+']:
        options['basis_set'] = "6-311+g(d,p)"
        filename_options['basis_set'] = ""
    elif args.basis_set == ['++']:
        options['basis_set'] = "6-311++g(d,p)"
        filename_options['basis_set'] = "_++"
    elif args.basis_set == ['pvdz']:
        options['basis_set'] = "AUG-cc-pVDZ"  
        filename_options['basis_set'] = "_pVDZ" 
    elif args.basis_set == ['pvtz']:
        options['basis_set'] = "AUG-cc-pVTZ"
        filename_options['basis_set'] = "_pVTZ"
    elif args.basis_set == ['pvqz']:
        options['basis_set'] = "AUG-cc-pVQZ"
        filename_options['basis_set'] = "_pVQZ"
    return(options, filename_options)

def get_frequency(args, options, filename_options):
    if args.nofreq == True:
        options['frequency'] = "Freq"
        options['hindered_rotor'] = ""
        filename_options['frequency'] = "_freq"
    elif args.nofreq == False and args.hindered_rotor == False:
        options['frequency'] = ""
        options['hindered_rotor'] = ""
        filename_options['frequency'] = ""
    elif args.hindered_rotor == True:
        options['frequency'] = "Freq=(HinderedRotor,ReadHinderedRotor)"
        options['hindered_rotor'] = "1.0\n1 2 3 1 1\n\n"
        filename_options['frequency'] = "_hr"
    return(options, filename_options)

def get_memory_and_processors(args, options):
    if args.memory == None:
        options['memory'] = "24"
        options['processors'] = "32"
    if args.memory != None:
        options['memory'] = args.memory[0]
        options['processors'] = args.memory[1]
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
    if args.qst3 == None and args.transition_state == False:
        if filename.__contains__("product"):
            filename_options['molecule_type'] = "_product"
        elif filename.__contains__("reactant"):
            filename_options['molecule_type'] = "_reactant"
        elif filename.__contains__("epoxide"):
            filename_options['molecule_type'] = "_epoxide"
        else:
            filename_options['molecule_type'] = "_other"
    elif args.qst3 == None and args.transition_state == True:
        filename_options['molecule_type'] = "_ts"
    elif args.qst3 != None:   
        filename_options['molecule_type'] = "_qst3"
   
    ## adds aldehyde and ketone to the name of the file
    ## in the case of QST3 files, adds aldehyde and ketone to the file    
    if filename.__contains__("ketone"):
        filename_options['pathway'] = "_ketone"
        options['product_type'] = "Ketone"
    elif filename.__contains__("aldehyde"):
        filename_options['pathway'] = "_aldehyde"
        options['product_type'] = "Aldehyde"
    elif filename.__contains__("epoxide"):
        filename_options['pathway'] = "_epoxide"
        options['product_type'] = "Epoxide"
    else:
        filename_options['pathway'] = ""
        options['product_type'] = "Other"
        
    if filename.__contains__("_s_"):
        filename_options['enantiomer'] = "_s"
    elif filename.__contains__("_r_"):
        filename_options['enantiomer'] = "_r"
    else:
        filename_options['enantiomer'] = ""
        
    if filename.__contains__("_longarm1_"):
        filename_options['reactant_type'] = "_longarm1"
    elif filename.__contains__("_longarm2_"):
        filename_options['reactant_type'] = "_longarm2"
    elif filename.__contains__("_longarm_"):
        filename_options['reactant_type'] = "_longarm"    
    elif filename.__contains__("_shortarm_"):
        filename_options['reactant_type'] = "_shortarm"
    elif filename.__contains__("_shortarm1_"):
        filename_options['reactant_type'] = "_shortarm1"
    elif filename.__contains__("_shortarm2_"):
        filename_options['reactant_type'] = "_shortarm2"
    elif filename.__contains__("_shortarm_"):
        filename_options['reactant_type'] = "_shortarm"
    else:
        filename_options['reactant_type'] = ""
        
    if args.comment != None:
        comment = args.comment[0]
        filename_options['comment'] = "_" + comment
    elif args.comment == None:
        filename_options['comment'] = ""
    return(filename_options, options)

def get_dot_com(args, options, filename_options):
    #JUST reactant or product
    if args.spc == False and args.qst3 == None and args.transition_state == False:
        syntax="%mem={memory}GB\n%NProcShared={processors}\n%chk=min.chk\n#n opt=Z-Matrix NoSymm {density_functional}/{basis_set} {empirical_dispersion} {field_strength} {solvent} {frequency}\n\n {filename}\n\n{coordinates}\n{hindered_rotor}".format_map(options)
        completed_filename=("cis-stilbene_oxide{molecule_type}{density_functional}{basis_set}{solvent}{pathway}{enantiomer}{reactant_type}{field_strength}{frequency}{comment}.com").format_map(filename_options)        
    elif args.spc == True:
        if args.freq == False and args.hindered_rotor == False:
            syntax="%mem={memory}GB\n%NProcShared={processors}\n%chk=min.chk\n#n NoSymm {density_functional}/{basis_set} {empirical_dispersion} {field_strength} {solvent} Polar\n\n {filename}\n\n{coordinates}\n".format_map(options)
            completed_filename=("cis-stilbene_oxide{molecule_type}{density_functional}{basis_set}{solvent}{pathway}{enantiomer}{reactant_type}{field_strength}_spc_polar{comment}.com").format_map(filename_options)
        elif args.freq == True or args.hindered_rotor == True:
            syntax="%mem={memory}GB\n%NProcShared={processors}\n%chk=min.chk\n#n NoSymm {density_functional}/{basis_set} {empirical_dispersion} {field_strength} {solvent} {frequency}\n\n {filename}\n\n{coordinates}\n{hindered_rotor}".format_map(options)
            completed_filename=("cis-stilbene_oxide{molecule_type}{density_functional}{basis_set}{solvent}{pathway}{enantiomer}{reactant_type}{field_strength}_spc_{frequency}{comment}.com").format_map(filename_options)
    elif args.qst3 != None:
        syntax="%mem={memory}GB\n%NProcShared={processors}\n%chk=min.chk\n#n opt=(Z-Matrix,QST3,calcfc) NoSymm {density_functional}/{basis_set} {empirical_dispersion} {field_strength} {solvent} {frequency}\n\nStarting Material\n\n{starting_material}\n{product_type} Product\n\n{product}\nSaddle Point Guess\n\n{saddle_point}\n{hindered_rotor}".format_map(options)
        completed_filename=("cis-stilbene_oxide{molecule_type}{density_functional}{basis_set}{solvent}{pathway}{enantiomer}{reactant_type}{field_strength}{frequency}{comment}.com").format_map(filename_options)
    elif args.transition_state == True:
        syntax="%mem={memory}GB\n%NProcShared={processors}\n%chk=min.chk\n#n opt=(Z-Matrix,TS,calcfc) NoSymm {density_functional}/{basis_set} {empirical_dispersion} {field_strength} {solvent} {frequency}\n\n {filename}\n\n{coordinates}\n{hindered_rotor}".format_map(options)
        completed_filename=("cis-stilbene_oxide{molecule_type}{density_functional}{basis_set}{solvent}{pathway}{enantiomer}{reactant_type}{field_strength}{frequency}{comment}.com").format_map(filename_options)        
    print("Prepared .com file name:", completed_filename)     
    with open(completed_filename, 'w') as f:
        f.write(syntax)

def aligned_to_com(aligned_filename, args):
    options = dict()
    filename_options = dict()
    get_density_functional(args, options, filename_options)
    get_field_strength(args, options, filename_options)
    get_solvent(args, options, filename_options)
    get_basis_set(args, options, filename_options)
    get_frequency(args, options, filename_options)
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
