#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 10:39:20 2024

@author: samantha.bealldenn
"""

from AaronTools.comp_output import CompOutput
from AaronTools.fileIO import FileReader

def get_entropy_corrected_G(filename, temperature=None, w0=100.):
    infile = FileReader(filename, just_geom=False)
    co = CompOutput(infile)

    nrg = co.energy

    qrrho_dG = co.calc_G_corr(v0=w0, temperature=temperature, method="QRRHO")
    qharm_dG = co.calc_G_corr(v0=w0, temperature=temperature, method="QHARM")

    quasi_rrho_G = (nrg + qrrho_dG) * 2625.5
    quasi_harmonic_G = (nrg + qharm_dG) * 2625.5
    
    print("Quasi-Rho G:", quasi_rrho_G)
    print("Quasi-Harmonic G", quasi_harmonic_G)
    
    entropy_corrected_G = [quasi_rrho_G, quasi_harmonic_G]
    
    return entropy_corrected_G