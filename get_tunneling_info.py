#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 14:15:48 2024

@author: samantha.bealldenn
"""

import numpy as np

def get_tunneling_coefficients(frequency):
    T = 298
    k_b = 1.381*10**(-23)
    h = 6.626*10**(-34)
    c = 2.998*10**(10)
    beta = 1/(T*k_b)
    h_bar_omega = h * c * frequency
    
    standard_tunneling_coefficient = (beta*h_bar_omega*0.5)/np.sin(beta*h_bar_omega*0.5)
    wigner_tunneling_coeff = 1+(1/24)*(beta*h_bar_omega)**2
    
    
    print("Standard Wigner Tunneling Coefficient:", standard_tunneling_coefficient)
    print("Truncated Wigner Tunneling Coefficient:", wigner_tunneling_coeff)
    tunneling_info = [h_bar_omega, standard_tunneling_coefficient, wigner_tunneling_coeff]
    return(tunneling_info)
    

def get_crossover_temp(h_bar_omega):
    k_b = 1.381*10**(-23)
    T = 298
    crossover_temp = h_bar_omega/(2*np.pi*k_b)
    if crossover_temp > (2/3)*T:
        print("Warning! Crossover temperature is too close to room temperature!")

    
def get_tunneling_information(frequency):
    tunneling_info=get_tunneling_coefficients(frequency)
    get_crossover_temp(tunneling_info[0])
    
    
