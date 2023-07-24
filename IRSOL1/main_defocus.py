# Lecture et traitement d'image en lucky imaging.
# 
# On first run, do:
#   - "python -m pip install numpy"
#   - "python -m pip install matplotlib"
#   - "python -m pip install from astropy"
#   - "python -m pip install from scipy"
# 
#
# Run from command line "python main.py arg1".
# Put no argument for help
#
# 2023.06.26    -   Dominique Humbert   -   Version initiale


#------------------LIBRARIES
from astropy.io import fits
from astropy import stats
import matplotlib.pyplot as plt
from matplotlib import cm
import functions as fnc
import numpy as np
from matplotlib.ticker import LinearLocator
from scipy import ndimage



import os
import time
import sys
from customClasses import eErrors


#----------------------Global variables----------
global eError; eError = eErrors.E_all_fine

#------------------------------------------------

k = len(sys.argv)
match k:
    case 1:
        a = 0
        synt_BCG = True
        eError = eErrors.E_arg_error
    case 2:
        a = int(sys.argv[1])
        synt_BCG = False
        eError = eErrors.E_all_fine
    case 3:
        a = int(sys.argv[1])
        synt_BCG = (sys.argv[2]=='True')
        eError = eErrors.E_all_fine
    case _:
        eError = eErrors.E_arg_error
        a=0



#-----------------KEEP HERE--------------------

fnc.lucky_process_defocus(a,synt_BCG,eError)
