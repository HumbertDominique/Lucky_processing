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




import sys
from customClasses import eErrors
import functions as fnc


#----------------------Global variables----------
global eError; eError = eErrors.E_all_fine

#------------------------------------------------

k = len(sys.argv)
match k:
    case 1:
        a = 1
        b = None
        r = None
        synt_BCG = False
        eError = eErrors.E_all_fine
    case 2:
        try: a = int(sys.argv[1])
        except: a = sys.argv[1]
        b = None
        r = None
        synt_BCG = False
        eError = eErrors.E_all_fine
        if isinstance(a,str):
            eError = eErrors.E_arg_error
    case 3:
        a = int(sys.argv[1])
        b = int(sys.argv[2])
        r = None
        synt_BCG = False
        eError = eErrors.E_all_fine
    case 4:
        a = int(sys.argv[1])
        b = int(sys.argv[2])
        r = float(sys.argv[3])
        synt_BCG = False
        eError = eErrors.E_all_fine
    case 5:
        a = int(sys.argv[1])
        b = int(sys.argv[2])
        r = float(sys.argv[3])
        synt_BCG = (sys.argv[4]=='True')
        eError = eErrors.E_all_fine
    case _:

        eError = eErrors.E_arg_error
        a=0

#-----------------KEEP HERE--------------------

fnc.lucky_process_focus(a,b,r,synt_BCG,eError)
