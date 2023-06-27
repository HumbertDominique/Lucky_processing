# Lecture d'image en lucky imaging.
# 
#
# 2023.06.26    -   Dominique Humbert   -   Version initiale
from astropy.io import fits
import numpy as np
from customClasses import eErrors



def readData(n, fname):

    #data = fits.getdata(fname+'001.fits', ext=0)



    for i in range(n):
        filename = fname+'{:03d}'.format(i+1)+'.fits'
        if i==0:                    # Gets the pic size on firts pass
            [k, j] = fits.getdata(fname+'001.fits', ext=0).shape
            data = np.zeros((k,j,n))

        data[:,:,i] = fits.getdata(filename, ext=0)
    return data

def eError_handling(code):
    match code:
        case eErrors.E_all_fine:
            print("Programm terminated with success")
        case eErrors.E_arg_error:
            print("--------------------------")
            print("Too many arguments")
            print("arg1 = process\n1: Whole process\n2: ROI\n3: Sigma filter\n10: plot")
            print("arg2 = plot type")
            print("--------------------------")


