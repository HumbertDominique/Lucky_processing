# Lecture et traitement d'image en lucky imaging.
# 
#
# 2023.06.26    -   Dominique Humbert   -   Version initiale

from astropy.io import fits
from astropy import stats
import matplotlib.pyplot as plt
import functions as fnc
import numpy as np
from scipy.interpolate import RegularGridInterpolator

from PIL import Image

import os
import time
import sys
from customClasses import eErrors


global eError
eError = eErrors.E_all_fine

k = len(sys.argv)
match k:
    case 1:
        print("here")
        a=1
        eError = eErrors.E_arg_error
    case 2:
        a = int(sys.argv[1])
        plotType = None
        eError = eErrors.E_all_fine
    case 3:
        a = int(sys.argv[1])
        plotType = np.array(sys.argv[2])
        eError = eErrors.E_all_fine
    case _:
        eError = eErrors.E_arg_error
        a=11


def main(a, eError, plotType):
    if (eError == eErrors.E_all_fine):
        st = time.time()
        # ------------------------------Variables------------------------------

        folder = 'Attempt2'
        filename = '/Attempt2_'
        path = folder+filename
        [k, j] = fits.getdata(path+'001.fits', ext=0).shape

        # -----------------------------count files-----------------------------
        # folder path
        dir_path = folder
        n = 0

        # Iterate directory
        for path in os.listdir(dir_path):
            # check if current path is a file
            if os.path.isfile(os.path.join(dir_path, path)):
                n += 1

        path = folder+filename


        data_raw = None
        data_roi = None
        data_sigma = None
        data_closed = None

        # -------------------------------program--------------------------------
        while eError==eErrors.E_all_fine:
            match a:
                case 0: # only read and save raw data
                    print('Read data')
                    data_raw = fnc.readData(n,path)
                    eError=eErrors.E_end_programm
                case 1:
                    print('Read data')
                    data_raw = fnc.readData(n,path)
                    np.save('data_raw',data_raw, allow_pickle=True, fix_imports=True)
                case 2:
                    print('ROI')
                    fnc.ROI(n,data_raw)
                case 3:
                    print("Dead pixels")
                    # https://www.astropy.org/ccd-reduction-and-photometry-guide/v/dev/notebooks/08-01-Identifying-hot-pixels.html
                    # for now, brute force
                    if data_roi is None:
                        print('loading ROI data')
                        data_roi = np.load('data_roi.npy')
                    #[k, j] = data_roi[:,:,1].shape
                    #data_sigma = np.zeros((k,j,n))
                case 4:
                    print('Sigma filter')
                    data_sigma = fnc.sigma(n,data_roi)

                case 5:
                    print("Speckles")
                    if data_sigma is None:
                        print('loading sigma data')
                        data_sigma = np.load('data_sigma.npy')
                    speckles = np.zeros((2,n))
                    for i in range(n):
                        speckles[:,i] = np.unravel_index(data_sigma[:,:,i].argmax(), data_sigma[:,:,i].shape)
                    print(speckles[0,0],speckles[1,0])
                    eError=eErrors.E_end_programm

                case 6:
                    print("mask")
                    if data_sigma is None:
                        print('loading sigma data')
                        data_sigma = np.load('data_sigma.npy')
                    maskR = 30  # [px]
                    mask = np.zeros((2*maskR,2*maskR))

                    k = data_sigma.shape[1]
                    x = np.linspace(0,k-1,num=k)
                    y = np.linspace(0,k-1,num=k)
                    interpolation = RegularGridInterpolator((x,y),data_sigma[:,:,0])

                    print(interpolation((x,y)).shape)
                    plt.figure()
                    plt.imshow()
                    plt.show()

                    eError=eErrors.E_end_programm
                case 7:
                    if data_sigma is None:
                        print('loading sigma data')
                        data_sigma = np.load('data_sigma.npy')
                    mask_radius = 100  # Rayon du cercle à exclure
                    approximated_image = fnc.polynomial_mask(data_sigma[:,:,1], mask_radius, 10,'nearest')
                    plt.figure()
                    plt.imshow(approximated_image, cmap='gray')
                    plt.show()
                case 10:
                    print('display')
                    fnc.display_images(plotType)
                    eError=eErrors.E_end_programm
                case _:
                    eError=eErrors.E_end_programm


            a += 1
            printTemp = 'eError:{}'.format(eError)
            #print(printTemp)



        et = time.time()
        # get the execution time
        elapsed_time = et - st
        print('Execution time:', elapsed_time, 'seconds')
    fnc.eError_handling(eError)


#-----------------KEEP HERE--------------------
main(a, eError, plotType)