# Lecture et traitement d'image en lucky imaging.
# 
# On first run, do:
#   - "python -m pip install numpy"
#   - "python -m pip install matplotlib"
#   - "python -m pip install from astropy"
#   - "python -m pip install from scipy"
#
# Run from command line "python main.py arg1 arg2".
# Put no argument for help
#
# 2023.06.26    -   Dominique Humbert   -   Version initiale


#------------------LIBRARIES (temporary)
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
        a=7
        plot_type = None
        eError = eErrors.E_arg_error
    case 2:
        a = int(sys.argv[1])
        plot_type = None
        eError = eErrors.E_all_fine
    case 3:
        a = int(sys.argv[1])
        plot_type = np.array(sys.argv[2])
        eError = eErrors.E_all_fine
    case _:
        eError = eErrors.E_arg_error
        a=11


def main(a, eError=eErrors.E_all_fine):
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
        # -----------------------------Debug-----------------------------
        # -----------------------------Debug-----------------------------
        # -----------------------------Debug-----------------------------
        # -----------------------------Debug-----------------------------
        n = 5
        # -----------------------------Debug-----------------------------
        # -----------------------------Debug-----------------------------
        # -----------------------------Debug-----------------------------
        # -----------------------------Debug-----------------------------

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
                case 2:
                    print('ROI')
                    fnc.ROI(n,data_raw)
                    del data_raw
                case 3:
                    print("Dead pixels")
                    # https://www.astropy.org/ccd-reduction-and-photometry-guide/v/dev/notebooks/08-01-Identifying-hot-pixels.html
                    # for now, brute force
                    if data_roi is None:
                        print('loading ROI data')
                        data_roi = np.load('temp/data_roi.npy')
                    #[k, j] = data_roi[:,:,1].shape
                    #data_sigma = np.zeros((k,j,n))
                case 4:
                    print('Sigma filter')
                    data_sigma = fnc.sigma(n,data_roi)
                    del data_roi
                case 5:
                    print("Speckles")
                    if data_sigma is None:
                        print('loading sigma data')
                        data_sigma = np.load('temp/data_sigma.npy')
                    speckles = np.zeros((2,n))
                    for i in range(n):
                        speckles[:,i] = np.unravel_index(data_sigma[:,:,i].argmax(), data_sigma[:,:,i].shape)

                case 6:
                    print("mask")
                    if data_sigma is None:
                        print('loading sigma data')
                        data_sigma = np.load('temp/data_sigma.npy')
                    mask_radius = 150  # Rayon du cercle à exclure


        
                    approx_data = fnc.polynomial_mask_py(data_sigma[:,:,0], mask_radius,'nearest')
                    data_no_bcg = data_sigma[:,:,0] - approx_data
                    np.save('temp/data_no_bcg',data_no_bcg, allow_pickle=True, fix_imports=True)

                    plt.figure()
                    plt.pcolor(data_sigma[:,:,0],vmin=0, vmax=4095)
                    plt.title('sigma filter')
                    plt.colorbar()

                    plt.figure()
                    plt.pcolor(approx_data,vmin=0, vmax=4095//2)
                    plt.title('approx bcg')
                    plt.colorbar()

                    plt.figure()
                    plt.pcolor(data_no_bcg,vmin=0, vmax=4095)
                    plt.title('no bcg')
                    plt.colorbar()

                case 7:
                    if data_sigma is None:
                        print('loading sigma data')
                        data_sigma = np.load('temp/data_sigma.npy')
                    mask_radius = 150  # Rayon du cercle à exclure
                    
                    mask = fnc.buid_mask(data_sigma[:,:,0],mask_radius)
                    # plt.figure()
                    # plt.pcolor(data_sigma[:,:,4], cmap='gray',vmin=0, vmax=4095)
                    # plt.title('before')
                    # plt.colorbar()

                    plt.figure()
                    plt.imshow(data_sigma[:,:,4])
                    plt.title('before')
                    plt.colorbar()

                    plt.figure()
                    plt.imshow(mask)
                    plt.title('mask')
                    plt.colorbar()

                    model = fnc.polynomial_mask(data_sigma[:,:,4],mask,4)
                    data_sigma[:,:,4] = data_sigma[:,:,4] - model
                    
                    plt.figure()
                    plt.imshow(model)
                    plt.title('model')
                    plt.colorbar()

                    plt.figure()
                    plt.imshow(data_sigma[:,:,4])
                    plt.title('after')
                    plt.colorbar()
                    plt.show()

                case 8:
                    a
                case 9:
                    a

                case 10:
                    et = time.time()
                    # get the execution time
                    elapsed_time = et - st
                    print('Execution time:', elapsed_time, 'seconds')
                    
                    print('display')
                    
                    #fnc.display_images(data_sigma[:,:,0:2])
                    eError=eErrors.E_end_programm
                case _:
                    et = time.time()
                    # get the execution time
                    elapsed_time = et - st
                    print('Execution time:', elapsed_time, 'seconds')

                    eError=eErrors.E_end_programm
                

            a += 1
            printTemp = 'eError:{}'.format(eError)
            #print(printTemp)

    fnc.eError_handling(eError)


#-----------------KEEP HERE--------------------

main(a)