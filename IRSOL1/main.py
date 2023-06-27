# Lecture et traitement d'image en lucky imaging.
# 
#
# 2023.06.26    -   Dominique Humbert   -   Version initiale

from astropy.io import fits
from astropy import stats
import matplotlib.pyplot as plt
import functions as fnc
import numpy as np
import os
import time
import sys
from customClasses import eErrors


eError = eErrors.E_all_fine

k = len(sys.argv)
match k:
    case 1:
        print("here")
        a=1
        eError = eError.E_arg_error
    case 2:
        a = int(sys.argv[1])
    case 3:
        a = int(sys.argv[1])
        plotType = np.array(sys.argv[2])
if k>=4:
    eError = eError.E_arg_error
    a=11
def main(a):
    if eError == eError.E_all_fine:
        print("here")
        st = time.time()
        # ------------------------------Variables------------------------------
        state = 0

        folder = 'Attempt1'
        filename = '/Attempt1_'
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

    # -----------------------------Debug-----------------------------
    # -----------------------------Debug-----------------------------
    # -----------------------------Debug-----------------------------

        n = 10                          # limite du fichiers pour debug

    # -----------------------------Debug-----------------------------
    # -----------------------------Debug-----------------------------
    # -----------------------------Debug-----------------------------

        data_raw = None
        data_roi = None
        data_sigma = None

        # -------------------------------program--------------------------------
        while state==0:
            match a:
                case 0: # only read and save raw data
                    print('Read data')
                    data_raw = fnc.readData(n,path)
                    np.save('data_raw',data_raw, allow_pickle=True, fix_imports=True)
                    state=1
                case 1:
                    print('Read data')
                    data_raw = fnc.readData(n,path)
                    np.save('data_raw',data_raw, allow_pickle=True, fix_imports=True)
                case 2:
                    print('ROI')
                    if data_raw is None:
                        print('loading raw data')
                        data_raw = np.load('data_raw.npy')
                    data_roi = data_raw[424-295:424+294,600-295:600+294,:]     # Manually done
                    np.save('data_roi',data_roi, allow_pickle=True, fix_imports=True)
                case 3:
                    print('Sigma filter')
                    if data_roi is None:
                        print('loading ROI data')
                        data_roi = np.load('data_roi.npy')
                    [k, j] = data_roi[:,:,1].shape
                    data_sigma = np.zeros((k,j,n))
                    for i in range(n):
                        data_sigma[:,:,i] = stats.sigma_clip(data_roi[:,:,i])
                    np.save('data_sigma',data_sigma, allow_pickle=True, fix_imports=True)

                    state=1

                case 10:
                    print('display')
                    if plotType == "raw":
                        if data_raw is None:
                            print('loading raw data')
                            data_raw = np.load('data_raw.npy')
                        plt.figure()
                        plt.imshow(data_raw[:,:,0], cmap='gray')
                        plt.show()
                    state=1
                    if plotType == "roi":
                        if data_roi is None:
                            print('loading ROI data')
                            data_roi = np.load('data_roi.npy')
                        plt.figure()
                        plt.imshow(data_roi[:,:,0], cmap='gray')
                        plt.show()
                    state=1


            a += 1
            printTemp = 'state:{:d}'.format(state)
            print(printTemp)

            #data = fnc.readData(2,filename)




        et = time.time()
        # get the execution time
        elapsed_time = et - st
        print('Execution time:', elapsed_time, 'seconds')
    fnc.eError_handling(eError)


#-----------------KEEP HERE--------------------
main(a)