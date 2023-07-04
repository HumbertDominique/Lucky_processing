# Lecture et traitement d'image en lucky imaging.
# 
# On first run, do:
#   - "python -m pip install numpy"
#   - "python -m pip install matplotlib"
#   - "python -m pip install from astropy"
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
#global r; r = 0.2     # [%]
#------------------------------------------------

k = len(sys.argv)
match k:
    case 1:
        a = 0
        b = 0
        eError = eErrors.E_arg_error
    case 2:
        a = int(sys.argv[1])
        b = None
        eError = eErrors.E_all_fine
    case 3:
        a = int(sys.argv[1])
        b = int(sys.argv[2])
        eError = eErrors.E_all_fine
    case 4:
        a = int(sys.argv[1])
        b = int(sys.argv[2])
        r = float(sys.argv[3])
        eError = eErrors.E_all_fine
    case _:
        eError = eErrors.E_arg_error
        a=0


def lucky_process(a, b, r=0.2, eError=eErrors.E_arg_error):
    '''------------------------------------------------------------------------------------
    # main(a, b, eError=eErrors.E_arg_error)
    Processes lucky imaging frames
    ## in:
      a: number of images
      b: 
      c:
    ## out:
      data: 3D array with pixel data
    ------------------------------------------------------------------------------------'''

    if (eError == eErrors.E_all_fine):
        st = time.time()
        # ------------------------------Variables------------------------------
        data = None
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

        # ----------------------------------------------------------
        if b is None:       # qty of frames to use
            b = 20
        else:
            n = b
        # ----------------------------------------------------------

        path = folder+filename

        # -------------------------------program--------------------------------
        while eError==eErrors.E_all_fine:
            match a:
                case 0: # only read and save raw data
                    print('Read data')
                    data = fnc.readData(n,path)
                    mean_image_raw = np.sum(data,2)/n
                    print('Saving intermediate data')
                    np.save('temp/mean_image_raw',mean_image_raw, allow_pickle=True, fix_imports=True)
                    print('Done')  
                    eError=eErrors.E_end_programm
                case 1:
                    print('Read data')
                    data = fnc.readData(n,path)
                    mean_image_raw = np.sum(data,2)/n
                    print('Saving intermediate data')
                    np.save('temp/mean_image_raw',mean_image_raw, allow_pickle=True, fix_imports=True)
                    print('Done')
                    # plt.figure()
                    # plt.imshow(data[:,:,0],cmap='gray',vmin=0,vmax=4095)
                    # plt.title('Raw image')
                    # plt.figure()
                    
                    # plt.imshow(mean_image_raw,cmap='gray',vmin=0,vmax=4095)
                    # text_temp = 'Raw mean on {n:d} samples'
                    # plt.title(text_temp.format(n=n))
                    #plt.show()
                case 2:
                    print('ROI')
                    data = fnc.ROI(n,data)      # new variable because of region of interest
                    mean_image_roi = np.sum(data,2)/n
                    # plt.figure()
                    # plt.imshow(data[:,:,0],cmap='gray',vmin=0,vmax=4095)
                    # plt.title('ROI image')

                    # plt.figure()
                    # plt.imshow(mean_image_roi,cmap='gray')
                    # text_temp = 'ROI mean on {n:d} samples'
                    # plt.title(text_temp.format(n=n))
                    # plt.show()
                    
                case 3:
                    print('Sigma filter')
                    if data is None:
                        print('Loading ROI data')
                        data = np.load('temp/data_roi.npy')
                    
                    status = 0
                    for i in range(0,n):
                        data[0,0,i] = stats.sigma_clip(data[:,:],maxiters=None)
                        if 100*i//n == status*10:
                            print('Sigma filtering: ',status*10,'% done')
                            status += 1
                    np.save('temp/data_sigma',data, allow_pickle=True, fix_imports=True)
                    print('Done')
                    mean_image_sigma = np.sum(data,2)/n
                    # plt.figure()
                    # plt.imshow(data[:,:,0],cmap='gray',vmin=0,vmax=4095)
                    # plt.title('Sigma image')

                    # plt.figure()
                    # plt.imshow(mean_image_sigma,cmap='gray',vmin=0,vmax=4095)
                    # text_temp = 'Sigma mean on {n:d} samples'
                    # plt.title(text_temp.format(n=n))



                    
                case 4:
                    print('Normalisation')
                    if data is None:
                        print('Loading sigma data')
                        data = np.load('temp/data_sigma.npy')

                    status = 0
                    for i in range(0,n):
                        data[:,:,i] = fnc.normalise(data[:,:,i],12)
                        if 100*i//n == status*10:
                            print('Sigma filtering: ',status*10,'% done')
                            status += 1

                    print('Saving intermediate data')
                    np.save('temp/data_norm',data, allow_pickle=True, fix_imports=True)
                    print('Done')
                    mean_image_norm = np.sum(data,2)/n

                
                case 5:
                    print("Background")
                    if data is None:
                        print('loading sigma data')
                        data = np.load('temp/data_sigma.npy')
                        
                    mask_radius = 150  # Rayon du cercle à exclure
                    
                    print("Speckles location")
                    speckles = np.zeros((2,n))
                    status = 0
                    for i in range(n):
                        speckles[:,i] = np.unravel_index(data[:,:,i].argmax(), data[:,:,i].shape)
                        if 100*i//n == status*10:
                            print('Locating brightest speckles: ',status*10,'% done')
                            status += 1


                    mask = np.zeros((data.shape))  # Remove if using mask and model at the same time as data takes too much
                    model = np.zeros(data.shape)

                    status = 0
                    for i in range(0,n):
                        mask[:,:,i] = fnc.buid_mask(data[:,:,0],mask_radius,speckles[:,i])
                        model[:,:,i] = fnc.polynomial_mask(data[:,:,i],mask[:,:,i],4)
                        data[:,:,i] = data[:,:,i] - model[:,:,i]
                        #mask = fnc.buid_mask(data[:,:,0],mask_radius)
                        #model = fnc.polynomial_mask(data[:,:,i],mask,4)
                        #data[:,:,i] = data[:,:,i] - model
                        if 100*i//n == status*10:
                            print('background removal: ',status*10,'% done')
                            status += 1

                    print('Saving intermediate data')
                    np.save('temp/data_mask',mask, allow_pickle=True, fix_imports=True)
                    del mask
                    np.save('temp/data_model',model, allow_pickle=True, fix_imports=True)
                    del model
                    np.save('temp/data_bcg',data, allow_pickle=True, fix_imports=True)
                    print('Done')

                    mean_image_bcg = np.sum(data,2)/n 
                    np.save('temp/mean_image_bcg',mean_image_bcg, allow_pickle=True, fix_imports=True)
                    print('Done')
                    # plt.figure()
                    # plt.imshow(data[:,:,0])
                    # plt.title('Bcg removed image')
                    # plt.show()
                    # plt.figure()
                    # plt.imshow(mean_image_bcg,cmap='gray',vmin=0,vmax=4095)
                    # text_temp = 'Bcg removed mean on {n:d} samples'
                    # plt.title(text_temp.format(n=n))
                    #plt.show()
                case 6:
                    print('Noise cut')
                    # if data is None:
                    #     print('Loading bcg data')
                    #     data = np.load('temp/data_sigma.npy')
                    #     k, j, n = data.shape
                    #     data = np.zeros((k,j,n))
                    # data[:,:,0] = fnc.noise_cut(data[:,:,0])
                    # plt.figure()
                    # plt.imshow(data[:,:,0])
                    # text_temp = 'Noise cut'
                    # plt.title(text_temp)
                    #plt.show()
                case 7:
                    a    
                case 8:
                    a
                    # print('Noise filtering')
                    # if data is None:
                    #     print('Loading some data')
                    #     data = np.load('temp/data_bcg.npy')
                    #     k, j, n = data.shape

                    # unfiltered = data[:,:,0]
                    # plt.figure()
                    # plt.imshow(unfiltered,cmap='gray')
                    # text_temp = 'Unfiltered'
                    # plt.title(text_temp)

                    # for i in range(0,n):
                    #     data[:,:,0] = fnc.noise_filter(data[:,:,0],1)
                    # print('Saving intermediate data')
                    # np.save('temp/data_filtered',data, allow_pickle=True, fix_imports=True)
                    # print('Done')
                    # plt.figure()
                    # plt.imshow(data[:,:,0],cmap='gray')
                    # text_temp = 'Noise filtered'
                    # plt.title(text_temp)
                    # plt.show()
                case 9:
                    print('Dust spots')
                    NdustSpots = 4
                    spot_model_radius = 25  
                    status = 0;
                    if data is None:
                        print('Loading filtered data')
                        data = np.load('temp/data_bcg.npy')

                    try: mean_image_bcg
                    except NameError: mean_image_bcg = np.load('temp/mean_image_bcg.npy')
                    # plt.figure()
                    # plt.imshow(mean_image_bcg)
                    # text_temp = 'Shadows everywhere!!!'
                    # plt.title(text_temp)
                    # find a dust spot for the model
                    index_model = np.unravel_index(mean_image_bcg[:,:].argmin(), mean_image_bcg[:,:].shape)
                    dust_model = mean_image_bcg[index_model[0]-spot_model_radius:index_model[0]+spot_model_radius,index_model[1]-spot_model_radius:index_model[1]+spot_model_radius]
                    dust_model /= np.min(dust_model)
                    
                    for i in range(0,n):
                        data[:,:,i] = fnc.image_dusting(data[:,:,i],dust_model,4)
                        if 100*i//n == status*10:
                            print('Dust shadow removal: ',status*10,'% done')
                            status += 1
                    
                    print('Saving intermediate data')
                    np.save('temp/data_dustfree',data, allow_pickle=True, fix_imports=True)

                    mean_image_dustfree = np.sum(data,2)/n 
                    np.save('temp/mean_image_dustfree',mean_image_dustfree, allow_pickle=True, fix_imports=True)
                    print('Done')

                    # plt.figure()
                    # plt.imshow(mean_image_dustfree)
                    # text_temp = 'Dust shadow removed'
                    # plt.title(text_temp)
                    # plt.show()
                case 10:
                    print('tracking')  
                    if data is None:
                        print('Loading shadowless data')
                        data = np.load('temp/data_dustfree.npy') 
                    status = 0

                    for i in range(n):
                        speckles = np.unravel_index(data[:,:,i].argmax(), data[:,:,i].shape)
                        if i == 0:
                            ref = speckles
                        Dx = ref[0]-speckles[0] 
                        Dy = ref[1]-speckles[1] 
                        data[:,:,i] = np.roll(data[:,:,i],(Dy,Dx),axis=(1,0))
                        if 100*i//n == status*10:
                            print('Tracking: ',status*10,'% done')
                            status += 1

                    print('Saving intermediate data')
                    np.save('temp/data_Stracked',data, allow_pickle=True, fix_imports=True)

                    mean_image_Stracking = np.sum(data,2)/n
                    np.save('temp/mean_image_Stracking',mean_image_Stracking, allow_pickle=True, fix_imports=True)
                    print('Done')
                    # plt.figure()
                    # plt.imshow(mean_image_Stracking)
                    # text_temp = 'Speckle tracking'
                    # plt.title(text_temp)
                    # plt.show()

                case 11:
                    print('Selection')
                    if data is None:
                        print('Loading Stracked data')
                        data = np.load('temp/data_Stracked.npy') 

                    #r = 0.2 # [%] # now a global variable
                    text_temp = 'Selecting the {r:.2f}% best frames'
                    print(text_temp.format(r=100*r))

                    # Find the most luminous pixel on each frame
                    max_pixel_light = np.zeros((n))
                    for i in range(0,n):
                        max_pixel_light[i] = np.amax(data[:,:,i])
                    order = np.argsort(max_pixel_light)

                    # sort Luminosity  the info 
                    n_needed = r*n
                    n_needed = int(-(-n_needed//1))   # Ceil without math module. We need a finite number of frame
                    selected_indices = order[n-n_needed:n]

                    # Select only the best
                    data_best  = data[:,:,selected_indices]

                    print('Saving intermediate data')
                    del data
                    np.save('temp/data_best',data_best, allow_pickle=True, fix_imports=True)

                    mean_lucky_Stracking = np.sum(data_best,2)/n
                    np.save('temp/mean_lucky_Stracking',mean_lucky_Stracking, allow_pickle=True, fix_imports=True)

                    print('Done')
                    
                    et = time.time()
                    # get the execution time
                    elapsed_time = et - st
                    print('Execution time:', elapsed_time, 'seconds')

                    plt.figure()
                    plt.imshow(mean_lucky_Stracking)
                    text_temp = 'Mean {r:.2f}% speckle tracked frames'
                    plt.title(text_temp.format(r=100*r))

                    mean_lucky_Stracking_norm = fnc.normalise(mean_lucky_Stracking,12)
                    # max = np.unravel_index(mean_lucky_Stracking_norm.argmax(), mean_lucky_Stracking_norm.shape)
                    # plt.figure()
                    # plt.plot(mean_lucky_Stracking_norm[max[0],:])
                    # plt.grid()
                    # text_temp = 'Vertical cut'
                    # plt.title(text_temp)
                    # plt.figure()
                    # plt.plot(mean_lucky_Stracking_norm[:,max[1]])
                    # plt.grid()
                    # text_temp = 'Horizontal cut'
                    # plt.title(text_temp)

                    #---------------------3D plot
                    # fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
                    # # Make data.
                    # widht, height = mean_lucky_Stracking.shape
                    # X = np.arange(-widht//2, widht//2,1)
                    # Y = np.arange(-height//2, height//2,1)
                    # X, Y = np.meshgrid(X, Y)
                    # R = np.sqrt(X**2 + Y**2)
                    # Z = np.sin(R)
                    # # Plot the surface.
                    # surf = ax.plot_surface(X, Y, mean_lucky_Stracking_norm, cmap=cm.coolwarm, linewidth=0, antialiased=False)
                    # # Customize the z axis.
                    # ax.set_zlim(0, 1)
                    # ax.zaxis.set_major_locator(LinearLocator(10))
                    # ax.zaxis.set_major_formatter('{x:.02f}')
                    # fig.colorbar(surf, shrink=0.5, aspect=5)
                    # text_temp = 'Mean {r:.2f}% speckle tracked frames'
                    # plt.title(text_temp.format(r=100*r))
                    #----------------------


                    mean_lucky_Stracking_rot = ndimage.rotate(mean_lucky_Stracking, 35,reshape=False)
                  


                    mean_lucky_Stracking_rot_norm = fnc.normalise(mean_lucky_Stracking_rot,0)
                    mean_lucky_Stracking_rot_norm = fnc.normalise(mean_lucky_Stracking_rot_norm,0)

                    plt.figure()
                    plt.imshow(mean_lucky_Stracking_rot_norm)
                    
                    max = np.unravel_index(mean_lucky_Stracking_rot_norm.argmax(), mean_lucky_Stracking_rot_norm.shape)

                    cut_1 = mean_lucky_Stracking_rot_norm[max[0],:]
                    cut_1 = ((cut_1-np.min(cut_1))/(np.max(cut_1)-np.min(cut_1)))
                    cut_2 = mean_lucky_Stracking_rot_norm[:,max[1]]
                    cut_2 = ((cut_2-np.min(cut_2))/(np.max(cut_2)-np.min(cut_2)))

                    plt.figure()
                    plt.plot(cut_1)
                    plt.grid()
                    text_temp = 'Vertical cut'
                    plt.title(text_temp)
                    plt.figure()
                    plt.plot(cut_2)
                    plt.grid()
                    text_temp = 'Horizontal cut'
                    plt.title(text_temp)                                     
                    plt.show()
                    eError=eErrors.E_end_programm   

                case 12:
                    print('display')
                    mean_raw = np.load('temp/mean_image_raw.npy')
                    mean_bcg = np.load('temp/mean_image_bcg.npy')
                    mean_dustfree = np.load('temp/mean_image_dustfree.npy')
                    mean_image_Stracking = np.load('temp/mean_image_Stracking.npy')
                    mean_lucky_Stracking = np.load('temp/mean_lucky_Stracking.npy')
                    
                    # plt.figure()
                    # plt.imshow(mean_raw)
                    # plt.title('Mean raw data')
                    # plt.figure()
                    # plt.imshow(mean_bcg)
                    # plt.title('Mean Bcg removed data')
                    # plt.figure()
                    # plt.imshow(mean_dustfree)
                    # plt.title('Mean dustfree data')
                    # plt.figure()
                    # plt.imshow(mean_image_Stracking)
                    # plt.title('Mean speckle tracked data')
                    # plt.figure()
                    # plt.imshow(mean_lucky_Stracking)
                    # plt.title('Mean speckle tracked data (20% best)')

                    

                    plt.show()


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

lucky_process(a,b,r,eError)