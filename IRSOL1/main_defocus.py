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
global folder; folder = 'C:/Users/ADM/OneDrive - HESSO/Dominique/07_mesures/07_aberation_disk'
global filename; filename = '/AD_'
global ROIx; ROIx = 478
global ROIy; ROIy = 400
global ROIdx; ROIdx = 256
global ROIdy; ROIdy = 256
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


def lucky_process_defocuss(a, synt_BCG = False, eError=eErrors.E_arg_error):
    '''------------------------------------------------------------------------------------
    # main(a, b, eError=eErrors.E_arg_error)
    Processes lucky imaging frames
    ## in:
      a: number of images
      b: starting point
      r: ratio to stack
    ## out:
      -
    ------------------------------------------------------------------------------------'''
    if (eError == eErrors.E_all_fine):

        st = time.time()
        # ------------------------------Variables------------------------------
        data = None
        path = folder+filename
        temp_path = folder+'/temp/'
        [k, j] = fits.getdata(path+'001.fits', ext=0).shape

        path = folder+filename
        
        # -------------------------------program--------------------------------
        while eError==eErrors.E_all_fine:
            match a:
                case 0: # only read and save raw data
                    1
                case 1:
                    print('Read data')
                    temp_name = temp_path+'selected_indices.npy'
                    indices = np.load(temp_name)
                    status = 0
                    n = np.prod(indices.shape)
                    for i in range(0, n):
                        Imfilename = path+'{:03d}'.format(indices[i])+'.fits'
                        if i==0:                    # Gets the pic size on first pass
                            [k, j] = fits.getdata(path+'001.fits', ext=0).shape
                            data = np.zeros((k,j,n))
                        data[:,:,i] = fits.getdata(Imfilename, ext=0)

                        if 100*(i+1)//n == status*10:
                            print('Reading data: ',status*10,'% done')
                            status += 1

                    mean_image_raw = np.sum(data,2)/n
                    print('Saving intermediate data')
                    temp_name = temp_path+'mean_image_raw_def'
                    np.save(temp_name,mean_image_raw, allow_pickle=True, fix_imports=True)
                    print('Done')

                    # eError = eErrors.E_end_programm
                case 2:
                    print('ROI')
                    if data is None:
                        temp_name = temp_path+'data_roi_def.npy'
                        print('Loading ROI data')
                        data = np.load(temp_name)
                        
                    data = data[ROIx-ROIdx:ROIx+ROIdx,ROIy-ROIdy:ROIy+ROIdy,:]
                    mean_image_roi = np.sum(data,2)/n

                    print('Saving intermediate data')
                    temp_name = temp_path+'data_roi_def'
                    np.save(temp_name,data, allow_pickle=True, fix_imports=True)
                    temp_name = temp_path+'mean_image_roi_def.fits'
                    hdu = fits.PrimaryHDU(mean_image_roi)
                    hdu.writeto(temp_name,overwrite=True)
                    print('Done')

                    # plt.figure()
                    # plt.imshow(data[:,:,0])
                    # plt.show()
                    # eError = eErrors.E_end_programm
                case 3:
                    print('Sigma filter')
                    if data is None:
                        temp_name = temp_path+'data_roi_def.npy'
                        print('Loading ROI data')
                        data = np.load(temp_name)
                    
                    status = 0
                    for i in range(0,n):
                        data[:,:,i] = stats.sigma_clip(data[:,:,i],maxiters=None)
                        if 100*i//n == status*10:
                            print('Sigma filtering: ',status*10,'% done')
                            status += 1
                    temp_name = temp_path+'data_sigma_def'
                    np.save(temp_name,data, allow_pickle=True, fix_imports=True)
                    print('Done')
                    mean_image_sigma = np.sum(data,2)/n
                    temp_name = temp_path+'mean_image_sigma_def.fits'
                    hdu = fits.PrimaryHDU(mean_image_sigma)
                    hdu.writeto(temp_name,overwrite=True)
                    
                case 4:
                    print('Normalisation')
                    if data is None:
                        print('Loading sigma data')
                        temp_name = temp_path+'data_sigma_def.npy'
                        data = np.load(temp_name)

                    status = 0
                    for i in range(0,n):
                        data[:,:,i] = fnc.normalise(data[:,:,i],12)
                        if 100*i//n == status*10:
                            print('Normalisation: ',status*10,'% done')
                            status += 1

                    print('Saving intermediate data')
                    temp_name = temp_path+'data_norm_def'
                    np.save(temp_name,data, allow_pickle=True, fix_imports=True)
                    print('Done')
                    mean_image_norm = np.sum(data,2)/n
                
                case 5:
                    if synt_BCG:
                        print("Background")
                        if data is None:
                            print('loading norm data')
                            temp_name = temp_path+'data_norm_def.npy'
                            data = np.load(temp_name)
   
                        mask_radius = 35  # Rayon du cercle à exclure
                        
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

                        # -----------    Build grids
                        height, width = mask[:,:,0].shape
                        x = np.linspace(0, width-1, width).astype(np.uint16)
                        y = np.linspace(0, height-1, height).astype(np.uint16)
                        grid_x_UINT, grid_y_UINT = np.meshgrid(x, y)

                        x = np.linspace(-1,1,width)
                        y = np.linspace(-1,1,height)
                        grid_x, grid_y = np.meshgrid(x, y)
                        # -----------    Build grids


                        status = 0
                        for i in range(0,n):
                            mask[:,:,i] = fnc.buid_mask(data[:,:,i],mask_radius,[grid_x_UINT, grid_y_UINT],speckles[:,i])
                            model[:,:,i], poly_error_flag = fnc.polynomial_mask(data[:,:,i],mask[:,:,i],[grid_x, grid_y],4)
                            if poly_error_flag:
                                data[:,:,i] -= data[:,:,i]            # if the model could not be build, the data is set to 0. thus not beeing used 
                            else:
                                data[:,:,i] = data[:,:,i] - model[:,:,i]

                            if 100*i//n == status*10:
                                print('background removal: ',status*10,'% done')
                                status += 1

                        print('Saving intermediate data')
                        temp_name = temp_path+'data_mask_def'
                        np.save(temp_name,mask, allow_pickle=True, fix_imports=True)
                        # plt.figure()
                        # plt.imshow(mask[:,:,0])
                        # plt.title('Mask')
                        # plt.show()
                        del mask
                        temp_name = temp_path+'data_model_def'
                        np.save(temp_name,model, allow_pickle=True, fix_imports=True)
                        del model
                        temp_name = temp_path+'data_bcg_def'
                        np.save(temp_name,data, allow_pickle=True, fix_imports=True)
                        print('Done')
                        status = 0
                        for i in range(0,n):
                            data[:,:,i] = fnc.normalise(data[:,:,i],12)
                            print(np.min(data[:,:,i]))
                            if 100*i//n == status*10:
                                print('Normalisation: ',status*10,'% done')
                                status += 1
                        
                        mean_image_bcg = np.sum(data,2)/n 
                        temp_name = temp_path+'mean_image_bcg_def.fits'
                        hdu = fits.PrimaryHDU(mean_image_bcg)
                        hdu.writeto(temp_name,overwrite=True)
                        print('Done')


                case 6:
                     if not synt_BCG:
                        print("Dark/flat removal")
                        if data is None:
                            print('loading norm data')
                            temp_name = temp_path + 'data_norm_def.npy'
                            data = np.load(temp_name)
                            print('Done')
                        print('loading Dark/flat frames')
                        temp_name = folder + '/flat.fits'
                        flat = fits.getdata(temp_name, ext=0)
                        flat = flat[ROIx-ROIdx:ROIx+ROIdx,ROIy-ROIdy:ROIy+ROIdy]
                        temp_name = folder + '/dark.fits'
                        dark = fits.getdata(temp_name, ext=0)
                        dark = dark[ROIx-ROIdx:ROIx+ROIdx,ROIy-ROIdy:ROIy+ROIdy]
                        print('Done')

                        for i in range(0,n):
                            data[:,:,i] = (data[:,:,i]-dark)/flat
                            data[:,:,i] = fnc.normalise(data[:,:,i], 12)
                case 7:
                    1
                case 8:
                    print('Noise filtering')
                    if data is None:
                        print('Loading filtered data')
                        temp_name = temp_path+'data_bcg_def.npy'
                        data = np.load(temp_name)
                        print('Done')
                    height, width, n = data.shape
                    x = np.linspace(-1,1,height)
                    y = np.linspace(-1,1,width)
                    grid_x, grid_y = np.meshgrid(x, y)

                    status = 0
                    for i in range(0,n):
                        data[:,:,i] = fnc.noise_filter(data[:,:,i],0.2,[grid_x, grid_y],0)
                        if 100*i//n == status*10:
                            print('Noise filtering: ',status*10,'% done')
                            status += 1
                    print('Saving intermediate data')
                    temp_name = temp_path+'data_filtered_def'
                    np.save(temp_name,data, allow_pickle=True, fix_imports=True)

                    mean_image_filtered = np.sum(data,2)/n 
                    temp_name = temp_path+'mean_image_filtered_def.fits'
                    hdu = fits.PrimaryHDU(mean_image_filtered)
                    hdu.writeto(temp_name,overwrite=True)
                    print('Done')

                case 9:
                    print('Dust spots')
                    # NdustSpots = 7
                    # spot_model_radius = 25  
                    # status = 0;
                    # if data is None:
                    #     print('Loading filtered data')
                    #     temp_name = temp_path+'data_filtered_def.npy'
                    #     data = np.load(temp_name)

                    # try: mean_image_bcg
                    # except NameError: mean_image_bcg = np.load('temp/mean_image_bcg_def.npy')
                    # # plt.figure()
                    # # plt.imshow(mean_image_bcg)
                    # # text_temp = 'Shadows everywhere!!!'
                    # # plt.title(text_temp)
                    # # find a dust spot for the model
                    # index_model = np.unravel_index(mean_image_bcg[:,:].argmin(), mean_image_bcg[:,:].shape)
                    # dust_model = mean_image_bcg[index_model[0]-spot_model_radius:index_model[0]+spot_model_radius,index_model[1]-spot_model_radius:index_model[1]+spot_model_radius]
                    # dust_model /= np.min(dust_model)
                    
                    # for i in range(0,n):
                    #     data[:,:,i] = fnc.image_dusting(data[:,:,i],dust_model,4)
                    #     if 100*i//n == status*10:
                    #         print('Dust shadow removal: ',status*10,'% done')
                    #         status += 1
                    
                    # print('Saving intermediate data')
                    # temp_name = temp_path+'data_dustfree_def'
                    # np.save(temp_name,data, allow_pickle=True, fix_imports=True)

                    # mean_image_dustfree = np.sum(data,2)/n 
                    # temp_name = temp_path+'mean_image_dustfree_def.fits'
                    # hdu = fits.PrimaryHDU(mean_image_dustfree_def)
                    # hdu.writeto(temp_name,overwrite=True)
                    # print('Done')

                    # plt.figure()
                    # plt.imshow(mean_image_dustfree)
                    # text_temp = 'Dust shadow removed'
                    # plt.title(text_temp)
                    # plt.show()
                case 10:
                    print('tracking')  
                    if data is None:
                        print('Loading shadowless data')
                        temp_name = temp_path+'data_dustfree_def.npy'
                        data = np.load(temp_name) 
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
                    temp_name = temp_path+'data_Stracked_def'
                    np.save(temp_name,data, allow_pickle=True, fix_imports=True)

                    mean_image_Stracking = np.sum(data,2)/n
                    temp_name = temp_path+'mean_image_Stracking_def.fits'
                    hdu = fits.PrimaryHDU(mean_image_Stracking)
                    hdu.writeto(temp_name,overwrite=True)
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
                        temp_name = temp_path+'data_Stracked_def.npy'
                        data = np.load(temp_name) 

                    mean_lucky_Stracking = np.sum(data,2)/n
                    

                    print('saving intermediate data')
                    hdu = fits.PrimaryHDU(mean_lucky_Stracking)
                    hdu.writeto('Attempt1/temp/mean_lucky_Stracking_def.fits',overwrite=True)
                    print('done')
                  
                    print('saving intermediate data')
                    hdu = fits.PrimaryHDU(mean_lucky_Stracking)
                    hdu.writeto('Attempt1/temp/mean_lucky_Stracking_def.fits',overwrite=True)
                    print('done')

                    plt.figure()
                    plt.imshow(mean_lucky_Stracking)
                    text_temp = 'Mean lucky Strackng'
                    plt.title(text_temp)

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


                    # mean_lucky_Stracking_rot = ndimage.rotate(mean_lucky_Stracking, 35,reshape=False)



                    # mean_lucky_Stracking_rot_norm = fnc.normalise(mean_lucky_Stracking_rot,0)
                    # #mean_lucky_Stracking_rot_norm = fnc.normalise(mean_lucky_Stracking_rot_norm,0)

                    # hdu = fits.PrimaryHDU(mean_lucky_Stracking_rot_norm)
                    # hdu.writeto('Attempt1/temp/mean_lucky_Stracking_rot_norm.fits',overwrite=True)

                    # plt.figure()
                    # plt.imshow(mean_lucky_Stracking_rot_norm)
                    
                    # max = np.unravel_index(mean_lucky_Stracking_rot_norm.argmax(), mean_lucky_Stracking_rot_norm.shape)

                    # cut_1 = mean_lucky_Stracking_rot_norm[max[0],:]
                    # cut_1 = ((cut_1-np.min(cut_1))/(np.max(cut_1)-np.min(cut_1)))
                    # cut_2 = mean_lucky_Stracking_rot_norm[:,max[1]]
                    # cut_2 = ((cut_2-np.min(cut_2))/(np.max(cut_2)-np.min(cut_2)))

                    # plt.figure()
                    # plt.plot(cut_1)
                    # plt.grid()
                    # text_temp = 'Vertical cut'
                    # plt.title(text_temp)
                    # plt.figure()
                    # plt.plot(cut_2)
                    # plt.grid()
                    # text_temp = 'Horizontal cut'
                    # plt.title(text_temp)                                     
                    #plt.show()
                    #eError=eErrors.E_end_programm   

                    et = time.time()
                    # get the execution time
                    elapsed_time = et - st
                    print('Execution time:', elapsed_time, 'seconds')
                    

                case 12:
                    print('display')
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

lucky_process_defocuss(a,synt_BCG,eError)
