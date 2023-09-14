# Lecture d'image en lucky imaging.
# 
# Inspired and somtime simply translated form process2.pro (L.J)
#
# 2023.06.26    -   Dominique Humbert   -   Version initiale

from astropy.io import fits
import numpy as np
from customClasses import eErrors
import matplotlib.pyplot as plt
import os
import time
from customClasses import eErrors
from customClasses import roi
#------------------LIBRARIES
from astropy import stats
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from scipy import ndimage


global folder; folder = 'C:/Users/ADM/OneDrive - HESSO/Dominique/07_mesures/08_aberations_with_foen'
global filename; filename = '/AF_'

global Roi_focus; Roi_focus = roi(485, 995, 64, 64)                # center x, center y, x dist from center, y dist from center
global Roi_defocus; Roi_defocus = roi(500, 355, 256, 256)          # center x, center y, x dist from center, y dist from center

global noise_filter_param; noise_filter_param = 0.2
global NdustSpots; NdustSpots = 0

'''------------------------------------------------------------------------------------
    # eError_handling(code)
    Error handling routine and help
    ## in:
      code: errore code
    ## out:
        -    ------------------------------------------------------------------------------------'''
def eError_handling(code=eErrors.E_arg_error):
    match code:
        case eErrors.E_all_fine:
            print("Programm terminated with success")
        case eErrors.E_arg_error:
            print("--------------------------")
            print("Too many arguments")
            print("arg1 = process\n\
    0: Read and save data\n\
    1: Whole process\n\
    2: ROI and the rest\n\
    3: Sigma filter and the rest\n\
    4: Normalisation\n\
    5: Background removal (if syntBCG == True)\n\
    6: Dark/flat frame removal (if syntBCG == False)\n\
    7: Noise cut (if implemented)\n\
    8: Gausian noise filtering (if implemented)\n\
    9: Dust shadow removal\n\
    10: tracking\n\
    11: Selection\n\
    12: plots")
            print("arg2 = qty of frames to use")
            print("arg3 = fraction of frames to stack (ie. 0.2)")
            print("arg4 = synthetic background mode")
            print("arg5 = error code. insert 0 to run the function directly")
            print("--------------------------")

def lucky_process_focus(a=1, b=None, r=0.1, synt_BCG = False, eError=eErrors.E_all_fine):
    '''------------------------------------------------------------------------------------
    # main(a, b, eError=eErrors.E_arg_error)
    Processes lucky imaging frames
    ## in:
      a: number of images
      b: starting point
      r: ratio to stack
      synt_BCG: synthetic background mode
      eError: Error status
    ## out:
      -
    ------------------------------------------------------------------------------------'''
    if help:
        eError == eErrors.E_arg_error
    
    if (eError == eErrors.E_all_fine):

        st = time.time()
        # ------------------------------Variables------------------------------
        data = None
        path = folder+filename
        temp_path = folder+'/temp'
        [k, j] = fits.getdata(path+'001.fits', ext=0).shape
        # -----------------------------count files-----------------------------

        n = -1  # there also is a folder that should not be counted

        # Iterate directory
        for path in os.listdir(folder):
            # check if current path is a file
            if os.path.isfile(os.path.join(folder, path)):
                n += 1
                
        # ----------------------------------------------------------
        if b is None:       # qty of frames to use
            b = 20
        else:
            n = b
        if r is None:       # qty of frames to use
            r = 0.05

        # ----------------------------------------------------------

        path = folder+filename
        
        # -------------------------------program--------------------------------
        while eError==eErrors.E_all_fine:
            match a:
                case 0: # only read and save raw data
                    print('Read data')
                    data = readData(n,path)

                    mean_image_raw = np.sum(data,2)/n
                    print('Saving intermediate data')
                    temp_name = temp_path+'data_raw'
                    np.save(temp_name,data, allow_pickle=True, fix_imports=True)
                    print('Done')  
                    eError=eErrors.E_end_programm
                case 1:
                    print('Read data')
                    data = readData(n,path)
                    mean_image_raw = np.sum(data,2)/n
                    print('Saving intermediate data')
                    temp_name = temp_path+'mean_image_raw'
                    np.save(temp_name,mean_image_raw, allow_pickle=True, fix_imports=True)
                    print('Done')
                    # plt.figure()
                    # plt.imshow(data[:,:,0],vmin=0,vmax=4095)
                    # plt.title('Raw image')
                    # plt.figure()
                    
                    # plt.imshow(mean_image_raw**(1/4))
                    # text_temp = 'Raw mean on {n:d} samples'
                    # plt.title(text_temp.format(n=n))
                    # plt.show()
                    # eError=eErrors.E_end_programm
                case 2:
                    print('ROI')
                    if data is None:
                        print('Loading raw data')
                        temp_name = temp_path+'data_raw.npy'
                        data = np.load(temp_name)


                    data = data[Roi_focus.x-Roi_focus.dx:Roi_focus.x+Roi_focus.dx,Roi_focus.y-Roi_focus.dy:Roi_focus.y+Roi_focus.dy,:]
                    mean_image_roi = np.sum(data,2)/n

                    print('Saving intermediate data')
                    temp_name = temp_path+'data_roi'
                    np.save(temp_name,data, allow_pickle=True, fix_imports=True)
                    temp_name = temp_path+'mean_image_roi.fits'
                    hdu = fits.PrimaryHDU(mean_image_roi)
                    hdu.writeto(temp_name,overwrite=True)
                    print('Done')

                    # plt.figure()
                    # plt.imshow(data[:,:,0],vmin=0,vmax=4095)
                    # plt.title('ROI image')

                    # # plt.figure()
                    # # plt.imshow(mean_image_roi)
                    # # text_temp = 'ROI mean on {n:d} samples'
                    # # plt.title(text_temp.format(n=n))
                    # plt.show()
                    # eError = eErrors.E_end_programm
                case 3:
                    print('Sigma filter')
                    if data is None:
                        temp_name = temp_path+'data_roi.npy'
                        print('Loading ROI data')
                        data = np.load(temp_name)
                    
                    status = 0
                    for i in range(0,n):
                        data[:,:,i] = stats.sigma_clip(data[:,:,i],maxiters=None)
                        if 100*i//n == status*10:
                            print('Sigma filtering: ',status*10,'% done')
                            status += 1
                    temp_name = temp_path+'data_sigma'
                    np.save(temp_name,data, allow_pickle=True, fix_imports=True)
                    print('Done')
                    mean_image_sigma = np.sum(data,2)/n
                    temp_name = temp_path+'mean_image_sigma.fits'
                    hdu = fits.PrimaryHDU(mean_image_sigma)
                    hdu.writeto(temp_name,overwrite=True)
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
                        temp_name = temp_path+'data_sigma.npy'
                        data = np.load(temp_name)

                    status = 0
                    for i in range(0,n):
                        data[:,:,i] = normalise(data[:,:,i],12)
                        if 100*i//n == status*10:
                            print('Normalisation: ',status*10,'% done')
                            status += 1

                    print('Saving intermediate data')
                    temp_name = temp_path+'data_norm'
                    np.save(temp_name,data, allow_pickle=True, fix_imports=True)
                    print('Done')
                    mean_image_norm = np.sum(data,2)/n
                    # plt.figure()
                    # plt.imshow(data[:,:,0])
                    # plt.show()

                    # eError = eErrors.E_end_programm

                
                case 5:
                    if synt_BCG:
                        print("Background")
                        if data is None:
                            print('loading norm data')
                            temp_name = temp_path+'data_norm.npy'
                            data = np.load(temp_name)
                        # plt.figure()
                        # plt.imshow(data[:,:,0])
                        # plt.title('Data')
                        # plt.show()    
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
                            mask[:,:,i] = buid_mask(data[:,:,i],mask_radius,[grid_x_UINT, grid_y_UINT],speckles[:,i])
                            model[:,:,i], poly_error_flag = polynomial_mask(data[:,:,i],mask[:,:,i],[grid_x, grid_y],4)
                            if poly_error_flag:
                                data[:,:,i] -= data[:,:,i]            # if the model could not be build, the data is set to 0. thus not beeing used 
                            else:
                                data[:,:,i] = data[:,:,i] - model[:,:,i]

                            if 100*i//n == status*10:
                                print('background removal: ',status*10,'% done')
                                status += 1

                        print('Saving intermediate data')
                        temp_name = temp_path+'data_mask'
                        np.save(temp_name,mask, allow_pickle=True, fix_imports=True)
                        # plt.figure()
                        # plt.imshow(mask[:,:,0])
                        # plt.title('Mask')
                        # plt.show()
                        del mask
                        temp_name = temp_path+'data_model'
                        np.save(temp_name,model, allow_pickle=True, fix_imports=True)
                        del model
                        temp_name = temp_path+'data_bcg'
                        np.save(temp_name,data, allow_pickle=True, fix_imports=True)
                        print('Done')

                        # status = 0
                        # for i in range(0,n):
                        #     data[:,:,i] = (data[:,:,i]-np.min(data[:,:,i]))/(np.max(data[:,:,i])-np.min(data[:,:,i]))
                        #     if 100*i//n == status*10:
                        #         print('Normalisation: ',status*10,'% done')
                        #         status += 1
                        status = 0
                        for i in range(0,n):
                            data[:,:,i] = normalise(data[:,:,i],12)
                            # print(np.min(data[:,:,i]))
                            if 100*i//n == status*10:
                                print('Normalisation: ',status*10,'% done')
                                status += 1
                        
                        mean_image_bcg = np.sum(data,2)/n 
                        temp_name = temp_path+'mean_image_bcg.fits'
                        hdu = fits.PrimaryHDU(mean_image_bcg)
                        hdu.writeto(temp_name,overwrite=True)
                        print('Done')

                        
                        # plt.figure()
                        # plt.imshow(mean_image_bcg,cmap='gray',vmin=0,vmax=4095)
                        # text_temp = 'Bcg removed mean on {n:d} samples'
                        # plt.title(text_temp.format(n=n))
                        # plt.show()
                        # eError = eErrors.E_end_programm
                case 6:
                     if not synt_BCG:
                        print("Dark/flat removal")
                        if data is None:
                            print('loading norm data')
                            temp_name = temp_path + 'data_norm.npy'
                            data = np.load(temp_name)
                            print('Done')
                        print('loading Dark/flat frames')

                        temp_name = temp_path + '/flat_norm.fits'
                        flat = fits.getdata(temp_name, ext=0)
                        flat = flat[Roi_focus.x-Roi_focus.dx:Roi_focus.x+Roi_focus.dx,Roi_focus.y-Roi_focus.dy:Roi_focus.y+Roi_focus.dy]
                        temp_name = temp_path + '/dark.fits'
                        dark = fits.getdata(temp_name, ext=0)
                        dark = dark[Roi_focus.x-Roi_focus.dx:Roi_focus.x+Roi_focus.dx,Roi_focus.y-Roi_focus.dy:Roi_focus.y+Roi_focus.dy]
                        print('Done')

                        status = 0
                        for i in range(0,n):
                            data[:,:,i] = (data[:,:,i]-dark)/flat
                            data[:,:,i] = normalise(data[:,:,i], 12)
                            if 100*i//n == status*10:
                                print('Normalisation: ',status*10,'% done')
                                status += 1
                case 7: 
                    print('Noise cut')
                    # if data is None:
                    #     print('Loading bcg data')
                    #     temp_name = temp_path+'data_bcg.npy'
                    #     data = np.load(temp_name)
                    #     k, j, n = data.shape
                    #     data = np.zeros((k,j,n))
                    # data[:,:,0] = noise_cut(data[:,:,0])
                    # plt.figure()
                    # plt.imshow(data[:,:,0])
                    # text_temp = 'Noise cut'
                    # plt.title(text_temp)
                    #plt.show()
   
                case 8:
                    print('Noise filtering')
                    if data is None:
                        print('Loading filtered data')
                        temp_name = temp_path+'data_bcg.npy'
                        data = np.load(temp_name)
                        print('Done')
                    height, width, n = data.shape
                    x = np.linspace(-1,1,height)
                    y = np.linspace(-1,1,width)
                    grid_x, grid_y = np.meshgrid(x, y)

                    # plt.figure()
                    # plt.imshow(data[:,:,0]**(1/8))
                    # text_temp = 'Unfiltered'
                    # plt.title(text_temp)

                    status = 0
                    for i in range(0,n):
                        data[:,:,i] = noise_filter(data[:,:,i],noise_filter_param,[grid_x, grid_y],0)
                        if 100*i//n == status*10:
                            print('Noise filtering: ',status*10,'% done')
                            status += 1
                    print('Saving intermediate data')
                    temp_name = temp_path+'noise_filter_param'
                    np.save(temp_name,noise_filter_param, allow_pickle=True, fix_imports=True)
                    temp_name = temp_path+'data_filtered'
                    np.save(temp_name,data, allow_pickle=True, fix_imports=True)

                    mean_image_filtered = np.sum(data,2)/n 
                    temp_name = temp_path+'mean_image_filtered.fits'
                    hdu = fits.PrimaryHDU(mean_image_filtered)
                    hdu.writeto(temp_name,overwrite=True)
                    print('Done')
                    # plt.figure()
                    # plt.imshow(data[:,:,0]**(1/1),)
                    # text_temp = 'Noise filtered'
                    # plt.title(text_temp)
                    # plt.show()
                    # eError = eErrors.E_end_programm
                case 9:
                    print('Dust spots')
                    for k in range(0,NdustSpots):
                        spot_model_radius = 25  
                        status = 0
                        if data is None:
                            print('Loading filtered data')
                            temp_name = temp_path+'data_filtered.npy'
                            data = np.load(temp_name)

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
                            data[:,:,i] = image_dusting(data[:,:,i],dust_model,4)
                            if 100*i//n == status*10:
                                print('Dust shadow removal: ',status*10,'% done')
                                status += 1
                        
                        print('Saving intermediate data')
                        temp_name = temp_path+'data_dustfree'
                        np.save(temp_name,data, allow_pickle=True, fix_imports=True)

                        # mean_image_dustfree = np.sum(data,2)/n 
                        # temp_name = temp_path+'mean_image_dustfree.fits'
                        # hdu = fits.PrimaryHDU(mean_image_dustfree)
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
                        temp_name = temp_path+'data_dustfree.npy'
                        data = np.load(temp_name) 
                    status = 0
                    Speckle_pos = np.zeros((n,2))

                    for i in range(0, n):
                        speckles = np.unravel_index(data[:,:,i].argmax(), data[:,:,i].shape)
                        if i == 0:
                            ref = speckles
                        Dx = ref[0]-speckles[0] 
                        Dy = ref[1]-speckles[1]
                        Speckle_pos[i,0] = speckles[0] 
                        Speckle_pos[i,1] = speckles[1] 
                        data[:,:,i] = np.roll(data[:,:,i],(Dy,Dx),axis=(1,0))
                        if 100*i//n == status*10:
                            print('Tracking: ',status*10,'% done')
                            status += 1
                    print('Saving intermediate data')
                    temp_name = temp_path+'Speckle_pos'
                    np.save(temp_name,Speckle_pos, allow_pickle=True, fix_imports=True)

                    temp_name = temp_path+'data_Stracked'
                    np.save(temp_name,data, allow_pickle=True, fix_imports=True)

                    mean_image_Stracking = np.sum(data,2)/n
                    temp_name = temp_path+'mean_image_Stracking.fits'
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
                        temp_name = temp_path+'data_Stracked.npy'
                        data = np.load(temp_name) 

                    #r = 0.2 # [%] # now a global variable
                    text_temp = 'Selecting the {r:.2f}% best frames'
                    print(text_temp.format(r=100*r))

                    # Find the most luminous pixel on each frame
                    max_pixel_light = np.zeros((n))
                    for i in range(0,n):
                        max_pixel_light[i] = np.amax(data[:,:,i])
                    order = np.argsort(max_pixel_light)
                    #print(order)

                    # sort by Luminosity  in the brightest pixel 
                    n_needed = r*n
                    n_needed = int(-(-n_needed//1))   # Ceil without math module. We need a finite number of frame
                    selected_indices = order[n-n_needed:n]

                    print('saving intermediate data')
                    temp_name = temp_path+'selected_indices'
                    np.save(temp_name,selected_indices, allow_pickle=True, fix_imports=True)
                    #print(selected_indices)
                    print('done')

                    # Select only the best
                    data_best  = data[:,:,selected_indices]
                    del data

                    mean_lucky_Stracking = np.sum(data_best,2)/n
                    mean_lucky_Stracking = normalise(mean_lucky_Stracking, 12)

                    print('saving intermediate data')
                    hdu = fits.PrimaryHDU(mean_lucky_Stracking)
                    temp_name = temp_path+'(mean_lucky_Strackings'
                    hdu.writeto(temp_name,overwrite=True)
                    print('done')
                  
                    plt.figure()
                    plt.imshow(mean_lucky_Stracking)
                    text_temp = 'Mean {r:.2f}% speckle tracked frames'
                    plt.title(text_temp.format(r=100*r))


                    # mean_lucky_Stracking_norm = normalise(mean_lucky_Stracking,12)
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



                    # mean_lucky_Stracking_rot_norm = normalise(mean_lucky_Stracking_rot,0)
                    # #mean_lucky_Stracking_rot_norm = normalise(mean_lucky_Stracking_rot_norm,0)

                    # temp_name = temp_path+'mean_lucky_Stracking_rot_norm'
                    # hdu = fits.PrimaryHDU(mean_lucky_Stracking_rot_norm)
                    # hdu.writeto(temp_name,overwrite=True)

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
                    # mean_raw = np.load('temp/mean_image_raw.npy')
                    # mean_bcg = np.load('temp/mean_image_bcg.npy')
                    # mean_dustfree = np.load('temp/mean_image_dustfree.npy')
                    # mean_image_Stracking = np.load('temp/mean_image_Stracking.npy')

                    #data_roi = np.load('temp/data_roi.npy')
                    #mean_lucky_Stracking = np.load('temp/mean_lucky_Stracking.npy')
                    #plt.imsave('roi.jpg',data_roi[:,:,0])


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

    eError_handling(eError)

def lucky_process_defocus(a = 1, synt_BCG = False, eError=eErrors.E_all_fine):
    '''------------------------------------------------------------------------------------
    # main(a, b, eError=eErrors.E_arg_error)
    Processes lucky imaging frames
    ## in:
      a: number of images
      r: ratio to stack
      synt_BCG: synthetic background mode
      eError: Error status
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
                        
                    data = data[Roi_defocus.x-Roi_defocus.dx:Roi_defocus.x+Roi_defocus.dx,Roi_defocus.y-Roi_defocus.dy:Roi_defocus.y+Roi_defocus.dy,:]
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
                        data[:,:,i] = normalise(data[:,:,i],12)
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
                            mask[:,:,i] = buid_mask(data[:,:,i],mask_radius,[grid_x_UINT, grid_y_UINT]) # mask is buit in the center
                            model[:,:,i], poly_error_flag = polynomial_mask(data[:,:,i],mask[:,:,i],[grid_x, grid_y],4)
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
                            data[:,:,i] = normalise(data[:,:,i],12)
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
                        temp_name = temp_path + '/flat_norm.fits'
                        flat = fits.getdata(temp_name, ext=0)
                        flat = flat[Roi_defocus.x-Roi_defocus.dx:Roi_defocus.x+Roi_defocus.dx,Roi_defocus.y-Roi_defocus.dy:Roi_defocus.y+Roi_defocus.dy]
                        temp_name = temp_path + '/dark.fits'
                        dark = fits.getdata(temp_name, ext=0)
                        dark = dark[Roi_defocus.x-Roi_defocus.dx:Roi_defocus.x+Roi_defocus.dx,Roi_defocus.y-Roi_defocus.dy:Roi_defocus.y+Roi_defocus.dy]
                        print('Done')

                        for i in range(0,n):
                            data[:,:,i] = (data[:,:,i]-dark)/flat
                            data[:,:,i] = normalise(data[:,:,i], 12)
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
                        data[:,:,i] = noise_filter(data[:,:,i],0.2,[grid_x, grid_y],0)
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
                    for k in range(0,NdustSpots):
                        print('Dust spots')
                        spot_model_radius = 25  
                        status = 0
                        if data is None:
                            print('Loading filtered data')
                            temp_name = temp_path+'data_filtered_def.npy'
                            data = np.load(temp_name)

                        try: mean_image_bcg
                        except NameError: mean_image_bcg = np.load('temp/mean_image_bcg_def.npy')
                        # plt.figure()
                        # plt.imshow(mean_image_bcg)
                        # text_temp = 'Shadows everywhere!!!'
                        # plt.title(text_temp)
                        # find a dust spot for the model
                        index_model = np.unravel_index(mean_image_bcg[:,:].argmin(), mean_image_bcg[:,:].shape)
                        dust_model = mean_image_bcg[index_model[0]-spot_model_radius:index_model[0]+spot_model_radius,index_model[1]-spot_model_radius:index_model[1]+spot_model_radius]
                        dust_model /= np.min(dust_model)
                        
                        for i in range(0,n):
                            data[:,:,i] = image_dusting(data[:,:,i],dust_model,4)
                            if 100*i//n == status*10:
                                print('Dust shadow removal: ',status*10,'% done')
                                status += 1
                        
                        print('Saving intermediate data')
                        temp_name = temp_path+'data_dustfree_def'
                        np.save(temp_name,data, allow_pickle=True, fix_imports=True)

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
                    print('loading some data')
                    temp_name = temp_path+'Speckle_pos.npy'
                    speckles = np.load(temp_name)
                    status = 0
                    speckles = speckles[indices].astype(np.uint16)

                    for i in range(0, n):
                        Dx = speckles[0,0]-speckles[i,0] 
                        Dy = speckles[0,1]-speckles[i,1] 
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
                    
                    mean_lucky_Stracking = normalise(mean_lucky_Stracking, 12)
                    print('saving intermediate data')
                    temp_name = temp_path+'mean_lucky_Stracking_def.fits'
                    hdu = fits.PrimaryHDU(mean_lucky_Stracking)
                    hdu.writeto(temp_name,overwrite=True)
                    print('done')


                    plt.figure()
                    plt.imshow(mean_lucky_Stracking)
                    text_temp = 'Mean lucky Stracking_def'
                    plt.title(text_temp)

                    # mean_lucky_Stracking_norm = normalise(mean_lucky_Stracking,12)
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



                    # mean_lucky_Stracking_rot_norm = normalise(mean_lucky_Stracking_rot,0)
                    # #mean_lucky_Stracking_rot_norm = normalise(mean_lucky_Stracking_rot_norm,0)

                    # temp_name = temp_path+'mean_lucky_Stracking_rot_norm_def.npy'
                    # hdu = fits.PrimaryHDU(mean_lucky_Stracking_rot_norm)
                    # hdu.writeto(temp_name,overwrite=True)

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
                    # plt.show()
                    # eError=eErrors.E_end_programm   

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

    eError_handling(eError)

def readData(n, fname):
    '''------------------------------------------------------------------------------------
    # readData(n, fname)
    Reads n .fits files with the path fnameXXX.fits (3 or more X's)
    ## in:
      n: number of images
      fname: filename
    ## out:
      data: 3D array with pixel data
    ------------------------------------------------------------------------------------'''
    status = 0
    for i in range(n):
        filename = fname+'{:03d}'.format(i+1)+'.fits'
        if i==0:                    # Gets the pic size on first pass
            [k, j] = fits.getdata(fname+'001.fits', ext=0).shape
            data = np.zeros((k,j,n))
        data[:,:,i] = fits.getdata(filename, ext=0)

        if 100*(i+1)//n == status*10:
            print('Reading data: ',status*10,'% done')
            status += 1
    return data


def normalise(data, bits_depth = 12):
    '''------------------------------------------------------------------------------------
    # normalize(data,bits_depth)

    Normalise images 

    ## in:
        data : data: 2D array with pixel data
        bit_depth : 2^bits resolution (1<=bits<=16)
    ## out:
        data: 3D array with pixel data
    ------------------------------------------------------------------------------------'''
    

    if bits_depth<0 or bits_depth>16:
        bits_depth = 12

    #print('Normalising images in ',bits,'bits')
    if bits_depth == 0:
        data = ((data - np.min(data)) / (np.max(data) - np.min(data))).astype(np.double) 
    else:
        data = (data - np.min(data)) / (np.max(data) - np.min(data)).astype(np.double) 
        data = (data * 2**bits_depth).astype(np.uint16)

    return data


def buid_mask(data, radius, grid, center=[None,None]):
    '''#------------------------------------------------------------------------------------
        buid_maks(data,center=None)

        in:
            data : Array the size of the image for whoch the mask is built
            radius : mask radius
            center (optional) : x,y coordinate of the mask's center 

        out:
            mask : boolean array with True in the mask circle
        ------------------------------------------------------------------------------------
    '''

    if not (center[0] and center[1]):
        center = np.zeros((2))
        height, width = data.shape
        center[0] = width // 2
        center[1] = height // 2
                
    distance_from_center = np.sqrt((grid[0]- center[0])**2 + (grid[1] - center[1])**2)
    # Créer le masque pour exclure les pixels à l'intérieur du cercle
    mask = distance_from_center <= radius

    return mask



def polynomial_mask(data, mask, grid, order=4):
    '''------------------------------------------------------------------------------------
        polynomial_mask(data, mask_radius,order=4)
        according to "Analyse de données"

        in:
          data : image to interpolate from
          mask : mask boolean mask. True means to be removed
          Order : Order of the appoximation. Default = 4
        
        out: 
          approximated image
          error flag
         
        to do:
        ------------------------------------------------------------------------------------'''
    error_flag = False
    nmodes = int(((order+1)**2-(order+1))/2+(order+1))
    mask_dbl = np.zeros(mask.shape)
    mask_dbl[mask==False] = 1.   # met des 1. là ou le masque n'est pas présent
    height, width = data.shape
    model = np.zeros(mask.shape)
    basis = np.zeros((height,width,nmodes))


    j_mode = -1
    for i in range(order+1):
        for j in range(i+1):
            j_mode +=1
            basis[:,:,j_mode] =grid[0]**(i-j)*grid[1]**j
    A = np.zeros((nmodes,nmodes)).astype(np.double)

    for i in range(0,nmodes):
        for j in range(0,i+1): 
            A[i,j] = np.sum(basis[:,:,i]*basis[:,:,j]*mask_dbl)

    b = np.zeros(nmodes)
    for i in range(nmodes):
         b[i]=np.sum(data*basis[:,:,i]*mask_dbl)

    try:
        a = np.dot(np.linalg.inv(A),b)
    except np.linalg.LinAlgError:
        error_flag = True
        print('Linear algorithm error, model not build')

    if not error_flag:
        model = np.zeros((height,width))
        for i in range(0,nmodes):
            model += a[i]*basis[:,:,i] 

    return model, error_flag


def noise_cut(data, threshold):
    '''------------------------------------------------------------------------------------
        noise_cut(data, threshold=None)
        Noise cut filter

        in:
          data : image to interpolate from. square and dimmentions are ^2 (padding not implemented)
          threshold : cutting threshold
        
        out: 
          approximated image
         
        to do:
            faire fonctionner
        ------------------------------------------------------------------------------------'''

    # Perform Fourier transform
    data_tf = np.fft.fftshift(np.fft.fft2(data))
    
    # Apply filter in Fourier domain
    mask = np.abs(data_tf) > threshold
    data_tf_filtered = data_tf*mask
    
    # Inverse Fourier transform
    data = np.fft.ifft2(np.fft.ifftshift(data_tf_filtered)).real

    
    return data.astype(np.double)

def noise_filter(data, sigma, grid, power = 2):
    '''------------------------------------------------------------------------------------
        noise_filtering(data, sigma=1)
        gaussin noise filtering

        in:
          data : image to interpolate from. square and dimmentions are ^2 (padding not implemented)
          sigma : 
        
        out: 
          approximated image
         
        to do:
        ------------------------------------------------------------------------------------'''

    data_tf_filtered = np.fft.fftshift(np.fft.fft2(data))*np.exp(-(((grid[0]**2+grid[1]**2)/(sigma))**power))
        
        # Inverse Fourier transform
    data = np.fft.ifft2(np.fft.ifftshift(data_tf_filtered)).real


    return data.astype(np.double)



def image_dusting(data, model, NdustSpots=1):
    '''------------------------------------------------------------------------------------
        image_dusting(data, sigma=1)
        gaussin noise filtering

        in:
          data: image to remove shdows from.
          model: shadow model
          NdustSpots: Number of spots to remove
        
        out: 
          shadowless image
         
        to do:
        ------------------------------------------------------------------------------------'''
    height,width = data.shape
    spot_model_radius, trash = model.shape
    spot_model_radius //= 2

    bigger = np.zeros((height+2*spot_model_radius, width+2*spot_model_radius))  # Allows for removal of shadows that are not entierly on the image
    index = np.zeros(2)
    for i in range(NdustSpots):
                        
        index = np.unravel_index(data.argmin(), data.shape)
        aa = index[0]-spot_model_radius
        bb = index[0]+spot_model_radius
        cc = index[1]-spot_model_radius
        dd = index[1]+spot_model_radius
                        
        bigger[spot_model_radius:height+spot_model_radius,spot_model_radius:width+spot_model_radius] = data[:,:]
        bigger[index[0]:index[0]+2*spot_model_radius,index[1]:index[1]+2*spot_model_radius] -= model*data[index[0],index[1]]
        data[:,:] = bigger[spot_model_radius:height+spot_model_radius,spot_model_radius:width+spot_model_radius]

    return data
