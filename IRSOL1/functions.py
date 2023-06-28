# Lecture d'image en lucky imaging.
# 
#
# 2023.06.26    -   Dominique Humbert   -   Version initiale
from astropy.io import fits
import numpy as np
from customClasses import eErrors
from customClasses import plotTypes

import matplotlib.pyplot as plt
from astropy import stats
from PIL import Image
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter

def eError_handling(code):
    match code:
        case eErrors.E_all_fine:
            print("Programm terminated with success")
        case eErrors.E_arg_error:
            print("--------------------------")
            print("Too many arguments")
            print("arg1 = process\n\
                  1:Read and save data\n\
                  2: Whole process\n\
                  3: ROI\n\
                  4: Sigma filter\n\
                  10: plot")
            print("arg2 = plot type\n\
                  raw : raw images\n\
                  roi: ROI normalised\
                  sigma: sigma filtered image\
                  mask: appoximation mask\
                  approx: picture with approximated mask\
                  no_bcg: picture with backgound removed")
            print("--------------------------")

#------------------------------------------------------------------------------------
# readData(n, fname)
# Reads n .fits files with the path fnameXXX.fits (3 or more X's)
# in:
#   n : number of images
#   fname : filename
#out:
#   data: 3D array with pixel data
#------------------------------------------------------------------------------------
def readData(n, fname):
    #data = fits.getdata(fname+'001.fits', ext=0)

    for i in range(n):
        filename = fname+'{:03d}'.format(i+1)+'.fits'
        if i==0:                    # Gets the pic size on firts pass
            [k, j] = fits.getdata(fname+'001.fits', ext=0).shape
            data = np.zeros((k,j,n))
        data[:,:,i] = fits.getdata(filename, ext=0)
        np.save('temp/data_raw',data, allow_pickle=True, fix_imports=True)
    return data


#------------------------------------------------------------------------------------
# ROI(n, data)
# Applies a 2-D region of interest do a 3-D array 
# in:
#   n : number of images
#   data : data: 3D array with pixel data
#out:
#   data: 3D array with pixel data
#------------------------------------------------------------------------------------
def ROI(n,data):
    if data is None:
        print('loading raw data')
        data = np.load('temp/data_raw.npy')
    data = data[500-295:500+294,625-295:625+294,:]     # Manually done

    for i in range(n):
        data[:,:,i] = (data[:,:,i] - np.min(data[:,:,i])) / (np.max(data[:,:,i]) - np.min(data[:,:,i]))
        data[:,:,i] = (data[:,:,i] * 2**12).astype(np.uint16)

    np.save('temp/data_roi',data, allow_pickle=True, fix_imports=True)
    return data


#------------------------------------------------------------------------------------
# sigma(n, data)
# in:
#   n : number of images
#   data : data: 3D array with pixel data
#out:
#   data: 3D array with pixel data
#------------------------------------------------------------------------------------
def sigma(n,data):
    if data is None:
        print('loading ROI data')
        data = np.load('temp/data_roi.npy')
    [k, j] = data[:,:,1].shape
    data_sigma = np.zeros((k,j,n))
    for i in range(n):
        data_sigma[:,:,i] = stats.sigma_clip(data[:,:,i],maxiters=None)         # A développer
    np.save('temp/data_sigma',data_sigma, allow_pickle=True, fix_imports=True)
    return data_sigma

#------------------------------------------------------------------------------------
# display_images(plotType)
# in:
#   plotType : omage to plot
#out:
#   -
#  
# to do:
#   add image selection
#------------------------------------------------------------------------------------
def display_images(plot_type):
    data_raw = None
    data_roi = None
    data_sigma = None
    data_closed = None
    match plot_type:
        case 'raw':
            print('loading raw data')
            data_raw = np.load('temp/data_raw.npy')
            plt.figure()
            plt.imshow(data_raw[:,:,0], cmap='gray')

        case 'roi':
            print('loading ROI data')
            data_roi = np.load('temp/data_roi.npy')
            plt.figure()
            plt.imshow(data_roi[:,:,0], cmap='gray')

        case 'sigma':
            print('loading sigma data')
            data_sigma = np.load('temp/data_sigma.npy')
            plt.figure()
            plt.imshow(data_sigma[:,:,0], cmap='gray')

        case 'mask':
            print('loading mask data')
            data_mask = np.load('temp/data_mask.npy')
            plt.figure()
            plt.imshow(data_mask[:,:,0], cmap='gray')    

        case 'approx':
            print('loading approx data')
            data_approx = np.load('temp/data_approx.npy')
            plt.figure()
            plt.imshow(data_approx[:,:,0], cmap='gray')
        
        case 'no_bcg':
            print('loading data without background')
            data_no_bcg = np.load('temp/data_no_bcg.npy')
            plt.figure()
            plt.imshow(data_no_bcg[:,:,0], cmap='gray')
            
    plt.show()
    

#------------------------------------------------------------------------------------
# polynomial_mask(data, mask_radius,bits,methode)
# Build an approximated mask on an image. 'nearest' is the fastest methode
# in:
#   data : image to interpolate from
#   mask_radius : mask radius in px
#   bits : Pixel bit depth (8,16) no bits results to a depth of 12 bits in a 16 bits unsigned integer
#   methode : approximation methode ('nearest','linear','cubic')
#
# out: 
#   approximated image
#  
# to do:
#   add mask cendered on speckle or centroid
#------------------------------------------------------------------------------------
def polynomial_mask(data, mask_radius,bits,methode):
    
    # Créer un maillage régulier pour l'approximation polynomiale
    height, width = data.shape

    x = np.linspace(0, width-1, width)
    y = np.linspace(0, height-1, height)

    grid_x, grid_y = np.meshgrid(x, y)
    print('grid ok')

    # Créer les coordonnées pour le masque circulaire
    mask_center_x = width // 2
    mask_center_y = height // 2
    distance_from_center = np.sqrt((grid_x - mask_center_x)**2 + (grid_y - mask_center_y)**2)
    print('coo ok')
    # Créer le masque pour exclure les pixels à l'intérieur du cercle
    mask = distance_from_center <= mask_radius
    
    # Appliquer le masque sur les coordonnées du maillage
    masked_grid_x = grid_x[~mask]
    masked_grid_y = grid_y[~mask]
    masked_image_values = data[~mask]

    print('mask ok')

    match methode:
        case 'nearest':
            met = 'nearest'
        case 'linear':
            met = 'linear'
        case 'cubic':
            met = 'cubic'
        case _:
            met = 'nearest'

    # Effectuer l'approximation polynomiale
    data_approx = griddata((masked_grid_x.flatten(), masked_grid_y.flatten()),
                                  masked_image_values.flatten(),
                                  (grid_x, grid_y),
                                  method=met)
    
    print('approx ok')

    # plt.figure()
    # plt.imshow(data_approx, cmap='gray')
    # plt.colorbar()

    
    #blur
    print('blur')
    total_blured = gaussian_filter(data_approx,sigma=10)
    mask_blured = total_blured
    mask_blured[mask==False] = 0    # ne garde que les valeurs sur le mask



    # plt.figure()
    # plt.imshow(mask_blured, cmap='gray')
    # plt.colorbar()
    data_masked = data.copy()
    data_masked[mask==True] = 0 # ne garde que les valeurs hors du mask
    data_approx = data_masked + mask_blured

    # plt.figure()
    # #plt.imshow(data_approx, cmap='gray')
    # plt.pcolor(data_approx, cmap='gray',vmin=0, vmax=4095)
    # plt.title('approx')
    # plt.colorbar()



    np.save('temp/data_mask',mask_blured, allow_pickle=True, fix_imports=True)
    np.save('temp/data_approx',data_approx, allow_pickle=True, fix_imports=True)
    return data_approx
