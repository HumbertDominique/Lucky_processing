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
def display_images(data,n=1):

    if n > 1:
        for i in range(n):
            plt.figure()
            plt.imshow(data[:,:,i], cmap='gray')
    else:
        plt.imshow(data, cmap='gray')
    plt.show()
    

#------------------------------------------------------------------------------------
# buid_maks(data,center=None)
# in:
#   data : Array the size of the image for whoch the mask is built
#   radius : mask radius
#   center (optional) : x,y coordinate of the mask's center 
#out:
#   mask : boolean array with True in the mask circle
#------------------------------------------------------------------------------------
def buid_mask(data,radius,center=None):

    width, height = data.shape

    if center == None:
        center = np.zeros((2))
        center[0] = width // 2  # x
        center[1] = height // 2 # x

    x = np.linspace(0, width-1, width).astype(np.uint16)
    y = np.linspace(0, height-1, height).astype(np.uint16)
    grid_x, grid_y = np.meshgrid(x, y)
                    
    distance_from_center = np.sqrt((grid_x - center[0])**2 + (grid_y - center[1])**2)

    # Créer le masque pour exclure les pixels à l'intérieur du cercle
    mask = distance_from_center <= radius

    return mask

#------------------------------------------------------------------------------------
# polynomial_mask(data, mask_radius,order=4)
# according to "Analyse de données"
# in:
#   data : image to interpolate from
#   mask_radius : mask radius in px
#   Order : Order of the appoximation. Default = 4
#
# out: 
#   approximated image
#  
# to do:
#   add mask cendered on speckle or centroid
#------------------------------------------------------------------------------------
def polynomial_mask(data, mask,order=4):

    mask_dbl = np.zeros(mask.shape)
    mask_dbl[mask==False] = 1.   # met des 1. là ou le masque n'est pas présent
    width, height = data.shape

    a = np.zeros((width, height, order))
    b = a.copy()
    A = np.zeros((width, height, order)).astype(np.double)

    #Coordinate matrix
    xs = np.linspace(-1, 1, width).astype(np.double)
    ys = np.linspace(-1, 1, height).astype(np.double)

    index_grid_x = np.zeros((width, height)).astype(np.double)
    index_grid_y = np.zeros((width, height)).astype(np.double)

    for i in range(height):
        index_grid_x[:,i] = xs
    for i in range(width):
        index_grid_y[i,:] = ys

    o = -1   
    for i in range(order):
        o=o+1
        for j in range(0,i+1):
            A[:,:,o] = index_grid_x**(i-(j-1))*index_grid_y**(i)

    Gij = np.zeros((order,order))

    for o in range(order):
        for i in range(0,o+1):
            Gij[o,i] = np.sum(A[:,:,o]*A[:,:,i-1]*mask_dbl)
            Gij[o,i] = Gij[i,o]

    b = np.zeros((order))
    for o in range(order):      # Calcule b = [mi, xi*mi, yi*mi,xi^2*mi, xi*yi*mi, yi^2*mi, ...] pour chaque pixel
        b[o] = np.sum((data*A[:,:,o]*mask_dbl))  # outside of mask
    a=np.zeros(b.shape)
    for o in range(order):
        a = np.dot(b,np.linalg.inv(Gij)).copy()

    model = np.zeros((width, height))
    for i in range(order):
        model=model+a[i]*A[:,:,i]

    return model



    




#------------------------------------------------------------------------------------
# polynomial_mask_py(data, mask_radius,bits,methode)
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
def polynomial_mask_py(data, mask_radius,methode):
    
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
