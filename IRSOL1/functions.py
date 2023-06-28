# Lecture d'image en lucky imaging.
# 
#
# 2023.06.26    -   Dominique Humbert   -   Version initiale
from astropy.io import fits
import numpy as np
from customClasses import eErrors
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
            print("arg1 = process\n1: Whole process\n2: ROI\n3: Sigma filter\n10: plot")
            print("arg2 = plot type")
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
        np.save('data_raw',data, allow_pickle=True, fix_imports=True)
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
        data = np.load('data_raw.npy')
    data = data[500-295:500+294,625-295:625+294,:]     # Manually done
    np.save('data_roi',data, allow_pickle=True, fix_imports=True)
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
        data = np.load('data_roi.npy')
    [k, j] = data[:,:,1].shape
    data_sigma = np.zeros((k,j,n))
    for i in range(n):
        data_sigma[:,:,i] = stats.sigma_clip(data[:,:,i],maxiters=None)         # A développer
    np.save('data_sigma',data_sigma, allow_pickle=True, fix_imports=True)
    return data_sigma

#------------------------------------------------------------------------------------
# display_images(plotType)
# in:
#   plotType : omage to plot
#out:
#   -
#------------------------------------------------------------------------------------
def display_images(plotType):
    data_raw = None
    data_roi = None
    data_sigma = None
    data_closed = None

    if plotType == "raw":
        print('loading raw data')
        data_raw = np.load('data_raw.npy')
        plt.figure()
        plt.imshow(data_raw[:,:,0], cmap='gray')
        plt.show()
    if plotType == "roi":
        print('loading ROI data')
        data_roi = np.load('data_roi.npy')
        plt.figure()
        plt.imshow(data_roi[:,:,0], cmap='gray')
        plt.show()
    if plotType == "sigma":
        print('loading sigma data')
        data_sigma = np.load('data_sigma.npy')
        plt.figure()
        plt.imshow(data_sigma[:,:,0], cmap='gray')
        plt.show()
    

#------------------------------------------------------------------------------------
# polynomial_mask(data, mask_radius,bits,methode)
# Build an approximated mask on an image. 'nearest' is the fastest methode
# in:
#   data : image to interpolate from
#   mask_radius : mask radius in px
#   bits : Pixel bit depth (8,12,16)
#   methode : approximation methode ('nearest','linear','cubic')
#
#out: 
#   approximated image
#  
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
    approximated_image = griddata((masked_grid_x.flatten(), masked_grid_y.flatten()),
                                  masked_image_values.flatten(),
                                  (grid_x, grid_y),
                                  method=met)
    print('approx ok')

    # Normaliser les valeurs de l'image approximée
    approximated_image = (approximated_image - np.min(approximated_image)) / (np.max(approximated_image) - np.min(approximated_image)) # Décends toutes les valeurs du noise floor. puis normalise. On a donc des valeurs de 0 à 1
    print('norm ok')

     # Convertir les valeurs normalisées en entiers x bits
    match bits:
        case 8:
            approximated_image = (approximated_image * 2**bits).astype(np.uint8)
        case 12:
            approximated_image = (approximated_image * 2**bits).astype(np.uint12)
        case 16:
            approximated_image = (approximated_image * 2**bits).astype(np.uint16)
    print('cast ok')

    if met == 'nearest':
        #blur
        
        total_blured = gaussian_filter(approximated_image,sigma=10)
        mask_blured = np.dstack((total_blured, ~mask))



  
    
   

    return approximated_image
