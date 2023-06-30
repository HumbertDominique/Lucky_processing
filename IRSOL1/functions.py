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


def readData(n, fname):
    '''------------------------------------------------------------------------------------
    readData(n, fname)

    Reads n .fits files with the path fnameXXX.fits (3 or more X's)

    in:
      n : number of images
      fname : filename
    out:
      data: 3D array with pixel data
    ------------------------------------------------------------------------------------'''

    for i in range(n):
        filename = fname+'{:03d}'.format(i+1)+'.fits'
        if i==0:                    # Gets the pic size on firts pass
            [k, j] = fits.getdata(fname+'001.fits', ext=0).shape
            data = np.zeros((k,j,n))
        data[:,:,i] = fits.getdata(filename, ext=0)
        np.save('temp/data_raw',data, allow_pickle=True, fix_imports=True)
    return data



def ROI(n,data):
    '''------------------------------------------------------------------------------------
    ROI(n, data)

    Applies a 2-D region of interest do a 3-D array 

    in:
        n : number of images
        data : data: 3D array with pixel data
    out:
        data: 3D array with pixel data
    ------------------------------------------------------------------------------------'''
    if data is None:
        print('loading raw data')
        data = np.load('temp/data_raw.npy')
    data = data[580-295:580+294,555-295:555+294,:]     # Manually done

    for i in range(n):
        data[:,:,i] = (data[:,:,i] - np.min(data[:,:,i])) / (np.max(data[:,:,i]) - np.min(data[:,:,i]))
        data[:,:,i] = (data[:,:,i] * 2**12).astype(np.uint16)

    np.save('temp/data_roi',data, allow_pickle=True, fix_imports=True)
    return data



def sigma(n,data):
    '''    ------------------------------------------------------------------------------------
    sigma(n, data)
    in:
    n : number of images
    data : data: 3D array with pixel data
    out:
    data: 3D array with pixel data
    ------------------------------------------------------------------------------------'''
    if data is None:
        print('loading ROI data')
        data = np.load('temp/data_roi.npy')
    [k, j] = data[:,:,1].shape
    data_sigma = np.zeros((k,j,n))
    for i in range(n):
        data_sigma[:,:,i] = stats.sigma_clip(data[:,:,i],maxiters=None)         # A développer
    np.save('temp/data_sigma',data_sigma, allow_pickle=True, fix_imports=True)
    return data_sigma


def display_images(data,n=1):
    '''------------------------------------------------------------------------------------
    display_images(plotType)

    in:
      plotType : omage to plot
    out:
      -
     
    to do:
      add image selection
    ------------------------------------------------------------------------------------'''

    if n > 1:
        for i in range(n):
            plt.figure()
            plt.imshow(data[:,:,i], cmap='gray')
    else:
        plt.imshow(data, cmap='gray')
    plt.show()
    


def buid_mask(data,radius,center=None):
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


def polynomial_mask(data, mask,order=4):
    '''------------------------------------------------------------------------------------
        polynomial_mask(data, mask_radius,order=4)
        according to "Analyse de données"

        in:
          data : image to interpolate from
          mask_radius : mask radius in px
          Order : Order of the appoximation. Default = 4
        
        out: 
          approximated image
         
        to do:
          add mask cendered on speckle or centroid
        ------------------------------------------------------------------------------------'''

    nmodes = int(((order+1)**2-(order+1))/2+(order+1))
    mask_dbl = np.zeros(mask.shape)
    mask_dbl[mask==False] = 1.   # met des 1. là ou le masque n'est pas présent
    width, height = data.shape
    model = np.zeros(mask.shape)

    index_x = np.linspace(-1,1,width)
    index_y = np.linspace(-1,1,height)

    basis = np.zeros((width,height,nmodes))

    Mx = np.zeros((width,height))
    for i in range(height):
        Mx[i,:] = index_x
    My = np.zeros((width,height))
    for i in range(width):
        My[:,i] = index_y
 
    j_mode = -1
    for i in range(order+1):
        for j in range(i+1):
            j_mode +=1
            basis[:,:,j_mode] =Mx**(i-j)*My**j
    print(basis[0:5,0:5,11])
    A = np.zeros((nmodes,nmodes)).astype(np.double)

    for i in range(0,nmodes):
        for j in range(0,i+1): 
            A[i,j] = np.sum(basis[:,:,i]*basis[:,:,j]*mask_dbl)
            #A[j,i] = A[i,j]
            #A=A.T
            #print(i,j)
    print(A[1,0])

    b = np.zeros(nmodes)
    for i in range(nmodes):
         b[i]=np.sum(data*basis[:,:,i]*mask_dbl)

    a = np.dot(np.linalg.inv(A),b)
    model = np.zeros((width,height))
    for i in range(0,nmodes):
        model += a[i]*basis[:,:,i] 
    
    return model