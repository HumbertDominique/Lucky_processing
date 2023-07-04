# Lecture d'image en lucky imaging.
# 
# Inspired and somtime simply translated form process2.pro (L.J)
#
# 2023.06.26    -   Dominique Humbert   -   Version initiale

from astropy.io import fits
import numpy as np
from customClasses import eErrors
import matplotlib.pyplot as plt

'''------------------------------------------------------------------------------------
    # eError_handling(code)
    Error handling routine and help
    ## in:
      code: errore code
    ## out:
        -    ------------------------------------------------------------------------------------'''
def eError_handling(code):
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
    5: Background removal\n\
    6: Noise cut (if implemented)\n\
    7: -\n\
    8: Gausian oise filtering (if implemented)\n\
    9: Dust shadow removal\n\
    10: tracking\n\
    11: Selection\n\
    12: plots")
            print("arg2 = fraction of frames to use (ie. 0.2)")
            print("--------------------------")


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
        if i==0:                    # Gets the pic size on firts pass
            [k, j] = fits.getdata(fname+'001.fits', ext=0).shape
            data = np.zeros((k,j,n))
        data[:,:,i] = fits.getdata(filename, ext=0)

        if 100*(i+1)//n == status*10:
            print('Reading data: ',status*10,'% done')
            status += 1
    print('Saving intermediate data')
    np.save('temp/data_raw',data, allow_pickle=True, fix_imports=True)
    print('Done')
    return data



def ROI(n,data_in):
    '''------------------------------------------------------------------------------------
    # ROI(n, data)

    Applies a 2-D region of interest do a 3-D array 

    ## in:
        n: number of images
        data: data: 3D array with pixel data
        bit_depth: 2^bits resolution (1<=bits<=16)
    ## out:
        data: 3D array with pixel data
    ------------------------------------------------------------------------------------'''
    if data_in is None:
        print('Loading raw data')
        data_in = np.load('temp/data_raw.npy')

    data = data_in[530-256:530+256,600-256:600+256,:]     # Manually done
    print('Saving intermediate data')
    np.save('temp/data_roi',data, allow_pickle=True, fix_imports=True)
    print('Done')

    return data


def normalise(data, bits_depth = 12,):
    '''------------------------------------------------------------------------------------
    # normalize(data,bits_depth)

    Normalise images 

    ## in:
        data : data: 3D array with pixel data
        bit_depth : 2^bits resolution (1<=bits<=16)
    ## out:
        data: 3D array with pixel data
    ------------------------------------------------------------------------------------'''
    

    if bits_depth<0 or bits_depth>16:
        bits = 12
        #print('Normalising images in ',bits,'bits')
        data[:,:] = (data[:,:] - np.min(data[:,:])) / (np.max(data[:,:]) - np.min(data[:,:]))
        data[:,:] = (data[:,:] * 2**bits).astype(np.uint16)

    if bits_depth == 0:
        data[:,:] = ((data[:,:] - np.min(data[:,:])) / (np.max(data[:,:]) - np.min(data[:,:]))).astype(np.double)  

    return data


def buid_mask(data,radius,center=[None,None]):
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

    height,width = data.shape

    if not (center[0] and center[1]):
        center = np.zeros((2))
        center[0] = height // 2  # x
        center[1] = width // 2 # x

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
          mask : mask boolean mask. True means to be removed
          Order : Order of the appoximation. Default = 4
        
        out: 
          approximated image
         
        to do:
        ------------------------------------------------------------------------------------'''

    nmodes = int(((order+1)**2-(order+1))/2+(order+1))
    mask_dbl = np.zeros(mask.shape)
    mask_dbl[mask==False] = 1.   # met des 1. là ou le masque n'est pas présent
    height, width = data.shape
    model = np.zeros(mask.shape)

    index_x = np.linspace(-1,1,width)
    index_y = np.linspace(-1,1,height)

    basis = np.zeros((height,width,nmodes))

    Mx = np.zeros((height,width))
    My = np.zeros((height,width))
    for i in range(width):
        Mx[:,i] = index_x
    for i in range(height):    
        My[i,:] = index_y
 
    j_mode = -1
    for i in range(order+1):
        for j in range(i+1):
            j_mode +=1
            basis[:,:,j_mode] =Mx**(i-j)*My**j
    A = np.zeros((nmodes,nmodes)).astype(np.double)

    for i in range(0,nmodes):
        for j in range(0,i+1): 
            A[i,j] = np.sum(basis[:,:,i]*basis[:,:,j]*mask_dbl)
            #A[j,i] = A[i,j]
            #A=A.T
            #print(i,j)

    b = np.zeros(nmodes)
    for i in range(nmodes):
         b[i]=np.sum(data*basis[:,:,i]*mask_dbl)

    a = np.dot(np.linalg.inv(A),b)
    model = np.zeros((height,width))
    for i in range(0,nmodes):
        model += a[i]*basis[:,:,i] 
    
    return model


def noise_cut(data, threshold=None):
    '''------------------------------------------------------------------------------------
        noise_cut(data, threshold=None)
        Noise cut filter

        in:
          data : image to interpolate from. square and dimmentions are ^2 (padding not implemented)
          threshold : cutting thresholde
        
        out: 
          approximated image
         
        to do:
            faire fonctionner
        ------------------------------------------------------------------------------------'''

    if threshold is None:
        threshold = 0.5*1                  # à définir avec Nyquist

    # Perform Fourier transform
    data_tf = np.fft.fftshift(np.fft.fft2(data))
    
    # Apply filter in Fourier domain
    mask = np.abs(data_tf) > threshold
    data_tf_filtered = data_tf*mask
    
    # Inverse Fourier transform
    data = np.fft.ifft2(np.fft.ifftshift(data_tf_filtered)).real

    
    return data.astype(np.double)

def noise_filter(data,sigma=None):
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

    if sigma is None:
        print('Input image Std for gaussian noise filtering')
    else:
        height, width = data.shape

        index_x = np.linspace(-1,1,width)
        index_y = np.linspace(-1,1,height)

        Mx = np.zeros((height,width))
        My = np.zeros((height,width))
        for i in range(height):
            Mx[:,i] = index_x
        for i in range(width):    
            My[i,:] = index_y

        # Perform Fourier transform
        data_tf = np.fft.fftshift(np.fft.fft2(data))
        
        # Apply filter in Fourier domain

        
        G = (1/(2*np.pi*sigma)**0.5)*np.exp((Mx**2+My**2)/(2*sigma**2))

        data_tf_filtered = data_tf*G
        
        # Inverse Fourier transform
        data = np.fft.ifft2(np.fft.ifftshift(data_tf_filtered)).real


    return data.astype(np.double)

def image_dusting(data,model, NdustSpots=1):
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

