from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import functions as fnc


pic = fits.getdata('C:/Users/ADM/OneDrive - HESSO/Dominique/07_mesures/07_aberation_disk/AD_010.fits', ext=0)
dark = fits.getdata('C:/Users/ADM/OneDrive - HESSO/Dominique/07_mesures/07_aberation_disk/dark1.fits', ext=0)
flat = fits.getdata('C:/Users/ADM/OneDrive - HESSO/Dominique/07_mesures/07_aberation_disk/flat_norm.fits', ext=0)


temp = pic
temp = (temp-np.min(temp))/(np.max(temp)-np.min(temp))
plt.figure()
plt.imshow(temp**(1/2),vmin=0,vmax=1)
plt.title('Raw image normalized')

temp = pic-dark
temp = (temp-np.min(temp))/(np.max(temp)-np.min(temp))
plt.figure()
plt.imshow(temp**(1/2),vmin=0,vmax=1)
plt.title('(Raw - Dark) normalized')

temp = (pic-dark)/flat
temp = (temp-np.min(temp))/(np.max(temp)-np.min(temp))
plt.figure()
plt.imshow(temp**(1/2),vmin=0,vmax=1)
plt.title('(Raw - Dark)/flat_norm normalized')



#----------------- BCG removal (2 masked zones)--------------

# speckles = np.zeros((2,2))
# speckles[:,0] = [209, 481]
# speckles[:,1] = [910, 462]

# height, width = pic.shape
# mask = np.zeros((height, width,3))  # Remove if using mask and model at the same time as data takes too much
# model = np.zeros((height, width))


# # -----------    Build grids
# x = np.linspace(0, width-1, width).astype(np.uint16)
# y = np.linspace(0, height-1, height).astype(np.uint16)
# grid_x_UINT, grid_y_UINT = np.meshgrid(x, y)

# x = np.linspace(-1,1,width)
# y = np.linspace(-1,1,height)
# grid_x, grid_y = np.meshgrid(x, y)
# # -----------    Build grids

# mask_radius = 150
# mask[:,:,0] = fnc.buid_mask(pic,mask_radius,[grid_x_UINT, grid_y_UINT],speckles[:,0])
# mask_radius = 80
# mask[:,:,1] = fnc.buid_mask(pic,mask_radius,[grid_x_UINT, grid_y_UINT],speckles[:,1])
# mask[:,:,2] = mask[:,:,0] + mask[:,:,1]

# model, poly_error_flag = fnc.polynomial_mask(pic,mask[:,:,2],[grid_x, grid_y],4)
# result = pic-model
# result = (result-np.min(result))/(np.max(result)-np.min(result))

plt.show()
