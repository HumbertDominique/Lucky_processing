from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import functions as fnc


foc = fits.getdata('C:/Users/dominiqu.humbert/OptoLab Dropbox/OptoLabMAIN/10_IRSOL/telescopeCharacterisation/Phase_divsity_tool/07_data/09_Characterization_data/mesure_setup/Measure_setup_001.fits', ext=0)
#defdofc = fits.getdata('C:/Users/ADM/OneDrive - HESSO/Dominique/07_mesures/08_aberations_with_foen/temp/mean_image_Stracking_def.fits', ext=0)



plt.figure()
plt.imshow(foc**(1/1))
plt.title('Focused images Stracked',cmap='gray')


# plt.figure()
# plt.imshow(defdofc**(1/1))
# plt.title('Defocused images Stracked')



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
