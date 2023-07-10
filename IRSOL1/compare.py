from astropy.io import fits
import matplotlib.pyplot as plt

ref = fits.getdata('immoy7.fits', ext=0)
mine = fits.getdata('Attempt1/Attempt1_1146.fits', ext=0)

# plt.figure()
# plt.imshow(ref)

plt.figure()
plt.imshow(mine)

# test = ref[:,:] - mine

# plt.figure()
# plt.imshow(test)

plt.show()
