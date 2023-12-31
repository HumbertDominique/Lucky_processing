from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy import stats

flat = fits.getdata('C:/Users/ADM/OneDrive - HESSO/Dominique/07_mesures/08_aberations_with_foen/temp/flat.fits', ext=0)
dark = fits.getdata('C:/Users/ADM/OneDrive - HESSO/Dominique/07_mesures/08_aberations_with_foen/temp/dark.fits', ext=0)

plt.figure()
plt.imshow(flat)
plt.title('Flat field pre-normalisation')
plt.figure()
plt.plot(flat[450,:])
plt.title('Flat field pre-normalisation')

plt.figure()
plt.imshow(dark)
plt.title('Dark field')

temp = flat-dark
temp = (temp-np.min(temp))/(np.max(temp)-np.min(temp))
flat = temp

plt.figure()
plt.imshow((temp))
plt.title('(Flat - Dark) normalized')

plt.figure()
plt.plot(temp[450,:])
plt.title('(Flat - Dark) normalized')


plt.show()

hdu = fits.PrimaryHDU(flat)
hdu.writeto('C:/Users/ADM/OneDrive - HESSO/Dominique/07_mesures/08_aberations_with_foen/temp/flat_norm.fits',overwrite=True)