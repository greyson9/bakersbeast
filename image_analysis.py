# image_analysis.py
# This script analyzes microscopy data from our Tamoxifen perturbation expt's
# Channel 1 is DAPI (blue nucleic acid stain)
# Channel 2 is FITC (100 ms)
# Channel 3 is FITC (1 s)
# Channel 4 is Cy3 (red/orange cell wall stain)
# Channel 5 is bright field
from skimage.feature import canny
from skimage.filters import threshold_otsu, threshold_isodata
from skimage.data import text
import skimage.morphology as morph
from skimage import io, filters
import skimage as sk
import numpy as np
import scipy.ndimage as ndi
import matplotlib.pyplot as plt
from scipy import stats
import nd2reader as nd2r

stack = nd2r.Nd2('/Users/student/Projects/PlateImages/Plate000_WellC01_Seq0000.nd2')
test = stack[4].astype(np.float64) / np.amax(stack[4])
threshold_global_otsu = threshold_otsu(test)
print(threshold_global_otsu)
global_otsu = test >= threshold_global_otsu
# threshold_isodata = threshold_isodata(test)
# isodata = test >= threshold_isodata
skel_test = morph.skeletonize(global_otsu)
close_test = ndi.binary_closing(skel_test, iterations=7	)
fill_cells = ndi.binary_fill_holes(close_test)
print(global_otsu)
plt.figure()
plt.hist(stack[0].astype(np.float64) / np.amax(stack[0]), bins=500)
# plt.figure()
# io.imshow(global_otsu)
# plt.figure()
# io.imshow(fill_cells)
# plt.figure()
# io.imshow(skel_test)
# plt.figure()
# io.imshow(close_test)
# plt.figure()
# io.imshow(fill_cells)
plt.savefig('dapihist.png')

# fig, axes = plt.subplots(3, 2, figsize=(8, 5), sharex=True, sharey=True,
# 							subplot_kw={'adjustable': 'box-forced'})
# ax = axes.ravel()
# plt.tight_layout()

# fig.colorbar(ax[0].imshow(test, cmap=plt.cm.gray),
#              ax=ax[0], orientation='horizontal')
# ax[0].set_title('Original')
# ax[0].axis('off')

# fig.colorbar(ax[1].imshow(local_otsu, cmap=plt.cm.gray),
#              ax=ax[1], orientation='horizontal')
# ax[1].set_title('Local Otsu (radius=10)')
# ax[1].axis('off')

# ax[2].imshow(test >= local_otsu, cmap=plt.cm.gray)
# ax[2].set_title('Original >= Local Otsu')
# ax[2].axis('off')

# ax[3].imshow(global_otsu, cmap=plt.cm.gray)
# ax[3].set_title('Global Otsu (threshold = %f)' % threshold_global_otsu)
# ax[3].axis('off')

# ax[4].imshow(test >= yen_thresh, cmap=plt.cm.gray)
# ax[4].set_title('Original >= Yen (threshold = %f)' % yen_thresh)
# ax[4].axis('off')

# ax[5].imshow(edges, cmap=plt.cm.gray)
# ax[5].set_title('Canny edging, low thresh is 0.95 * global otsu')
# ax[5].axis('off')

# plt.show()
# plt.figure()
# io.imshow(out)
# plt.figure()
# plt.hist(out)
# fill_cells = ndi.binary_fill_holes(edges)
# print(dir(filters.thresholding))
# plt.figure()
# io.imshow(edges)
# print(mod[500][500:600])
# plt.figure()
# io.imshow(test)
# plt.figure()
# io.imshow(mod)
# plt.figure()
# io.imshow(edges)
# plt.figure()
# io.imshow(watershed)
# plt.figure()
# io.imshow(test)
# plt.figure()
# io.imshow(edges)
# io.show()
