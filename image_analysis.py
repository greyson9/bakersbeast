# image_analysis.py
# This script analyzes microscopy data from our Tamoxifen perturbation expt's
# Channel 1 is DAPI (blue nucleic acid stain)
# Channel 2 is FITC (100 ms)
# Channel 3 is FITC (1 s)
# Channel 4 is Cy3 (red/orange cell wall stain)
# Channel 5 is bright field
import skimage.morphology as morph
import skimage as sk
import numpy as np
import scipy.ndimage as ndi
import skimage.measure as measure
import matplotlib.pyplot as plt
from scipy import stats
import nd2reader as nd2r
import seaborn as sns
import skimage.feature as feature
import skimage.filters as filters
import skimage.io as io
import skimage.exposure as exposure


def double_otsu(img, bins):
	hist, bin_edges = np.histogram(img.flatten(), bins=bins)
	t1 = 0
	t2 = 0
	w0k = 0
	w1k = 0
	m0k = 0
	m1k = 0
	length = len(img.flatten())
	mt = [u * hist[u] / length for u in range(0, bins)]
	mt = sum(mt)
	maxBetweenVar = 0
	for i in range(0, bins):
		w0k += hist[i] / length
		m0k += i * hist[i] / length
		w1k = 0
		m1k = 0
		m0 = m0k / w0k
		for j in range(i, bins):
			w1k += hist[j] / length
			m1k += j * hist[j] / length
			m1 = m1k / w1k
			w2k = 1 - (w0k + w1k)
			m2k = mt - (m0k + m1k)

			if w2k <= 0:
				break

			m2 = m2k / w2k
			currVarB = w0k * (m0 - mt) ** 2 + w1k * (m1 - mt) ** 2 + w2k * (m2 - mt) ** 2

			if maxBetweenVar < currVarB:
				maxBetweenVar = currVarB
				t1 = i
				t2 = j
	return 1.0 * t1 / bins, 1.0 * t2 / bins

def normalize(img):
	high = np.amax(img)
	low = np.amin(img)
	return (img - low) / (high - low)


stack = nd2r.Nd2('/Users/student/Projects/PlateImages/Plate000_WellC01_Seq0000.nd2')
test = stack[4].astype(np.float64) / np.amax(stack[4])
strong_blur = filters.gaussian(test, 20)
no_back = test - strong_blur
no_back = normalize(no_back)
equalized_no_back = exposure.equalize_hist(no_back)
equalized_no_back = normalize(equalized_no_back)
edges_nb = feature.canny(equalized_no_back, sigma=4)
close_nb = ndi.binary_closing(edges_nb, structure=np.ones((3, 3)), iterations=1)
fill_close_nb = ndi.binary_fill_holes(close_nb)
open_fcnb = ndi.binary_opening(fill_close_nb, structure=np.ones((10, 10)))
border = morph.binary_dilation(open_fcnb) - open_fcnb
dist = ndi.distance_transform_cdt(open_fcnb, metric='chessboard')
dist_bins = np.amax(dist)
# dist = normalize(dist)
local_peaks = feature.peak_local_max(dist, min_distance=3, threshold_abs=2,
										labels=open_fcnb, indices=False)
hist, bin_edges = np.histogram(dist.flatten(), bins=dist_bins)
dist_otsu = filters.threshold_otsu(dist)
markers = ndi.label(local_peaks)[0]
labels = morph.watershed(-dist, markers, mask=no_back)
# print(local_peaks)
# laplace = filters.scharr(no_back)
# laplace = (laplace - np.amin(laplace)) / (np.amax(laplace) - np.amin(laplace))
# adapt = filters.threshold_adaptive(image=laplace, block_size=251)
# skel = morph.skeletonize(adapt)
# close = ndi.binary_closing(skel, iterations=3)
# fill = ndi.binary_fill_holes(skel)
# opened = ndi.morphology.binary_opening(fill, structure=np.ones((10, 10)))
# dist = ndi.distance_transform_cdt(opened, metric='chessboard')
# hist, bin_edges = np.histogram(dist.flatten(), bins=np.amax(dist))

# markers = filters.rank.gradient(no_back, morph.disk(5)) > 10
# markers = ndi.label(markers)[0]
# grad = filters.rank.gradient(no_back, morph.disk(3))
# labels = morph.watershed(grad, markers)
# print(markers)
# second = (second - np.amin(second)) / (np.amax(second) - np.amin(second))
# u = filters.threshold_otsu(laplace)
# print(u)
# print((int)(u * 1000.0), (int)(v * 1000.0))
# filt_real, filt_imag = filters.gabor(test, frequency=0.8)
plt.figure()
io.imshow(labels)
# plt.figure()
# sns.barplot(x=bin_edges[1:-1], y=hist[1:])
# plt.axvline(x=(int)(1000.0 * otsu), color='k', linestyle='--')
# plt.axvline(x=(int)(1000.0 * u), color='k', linestyle='--')
# plt.axvline(x=(int)(v * 1000.0), color='k', linestyle='--')
# plt.figure()
# io.imshow(border)
io.show()
'''
# no_back = no_back - np.amin(no_back)
print(np.amin(no_back))
u, v = double_otsu(test)
otsu = filters.threshold_otsu(test)
print(u, v)
equalized = exposure.equalize_hist(test)
equalized_no_back = exposure.equalize_hist(no_back)

hist, bin_edges = np.histogram(no_back.flatten(), bins=500)
plt.figure()
io.imshow(no_back)
plt.figure()
io.imshow(equalized_no_back)
plt.figure()
# io.imshow(strong_blur)
sns.barplot(x=bin_edges[:-1], y=hist)
# xpos = [(int)(u * 500.0), (int)(v * 500.0)]
# for xc in xpos:
# 	plt.axvline(x=xc, color='k', linestyle='--')
# plt.axvline(x=(int)(otsu*500.0), color='k', linestyle='--')
io.show()
threshold_global_otsu = features.threshold_otsu(test)
threshold_mean = np.mean(test.flatten())
print(threshold_global_otsu)
print(threshold_mean)
global_otsu = test >= threshold_global_otsu
mean = test >= threshold_mean
u, v = double_otsu(test)
print(u, v)
ot1 = (test >= u) & (test <= v)
ot2 = test >= v
print(ot2)
elevation = sobel(test)
el2 = laplace(test, ksize=3)
markers = np.zeros_like(test)
markers[test >= v] = 2
markers[test <= u] = 1
seg = morph.watershed(el2, markers)
fill = ndi.binary_fill_holes(seg)
blur1 = filters.gaussian(test, sigma=0.25)
blur2 = filters.gaussian(test, sigma=0.5)
blur3 = filters.gaussian(test, sigma=1.0)
blur4 = filters.gaussian(test, sigma=2.0)
b10 = threshold_otsu(blur1)
b11, b12 = double_otsu(blur1)
b20 = threshold_otsu(blur2)
b21, b22 = double_otsu(blur2)
b30 = threshold_otsu(blur3)
b31, b32 = double_otsu(blur3)
b40 = threshold_otsu(blur4)
b41, b42 = double_otsu(blur4)
plt.figure()
io.imshow(blur1 >= b11)
plt.figure()
io.imshow(blur2 >= b21)
plt.figure()
io.imshow(blur3 >= b31)
plt.figure()
io.imshow(blur4 >= b41)
# plt.figure()
# io.imshow(walls)
io.show()
'''
# doh = blob_doh(test)
# dog = blog_dog(test)
# log = blog_log(test)
# edges = canny(test, low_threshold=u, high_threshold=v)
# filled = ndi.binary_fill_holes(close_otsu)
# edges_q = canny(test, low_threshold=u, high_threshold=v, use_quantiles=True)
# plt.figure()
# io.imshow(test)
# plt.figure()
# io.imshow(test < u)
# plt.figure()
# io.imshow(ot1)
# plt.figure()
# io.imshow(ot2)
# plt.figure()
# io.imshow(elevation)
# plt.figure()
# io.imshow(el2)
# plt.figure()
# io.imshow(seg)
# plt.figure()
# io.imshow(fill)
# io.show()
# skel_otsu = morph.skeletonize(global_otsu)
# skel_mean = morph.skeletonize(mean)
# plt.figure()
# io.imshow(global_otsu)
# plt.figure()
# io.imshow(mean)
# io.show()
# threshold_isodata = threshold_isodata(test)
# isodata = test >= threshold_isodata
# skel_test = morph.skeletonize(global_otsu)
# close_mean = ndi.binary_closing(mean, iterations=2)
# plt.figure()
# io.imshow(close_otsu)
# plt.figure()
# io.imshow(close_mean)
# io.show()
# fill_cells = ndi.binary_fill_holes(close_test)

# hist, bin_edges = np.histogram(stack[0].astype(np.float64).flatten() /
# 					np.amax(stack[0]), bins=1000)
# print(bin_edges)
# plt.figure()
# sns.barplot(x=bin_edges[:-1], y=hist)
# plt.xlim(min(bin_edges), max(bin_edges))
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
# plt.show()
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
