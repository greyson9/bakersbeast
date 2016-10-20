
import nd2reader
import skimage
import skimage.filters as filters
from skimage.morphology import watershed, disk
from skimage.filters import sobel
import skimage.feature as feature
from skimage import exposure
from skimage import transform as tf
import skimage.segmentation as seg
import skimage.morphology as morph
from skimage.segmentation import slic, join_segmentations
from skimage import data,io,filters,img_as_float
from skimage.color import rgb2gray
from scipy import fftpack, stats
from skimage.filters.rank import median
from scipy import ndimage as ndi
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from os import listdir

from skimage.feature import peak_local_max

from skimage.morphology import disk
from skimage.filters import threshold_otsu, rank
from skimage.util import img_as_ubyte


# img1= nd2reader.Nd2('Plate000_WellC12_Seq0011.nd2')
#img2= nd2reader.Nd2()
# print img1

def segment_image(photo_matrix, background_threshold, foreground_threshold):
	edges = sobel(photo_matrix)
	markers = np.zeros_like(photo_matrix)
	foreground, background = 1, 2
	markers[photo_matrix < background_threshold] = background
	markers[photo_matrix > foreground_threshold] = foreground
	ws = watershed(edges, markers)
	segmentation_matrix = ndi.label(ws == foreground)[0]
	return segmentation_matrix

def normalize(img):
	high = np.amax(img)
	low = np.amin(img)
	return (img - low) / (high - low)

def find_cells(img):
	strong_blur = filters.gaussian(img, 20)
	no_back = img - strong_blur
	no_back = normalize(no_back)
	equalized_no_back = exposure.equalize_hist(no_back)
	equalized_no_back = normalize(equalized_no_back)
	edges_nb = feature.canny(equalized_no_back, sigma=5)
	close_nb = ndi.binary_closing(edges_nb, structure=np.ones((3, 3)), iterations=1)
	fill_close_nb = ndi.binary_fill_holes(close_nb)
	open_fcnb = ndi.binary_opening(fill_close_nb, structure=np.ones((10, 10)))
	open_bigger = ndi.morphology.binary_dilation(open_fcnb, iterations=5)
	border = morph.binary_dilation(open_fcnb) - open_fcnb
	dist = ndi.distance_transform_edt(open_fcnb)
	local_peaks = feature.peak_local_max(dist, min_distance=12, threshold_abs=4,
											labels=open_fcnb, indices=False)
	markers = ndi.label(local_peaks)[0]
	labels = morph.watershed(-dist, markers, mask=open_bigger)
	find_boundaries = seg.find_boundaries(labels)
	return labels

# test_image1=img1[4].astype(np.float64)/np.amax(img1[4])
# test_image2=img1[0].astype(np.float64)/np.amax(img1[0])
# test_image3=img1[2].astype(np.float64)/np.amax(img1[2])
# matrix=segment_image(test_image1,0.4,0.45)
# matrix2=segment_image(test_image2,0.5,0.5)
# matrix=find_cells(test_image1)
# matrix2=find_cells(test_image2)


def find_intensity_cell(matrix,image):
	avg_intensity=[]
	index=[]
	for i in range(1,np.amax(matrix)):
		index=zip(*np.where(matrix==i))
		intensity=0
		k=0
		for j in range(0,len(index)):
			intensity=intensity+image.item(index[j][0],index[j][1])
			k=k+1
		avg_intensity.append(intensity/k)
	return avg_intensity

def total_intensity(matrix,image):
	tot_intensity=[]
	index=[]
	for i in range(1,np.amax(matrix)):
		index=zip(*np.where(matrix==i))
		intensity=0
		for j in range(0,len(index)):
			intensity=intensity+image.item(index[j][0],index[j][1])
		tot_intensity.append(intensity)
	return tot_intensity

def find_intensity_in_nuc(matrix,gfp_image):
	avg_intensity=[]
	index=[]
	for i in range(1,np.amax(matrix)):
		index=zip(*np.where(matrix==i))
		intensity=0
		k=0
		for j in range(0,len(index)):
			intensity=intensity+gfp_image.item(index[j][0],index[j][1])
			k=k+1
		avg_intensity.append(intensity/k)
	return avg_intensity

# intensity1= find_intensity_cell(matrix,test_image3)
# intensity_nuc=find_intensity_in_nuc(matrix2,test_image3)
# tot_int=total_intensity(matrix,test_image3)
# tot_int_nuc=total_intensity(matrix2,test_image3)
#
# sns.distplot(intensity1, len(np.unique(intensity1)))
# sns.distplot(intensity_nuc, len(np.unique(intensity_nuc)))

# sns.distplot(tot_int, len(np.unique(tot_int)))
# sns.distplot(tot_int_nuc, len(np.unique(tot_int_nuc)))
#plt.axis([0, 0.7, 0, 100])
# cbar = plt.colorbar()
# plt.xlabel("Value")
# plt.ylabel("Frequency")
#plt.imshow(matrix2,cmap='jet')
#plt.imshow(test_image3,alpha=0.5)


# print 'cell mean', np.mean(intensity1)
# print 'cell var', np.var(intensity1)
# print 'nuc mean', np.mean(intensity_nuc)
# print 'nuc var', np.var(intensity_nuc)
# print 'nuc/cell', np.mean(intensity_nuc)/np.mean(intensity1)
# plt.show()

#plt.savefig(name.png)

for filename in listdir('PlateImages1/'):
	if filename == '.DS_Store':
		continue
	img1= nd2reader.Nd2('PlateImages1/' + filename)
	test_image1=img1[4].astype(np.float64)/np.amax(img1[4])
	test_image2=img1[0].astype(np.float64)/np.amax(img1[0])
	test_image3=img1[2].astype(np.float64)/np.amax(img1[2])
	matrix=find_cells(test_image1)
	matrix2=segment_image(test_image2,0.5,0.5)
	intensity= find_intensity_cell(matrix,test_image3)
	intensity_nuc=find_intensity_in_nuc(matrix2,test_image3)
	sns.distplot(intensity, len(np.unique(intensity)))
	sns.distplot(intensity_nuc, len(np.unique(intensity_nuc)))
	plt.xlabel("Value")
	plt.ylabel("Frequency")
	plt.savefig('PlateImages2/' + filename+ '.png')
	plt.close()
