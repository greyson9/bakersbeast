import skimage.morphology as morph
import skimage as sk
import numpy as np
import scipy.ndimage as ndi
import skimage.measure as measure
from skimage import color
import matplotlib.pyplot as plt
from scipy import stats
import nd2reader as nd2r
import seaborn as sns
import skimage.feature as feature
import skimage.filters as filters
import skimage.io as io
import skimage.exposure as exposure
import skimage.segmentation as seg
from os import listdir
import math
import re
import cPickle as pik 

# normalize imported images
def normalize(img):
	
	img = img.astype(np.float64)
	high = np.amax(img)
	low = np.amin(img)
	return (img - low) / (high - low)

# segment cells
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


# eliminate things that are too large to be cells
def remove_large_things(labels, img):
	regions = measure.regionprops(labels, img)
	areas = []
	perimeters = []
	too_big = []
	for region in regions:
		areas.append(region.area)
		perimeters.append(region.perimeter)
	areas = np.array(areas)
	perimeters = np.array(perimeters)
	area_hist, area_bins = np.histogram(areas, bins=50)
	area_bins = area_bins
	otsu = filters.threshold_otsu(areas) * 1.0
	for region in regions:
		if region.area > otsu:
			too_big.append(region.label)
	single_cells = labels.copy()
	for lab in too_big:
		single_cells[single_cells == lab] = 0
	single_cells, fw, inv = seg.relabel_sequential(single_cells)
	return single_cells

# find average fluorescence at normalized distances from the center of the cell
def find_distance_from_center(gfp, bf, signal_list):
	
	cells = find_cells(bf)
	blur = filters.gaussian(bf, 20)
	cells = remove_large_things(cells, blur)

	regions = measure.regionprops(cells)

	centers = np.array([r.centroid for r in regions])
	
	centers = np.round(centers).astype(np.int)
	seeds = np.ones(gfp.shape)
	seeds[centers[:,0], centers[:, 1]] = 0

	distances = ndi.distance_transform_edt(seeds)

	distances[cells==0] = 0

	for r in regions:
		mr, mc, xr, xc = r.bbox
		distances[mr:xr, mc:xc] = 9 * distances[mr:xr, mc:xc]/np.amax(distances[mr:xr, mc:xc])
		
	dist_range = 9

	x, y = (cells > 0).nonzero()

	for row, column in zip(x,y):
		index = int(math.floor(distances[row,column]))
		signal_list[index].append(gfp[row,column])


def process(file):

	stack = nd2r.Nd2('Images/' + file)

	for i in range(6):

		bf = normalize(stack[i*5+4])
		gfp = normalize(stack[i*5+2])

		find_distance_from_center(gfp, bf, signal)

	# with open('pickles/' + file + '.pkl', 'wb') as handle:
	# 	pik.dump(signal, handle)

	make_boxplot(signal, file)
	

def make_boxplot(data, name):

	plt.figure(figsize=(20,10), dpi=300)
	labels = np.linspace(10,100,10)
	plt.xlabel('% distance from center of cell')
	plt.ylabel('GFP signal')
	axes = plt.gca()
	axes.set_ylim([0,1])
	
	plt.boxplot(data)
	# plt.show()

	plt.savefig('boxplots/' + name + '.png')	
	plt.close()


group = re.compile('.*Well.01.*')
full = True

for filename in listdir('Images/'):
	if re.match(group, filename):
		signal = [[] for k in range(10)]
		
		if full:
			process(filename)

		else:
			with open('pickles/' + filename + '.pkl', 'rb') as handle:
				data = pik.load(handle)

			make_boxplot(data, filename)























	# if debug:
	# 	io.imsave("cells_on_bf.png", color.label2rgb(cells, image=bf, bg_label=0))
	# 	lbrgb = color.label2rgb(cells, image=gfp, bg_label=0)
	# 	r,c = np.where(~seeds.astype(np.bool))
	# 	lbrgb[r,c,:] = [1,0,0]
	# 	# lbrgb = (lbrgb * 255).astype(np.uint8)
	# 	io.imsave("cells_on_dapi_with_seeds.png", lbrgb)



	# for r in regions:
	# 	r.centroid  # Centroid tuple (row, col).
	# 	r.equivalent_diameter  # Equiv circle dia.
	# 	r.image  # Binary image in bounding box.
	# 	r.intensity_image  # Grayscale image in bounding box.
	# 	distances[r.coords[0], r.coords[1]]  # Vector of distances in the region.
	# 	distances[r.image]  # Same as above.
	# 	prc = np.percentile(distances[r.image], np.linspace(0,100,20))
	# 	print prc

	# x, y = (cells == 0).nonzero()

	# gfp_with_centers = gfp
	# gfp_no_bg = gfp

	# for row, column in zip(x, y):
	# 	gfp_no_bg[row, column] = -1

	# x, y = (centers == True).nonzero()

	# for row, column in zip(x, y):
	

	# x, y = (cells == 0).nonzero()

	# for row, column in zip(x, y):
	# 	distances[row, column] = 0

	# dist_range = int(math.floor(np.amax(distances[cells!=0])))






# cells = find_cells(test)

# x, y = (cells == 0).nonzero()

# gfp_no_bg = gfp_l

# for row, column in zip(x, y):
# 	gfp_no_bg[row, column] = 0

# distances = ndi.distance_transform_edt(gfp_no_bg)

# x, y = (cells > 0).nonzero()

# averages = [[] for k in range(int(math.floor(np.amax(distances))))]

# for row, column in zip(x,y):
# 	index = int(math.ceil(distances[row,column]))
# 	averages[index].append(gfp_l[row,column])









# ave = []

# for ls in averages:
# 	ave.append(np.mean(ls))

# plt.boxplot(averages)
# plt.show

# plt.plot(range(len(ave)),ave)
# plt.show()


# distances_flat = distances.flat
# gfp_flat = gfp_l.flat

# averages = [[]]*int(math.ceil(np.amax(distances_flat)))


# for number, distance in enumerate(distances_flat):
# 	print int(math.floor(distance))
# 	# averages[index].append(gfp_flat[number])

# print averages

# plt.plot(distances_flat, gfp_flat)
# plt.show()


# io.imshow(distances)
# io.show()

# 
# image = image.astype(np.float64)

# image2 = image - np.amin(image)
# image3 = image / np.amax(image2)



# image = filters.median(image, morphology.disk(2))

# def segment_image(photo_matrix, background_threshold, foreground_threshold):
# 	edges = sobel(photo_matrix)
# 	markers = np.zeros_like(photo_matrix)
# 	foreground, background = 1, 2
# 	markers[photo_matrix < background_threshold] = background
# 	markers[photo_matrix > foreground_threshold] = foreground
# 	ws = watershed(edges, markers)
# 	segmentation_matrix = ndi.label(ws == foreground)[0]
# 	return segmentation_matrix

# image = segment_image(image, 0.1, 0.9)

# io.imshow(image)
# io.show()




# import image and convert values to floats between 0 and 1 in a matrix
#photo_file = ndi.imread('images/C3-Practice.tif')
#img = sk.img_as_float(photo_file)

# med = filters.median(photo_matrix,morphology.disk(2))

#img = np.random.rand(1000,1000)

#img = filters.median(img, morphology.square(3))

#img = filters.laplace(img)

#io.imshow(img)
#io.show()