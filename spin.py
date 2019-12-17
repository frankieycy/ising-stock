import imageio,glob,re
import numpy as np
import scipy.misc as misc

# generate spin evolution with spin data files

def file_num(f):
	return int(re.split('[_ .]',f)[1])

# make images

files = glob.glob('out/spin_*.csv') # get csv files

for f in files:
	d = np.loadtxt(f,delimiter=',')
	img = misc.toimage(d).resize((256,256)) # convert to img
	img.save('plt/'+f[4:-4]+'.png')

# create gif ani

files = glob.glob('plt/spin_*.png') # get png files
images = []

for f in sorted(files,key=file_num):
	images.append(imageio.imread(f))

imageio.mimsave('plt/ani_spin.gif', images, duration=0.5)
