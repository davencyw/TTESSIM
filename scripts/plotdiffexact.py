import csv
import seaborn
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from os import listdir
from os.path import isfile, join

referencefile = "../data/reference/sol-exact-5.00000E+03.dat"

#read reference file
refdata = None
refmeshsize = 0
with open(referencefile, 'rb') as reffile:
	refreader = csv.reader(reffile, delimiter=";")
	data = list(refreader)
	refdata = np.transpose(np.array(data, dtype=np.float32))
	refmeshsize = len(refdata[0])


folder = "../data/out/"
materials = ['f','s']
materials_text = ["fluid","solid"]

for mat in range(0,2) :
	ending = materials[mat]+".csv"

	files = [f for f in listdir(folder) if f.endswith(ending)]

	for f in files :
		print("reading " + f)

		with open(join(folder,f), 'rb') as csvfile:
			##SORRY for that ugly code!
		    csvreader = csv.reader(csvfile, delimiter=';')
		    data = list(csvreader)
		    steps = int(len(data) - 1)
		    metadata = np.transpose(np.array(data[0],dtype=np.float32))
		    ops = metadata[2]
		    dt = metadata[1]
		    numcells = metadata[0]
		    print(metadata)
		    npdata = np.transpose(np.array(data[1:steps+1],dtype=np.float32))
		    diffcells = refmeshsize/numcells
		    newnpdata = refdata[mat+1,0::diffcells]
		    npdata = npdata[:,0]
		    print (newnpdata.shape)
		    print(npdata.shape)
		    newnpdata = abs(newnpdata-npdata)/newnpdata
		    print(newnpdata.shape)
		    print(np.max(newnpdata))
		    dt = float(data[0][1])
		    height = float(data[0][3])
		    dx = height / float(numcells)
		    endrange = dx * float(numcells) + (0.5 * dx)
		    x = np.arange(0.5*dx,endrange,dx)
		    plt.plot(x,newnpdata)
		    title = "error for " + materials_text[mat] +" phase at t = 5000s" + " with n = " + str(int(numcells))
		    plt.title(title)
		    plt.xlabel("x")
		    plt.ylabel("relative error")
		    plt.xlim([0.5*dx,endrange])
		    plt.show()
