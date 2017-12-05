import csv
import seaborn
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from os import listdir
from os.path import isfile, join

folder = "../data/out/"
materials = ['f','s']

for mat in range(0,2) :
	ending = materials[mat]+".csv"

	files = [f for f in listdir(folder) if f.endswith(ending)]

	for f in files :
		print("reading " + f)

		with open(join(folder,f), 'rb') as csvfile:
		    csvreader = csv.reader(csvfile, delimiter=';')
		    data = list(csvreader)
		    steps = int(len(data) - 1)
		    metadata = np.transpose(np.array(data[0],dtype=np.float32))
		    ops = metadata[2]
		    dt = metadata[1]
		    numcells = metadata[0]
		    print(ops)
		    print(metadata)
		    npdata = np.transpose(np.array(data[1:steps+1],dtype=np.float32))
		    for step in range (0,steps) :
			    dt = float(data[0][1])
			    height = float(data[0][3])
			    dx = height / float(numcells)
			    endrange = dx * float(numcells) + (0.5 * dx)
			    x = np.arange(0.5*dx,endrange,dx)


			    title = mat + " TTESSIM at t=" + (str(ops*dt*(step+1)))
			    plt.plot(x,npdata[:,step])
			    plt.xlim([0.5*dx,endrange])
			    plt.title(title)
			    plt.show()
