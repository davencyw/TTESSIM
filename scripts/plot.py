import csv
import seaborn
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from os import listdir
from os.path import isfile, join

folder = "../data/out/"
files = [f for f in listdir(folder) if f.endswith(".csv")]

for f in files :
	print("reading " + f)


	with open(join(folder,f), 'rb') as csvfile:
	    csvreader = csv.reader(csvfile, delimiter=';')
	    data = list(csvreader)

	    steps = int(len(data) - 1)
	    npdata = np.array(data[1:steps+1],dtype=np.float32)

	    numcells = float(data[0][0])
	    dt = float(data[0][1])
	    endrange = dt*steps
	    endrange = dt*steps
	    x = np.arange(0.0,endrange,dt)
	    
	    plt.plot(x,npdata)
	    plt.show()

