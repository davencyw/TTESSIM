import csv
import seaborn
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from os import listdir
from os.path import isfile, join

folder = "../data/testing/"
materials = ['f','s']

for mat in range(0,2) :
	ending = materials[mat]+".csv"

	files = [f for f in listdir(folder) if f.endswith(ending)]
	for f in files :
		print("reading " + f)

		with open(join(folder,f), 'rb') as csvfile:
		    csvreader = csv.reader(csvfile, delimiter=';')
		    data = list(csvreader)
		    steps = int(len(data))
		    print(steps)
		    npdata = np.transpose(np.array(data[0:steps],dtype=np.float32))
		    npdata = np.fliplr(npdata)
		    print(npdata[0])
		    print(npdata[1])
			#title = "TTESSIM at t=" + (str(ops*dt*(step+1)))
		    plt.loglog(npdata[0],npdata[1])
		    plt.xlim(npdata[0][0],npdata[0][-1])
		    #plt.title(title)
		    plt.xlabel("N")
		    plt.ylabel("Error")
		    plt.show()
