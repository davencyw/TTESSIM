import csv
import seaborn
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from os import listdir
from os.path import isfile, join

folder = "../data/testing/"
materials = ['f','s']
material_name = ["fluid","solid"]

for i in range(0,2) :
	ending = materials[i]+".csv"

	files = [f for f in listdir(folder) if f.endswith(ending)]
	for f in files :
		print("reading " + f)

		with open(join(folder,f), 'rb') as csvfile:
		    csvreader = csv.reader(csvfile, delimiter=';')
		    data = list(csvreader)
		    steps = int(len(data))/4
		    pecletnumbers = np.zeros(4)
		    print(steps)
		    for curve in range(0,4):
		    	npdata = np.transpose(np.array(data[curve*steps:(curve+1)*steps],dtype=np.float32))
		    	npdata = np.fliplr(npdata)
		        plt.loglog(npdata[0],npdata[2])
		        pecletnumbers[curve] = npdata[1][0]
		    print(pecletnumbers)
		    title = "OVS Plot for " + material_name[i] + " phase"
		    plt.xlim(npdata[0][0],npdata[0][-1])
		    plt.title(title)
		    plt.xlabel("N")
		    plt.ylabel("Error")
		    plt.legend(map(str, pecletnumbers))
		    plt.show()
