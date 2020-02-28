import numpy as np
import matplotlib.pyplot as plt

filename = "MC_composition.txt"
fl=open(filename,'r')
line= fl.readline();

temps=np.linspace(100, 1500, num=29)

xx=[]
line=fl.readline();

while len(line)>0:
    tmp=[]
    line=line.split()
    for i in range(len(line)):
        tmp.append(float(line[i]))
    xx.append(tmp)
    line=fl.readline()

fl.close()

#print len(mus)
print len(xx[:][1])
print len(xx)
#print (xx[:][1])
#print (xx[1][:])
#print(size(xx))

for j in range(21):
  #  tmp= [a[j] for a in xx]
    plt.plot(xx[j],temps)
plt.title("Phase Diagram (Vnn=30meV)\n100x100 Lattice")
plt.xlabel("X Composition")
plt.ylabel("Temperature")
plt.show()

