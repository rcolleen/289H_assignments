import numpy as np
import matplotlib.pyplot as plt

filename = "MC_composition.txt"
fl=open(filename,'r')
line= fl.readline();

mus=np.linspace(-.5, 0.5, num=21)

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

print len(mus)
print len(xx[:][1])
print len(xx)
#print (xx[:][1])
#print (xx[1][:])
#print(size(xx))

for j in range(len(xx[1][:])):
    tmp= [a[j] for a in xx]
    plt.plot(tmp, mus)
plt.title("Chemical Potential vs X\n 100x100 grid with Vnn=30meV")
plt.ylabel("Chemical Potential meV")
plt.xlabel("COmposition X")
plt.show()

