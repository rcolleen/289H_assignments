import numpy as np
import matplotlib.pyplot as plt

filename = "MC_cv.txt"
fl=open(filename,'r')
line= fl.readline();

temps=np.linspace(100,1500,num=29)

cv=[]
line=fl.readline();

while len(line)>0:
    tmp=[]
    line=line.split()
    for i in range(len(line)):
        tmp.append(float(line[i]))
    cv.append(tmp)
    line=fl.readline()

fl.close()
for j in range(21):
#    tmp=[a[j] for a in cv]
    plt.plot(temps, cv[j])
plt.figure(2)
plt.title("C_v at Chemical Potential=0.0meV")
plt.xlabel("Temperature")
plt.ylabel("Heat Capacity")
plt.plot(temps, cv[11])
plt.show()


