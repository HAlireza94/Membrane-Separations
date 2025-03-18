import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

numC = 1008
ionPairs = 11
numWater = 1465
totalParticles = numC + 2*ionPairs + numWater*3
id, typelabel, mol, type, q, x, y, z = [], [], [], [], [], [], [], []
with open('dump.water_prd.lammpstrj') as f:
    for line in f:
        p = line.split()
        if len(p) == 8 and p[0] != 'ITEM:':
            id.append(p[0])
            typelabel.append(p[1])
            mol.append(p[2])
            type.append(p[3])
            q.append(p[4])
            x.append(float(p[5]))
            y.append(p[6])
            z.append(p[7])


df = pd.DataFrame({'id':id, 'mol':mol, 'q':q, 'x':x, 'y': y, 'z': z})
numFrames = int(len(df['id'])/totalParticles)

# Membrane is in X direction, membrane is frozen, then one frame is enough to get boundaries
LS_min, LS_max = 0, 0
for i in range(1):
    frame = df['x'][i*totalParticles:(i+1)*totalParticles]
    LS_min = np.min((frame[:numC]))
    LS_max = np.max((frame[:numC]))
# water displacement in bulk
numwaterBulk, numwaterFeed = [], []
for i in range(0,numFrames):
    data1 = df[i*totalParticles:(i+1)*totalParticles].reset_index(drop=True)
    data2 = data1[numC:numC+numWater*3]
    data = data2['x'].reset_index(drop=True)
    c_b, c_f = 0, 0
    for j in range(0,len(data),3):
        if data[j] > LS_max-1:
            c_b += 1    
    for k in range(0,len(data),3):
        if data[k] < LS_min+1:
            c_f += 1
    numwaterBulk.append(c_b)
    numwaterFeed.append(c_f)

time = [(i) * 250/1000 for i in range(numFrames)]  

print(f"{'Time(ns)':<15}{'Bulk':<15}{'Feed':<15}")
for i in range(len(time)):
    print(f"{time[i]:<15}{int(numwaterBulk[i]):<15}{int(numwaterFeed[i]):<15}")




#plt.plot(time, numwaterBulk, '.k')
#plt.plot(time, numwaterFeed, '+b')
#plt.legend(['Bulk','Feed'],frameon=False)
#plt.xlabel('Time (ns)')
#plt.ylabel('Number of Water Molecules')
#plt.savefig('Equilibrium')
