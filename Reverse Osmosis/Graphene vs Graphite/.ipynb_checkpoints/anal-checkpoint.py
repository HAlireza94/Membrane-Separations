import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lmfit import Model
AA, BB, alpha1, alpha2 = 1.7937, 1.5088, 0.10009, 1.1422 * 0.01
a00, a11, a22, a33 = 0.99654, 5.0615*0.01, -5.3907*0.001, -6.3600 * 0.0001

numC = 1008
numRP = 336
ionPairs = 11
numWater = 1465
Mw_water = 18.01528
Mw_ions = 58.44
v = 2
R = 8.314462618 # j/K.mol
T = 300
totalParticles = numC + 2*ionPairs + numWater*3 + numRP

id, typelabel, mol, type, q, x, y, z = [], [], [], [], [], [], [], []
with open('dump.water_eql.lammpstrj') as f:
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


df = pd.DataFrame({'id':id, 'label':typelabel,'mol':mol,'type':type,'q':q, 'x':x, 'y': y, 'z': z})
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
    data2 = data1[numC+numRP:numC+numRP+numWater*3]
    data = data2['x'].reset_index(drop=True)
    c_b, c_f = 0, 0
    for j in range(0,len(data),3):
        if data[j] > LS_max:
            c_b += 1    
    for k in range(0,len(data),3):
        if data[k] < LS_min:
            c_f += 1
    numwaterBulk.append(c_b)
    numwaterFeed.append(c_f)

time = [(i) * 500/1000 for i in range(numFrames)]  
molality = [(1000*ionPairs)/(Mw_water*i) for i in numwaterFeed]





def water_density_p(m,a11,a22,a33):
    
    rho = a00 + a11*m + a22*m**1.5 + a33*m**2
    
    return rho


def d_rho_p(a11,a22,a33,m):
    
    drho = a11 * m + (3/2)*a22*m**1.5 + 2*a33*m**2
    
    return drho

rho_p = []
drho_p = []

for i in range(len(molality)):
    rho_p.append(water_density_p(molality[i],a11,a22,a33))
    drho_p.append(d_rho_p(a11,a22,a33,molality[i]))
    
partial_molar_volume = []
for i in range(len(rho_p)):
    partial_molar_volume.append((Mw_water/rho_p[i]) + (drho_p[i]/rho_p[i]**2) *
                                (Mw_water + (((Mw_ions*Mw_water) / 1000) * molality[i])))




## Fitting for Partial Molar Volume
def pmv(m,gamma0,gamma1,gamma2,gamma3):
#     y = gamma0 + gamma1*m + gamma2*m**1.5 + gamma3*m**2
    y = gamma0 + gamma1*m + gamma2*m**2 + gamma3*m**3
    
    return y

model = Model(pmv)
params4 = model.make_params(gamma0 = 0.01, gamma1 = 0.1, gamma2 = 0.02, gamma3 = 0.2,)
result4 = model.fit(partial_molar_volume, params4, m=molality)

g0 = np.array(result4.params)[0]
g1 = np.array(result4.params)[1]
g2 = np.array(result4.params)[2]
g3 = np.array(result4.params)[3]
PMV = []
for i in range(len(molality)):
    PMV.append(pmv(molality[i],g0,g1,g2,g3))


## Osmotic Pressure Molal Description

def osmotic_pressure_molal(m,alpha1,alpha2,AA,BB,partial_mv):
    
    u1 = (10*v*R*T*0.001*(Mw_water/partial_mv))
    u2 = -((AA/BB**2) * ((2 * np.sqrt(m) + BB * m) / (1 + BB * np.sqrt(m)))) 
    u3 = ((2*AA/BB**3)*np.log(1 + BB * np.sqrt(m)))
    u4 = m + ((1/2)*alpha1*m**2) + ((2/3)*alpha2*m**3)
    
    return u1*(u4+u2+u3)
    

os_pressur_molal = []
for i in range(len(molality)):
    os_pressur_molal.append(osmotic_pressure_molal(molality[i],alpha1,
                                                   alpha2,AA,BB,pmv(molality[i],g0,g1,g2,g3)))


## Mean ionic activity

def mean_ionic_activity(alpha1,alpha2,AA,BB,m):
    u1 = alpha1*m + alpha2*m**2
    u2 = (AA*m**0.5)/(1+(BB*m**0.5))
    
    return u1 - u2 
mia = []
for i in range(len(molality)):
    mia.append(mean_ionic_activity(alpha1,alpha2,AA,BB,molality[i]))



## activity of water

def water_activity(m,alpha1,alpha2,AA,BB):
    
    u1 = -v*R*T*0.001*(Mw_water/1000)
    u2 = -((AA/BB**2) * ((2 * np.sqrt(m) + BB * m) / (1 + BB * np.sqrt(m)))) 
    u3 = ((2*AA/BB**3)*np.log(1 + BB * np.sqrt(m)))
    u4 = m + ((1/2)*alpha1*m**2) + ((2/3)*alpha2*m**3)
    
    return u1*(u4+u2+u3)

aw = []
sqrt_m = []
for i in range(len(molality)):
    aw.append(water_activity(molality[i],alpha1,alpha2,AA,BB))
    sqrt_m.append(molality[i]**0.5)

ms_exp, rhos_exp, phis_exp, lngamma_exp, os_exp = [], [], [], [], []
with open('exeprimental data') as f:
    for line in f:
        p = line.split()
        ms_exp.append(float(p[0]))
        rhos_exp.append(float(p[1]))
        phis_exp.append(float(p[2]))
        lngamma_exp.append(np.log(float(p[3])))
        os_exp.append(float(p[4]))




### experimental fitting ############

def RHOS_Fit_ex(m, a_ex00, a_ex11, a_ex22, a_ex33):
    rho = a_ex00 + a_ex11 * m + a_ex22 * m**1.5 + a_ex33 * m**2 
    return rho


model5 = Model(RHOS_Fit_ex)
params5 = model5.make_params(a_ex00=0.997, a_ex11=0.1, a_ex22=0.02, a_ex33=0.2)
result5 = model5.fit(rhos_exp, params5, m=ms_exp)

rho_P_exp = [RHOS_Fit_ex(i,np.array(result5.params)[0],np.array(result5.params)[1],np.array(result5.params)[2],np.array(result5.params)[3]) for i in molality]
def os_Fit_ex(m,p1, p2):
    os_fit = p1*m**(1) + p2*m**2
    return os_fit


model6 = Model(os_Fit_ex)
params6 = model6.make_params(p1=1,p2=1)
result6 = model6.fit(os_exp, params6, m=ms_exp)
fit_os_exp = [os_Fit_ex(i,np.array(result6.params)[0], np.array(result6.params)[1]) for i in molality]


c = 0
titles,mss_exp, phi_exp, lnaw_exp, RTlnaw_exp = [], [], [], [], []
with open('activity') as f:
    for line in f:
        p = line.split()
        c += 1
        if c == 1:
            titles.append(p)
        else:
            mss_exp.append(float(p[0]))
            phi_exp.append(float(p[1]))
            lnaw_exp.append(float(p[2]))
            RTlnaw_exp.append(float(p[3]))


mss_exp = [np.sqrt(i) for i in mss_exp]

def aw_fit_exp(m,p1, p2,p3):
    aw_fit = p1*m**1 + p2*m**2 + p3*m**3
    return aw_fit




model7 = Model(aw_fit_exp)
params7 = model7.make_params(p1=1,p2=1,p3=1)
result7 = model7.fit(RTlnaw_exp, params7, m=mss_exp)
fit_aw_exp = [aw_fit_exp(i,np.array(result7.params)[0], np.array(result7.params)[1],np.array(result7.params)[2]) for i in molality]

# plt.plot(mss_exp,RTlnaw_exp,'or')
# plt.plot(mss_exp,[aw_fit_exp(i,np.array(result7.params)[0], np.array(result7.params)[1],np.array(result7.params)[2]) for i in mss_exp],'--k')
# print(result7.fit_report())

print(f"{'Time(ns)':<10}{'Bulk':<10}{'Feed':<10}{'molality':<10}{'ρs(g/L)':<10}{'ρs(g/L)_ex':<13}{'V(L/mol)':<10}{'fit_V(L/mol)':<14}{'ln(γs)':<13}{'RTTLn(aw)':<13}{'RTTLn(aw)_exp':<15}{'Π(bar)':<12}{'Π(bar)_exp':<10}")



for i in range(len(time)):
    print(f"{time[i]:<10}{numwaterBulk[i]:<10}{numwaterFeed[i]:<10}{molality[i]:<10.5f}{rho_p[i]:<10.5f}{rho_P_exp[i]:<13.5f}{partial_molar_volume[i]:<10.5f}{PMV[i]:<14.5f}{mia[i]:<13.5f}{aw[i]:<13.5f}{fit_aw_exp[i]:<15.5f}{os_pressur_molal[i]:<12.5f}{fit_os_exp[i]:<10.5f}")


