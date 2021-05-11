# -*- coding: utf-8 -*-
"""
Created on Sun Dec  6 14:09:42 2020

@group: Baha, Marsel, and Rhesa
"""
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_excel(r'Data adsobsi.xlsx',sheet_name=0,header=None, index_col=None)
data = np.array(data)
dataCo = data[2:,0]
#M1
dataCe1 = data[2:,1]
dataqe1 = data[2:,2]
#M2
dataCe2 = data[2:,3]
dataqe2 = data[2:,4]
#M3
dataCe3 = data[2:,5]
dataqe3 = data[2:,6]
#M4
dataCe4 = data[2:,7]
dataqe4 = data[2:,8]

dataqe = np.array([dataqe1,dataqe2,dataqe3,dataqe4])
dataCe = np.array([dataCe1,dataCe2,dataCe3,dataCe4])

#data Langmuir
datayl = np.array([dataCe1/dataqe1,dataCe2/dataqe2,dataCe3/dataqe3,dataCe4/dataqe4])
dataxl = np.array([dataCe1,dataCe2,dataCe3,dataCe4])


#data Freundlich
datayf = np.zeros((4,6))
dataxf = np.zeros((4,6))

for j in range (4):
    for i in range (6):
        datayf[j][i]=np.log(dataqe[j][i])
        dataxf [j][i] = np.log(dataCe[j][i])
        
#data Dubinin-Radushkevich
R = 8.314*10**-3
T = 35 +273.15
dataydr = np.zeros((4,6))
dataxdr = np.zeros((4,6))
for jj in range (4):
    for ii in range (6):
        dataydr[jj][ii] = np.log(dataqe[jj][ii])
        e = R*T*np.log(1+1/dataCe[jj][ii])
        dataxdr[jj][ii] = e**2
      

def linL (x,K,qm):
    return 1/(K*qm)+x/qm

def linF (x,K,n):
    return np.log(K)+x/n

def linDR (x,K,qs):
    return np.log(qs)-K*x

def reg (data1,data2,arg,i):
    if i == 0:
        calc = linL(data2,arg[0],arg[1])
    elif i == 1:
        calc = linF(data2,arg[0],arg[1])
    elif i == 2 :
        calc = linDR(data2,arg[0],arg[1])
    avg = np.sum(calc)/len(calc)
    return (np.sum((data1-avg)**2))/(np.sum((data1-avg)**2)+np.sum((data1-calc)**2))

#Linear Langmuir
solvl1 = curve_fit(linL,dataxl[0,1:],datayl[0,1:],bounds=(0,[0.5,6.5]))
solvl2 = curve_fit(linL,dataxl[1,1:],datayl[1,1:])
solvl3 = curve_fit(linL,dataxl[2,1:],datayl[2,1:],bounds=(1,[2.3,46]))
solvl4 = curve_fit(linL,dataxl[3,1:],datayl[3,1:])

plt.figure(0)
plt.plot(dataxl[0,:],datayl[0,:],'.k',label='M1')
plt.plot(np.linspace(0,dataxl[0][-1]+1,101),linL(np.linspace(0,dataxl[0][-1]+1,101),solvl1[0][0],solvl1[0][1]),'k')
plt.plot(dataxl[1,:],datayl[1,:],'.r',label='M2')
plt.plot(np.linspace(0,dataxl[0][-1]+1,101),linL(np.linspace(0,dataxl[1][-1]+1,101),solvl2[0][0],solvl2[0][1]),'r')
plt.plot(dataxl[2,:],datayl[2,:],'.g',label='M3')
plt.plot(np.linspace(0,dataxl[0][-1]+3,101), linL(np.linspace(0,dataxl[2][-1]+3,101),solvl3[0][0],solvl3[0][1]),'g')
plt.plot(dataxl[3,:],datayl[3,:],'.b',label='M4')
plt.plot(np.linspace(0,dataxl[0][-1]+1,101),linL(np.linspace(0,dataxl[3][-1]+1,101),solvl4[0][0],solvl4[0][1]),'b')
plt.ylim(0,12)
plt.xlim(0,60)
plt.legend()
plt.xlabel("Ce/qe (g/L)")
plt.ylabel("Ce (mg P/L)")

#Linear Freundlich
solvf1 = curve_fit(linF,dataxf[0,1:],datayf[0,1:],bounds=(0.,[1.,2.]))
solvf2 = curve_fit(linF,dataxf[1,:],datayf[1,:],bounds=(1.,[7,2.7]))
solvf3 = curve_fit(linF,dataxf[2,:],datayf[2,:])
solvf4 = curve_fit(linF,dataxf[3,:],datayf[3,:])

plt.figure(1)
plt.plot(dataxf[0,:],datayf[0,:],'.k',label='M1')
plt.plot(np.linspace(dataxf[0][0]-0.5,dataxf[0][-1]+1,101),linF(np.linspace(dataxf[0][0]-0.5,dataxf[0][-1]+1,101),solvf1[0][0],solvf1[0][1]),'k')
plt.plot(dataxf[1,:],datayf[1,:],'.r',label='M2')
plt.plot(np.linspace(dataxf[1][0]-1,dataxf[1][-1]+1,101),linF(np.linspace(dataxf[1][0]-1,dataxf[1][-1]+1,101),solvf2[0][0],solvf2[0][1]),'r')
plt.plot(dataxf[2,:],datayf[2,:],'.g',label='M3')
plt.plot(np.linspace(dataxf[2][0]-1,dataxf[2][-1]+1,101),linF(np.linspace(dataxf[2][0]-1,dataxf[2][-1]+1,101),solvf3[0][0],solvf3[0][1]),'g')
plt.plot(dataxf[3,:],datayf[3,:],'.b',label='M4')
plt.plot(np.linspace(dataxf[3][0]-1,dataxf[3][-1]+1,101),linF(np.linspace(dataxf[3][0]-1,dataxf[3][-1]+1,101),solvf4[0][0],solvf4[0][1]),'b')
plt.legend()
plt.ylabel("ln(qe)")
plt.xlabel("ln(Ce)")

#Linear Dubinin-Radushkevich
solvDR1 = curve_fit(linDR,dataxdr[0,:],datayf[0,:],bounds=(0.,[1.,3.9]))
solvDR2 = curve_fit(linDR,dataxdr[1,:],datayf[1,:])
solvDR3 = curve_fit(linDR,dataxdr[2,:],datayf[2,:],bounds=(0.,[1,45]))
solvDR4 = curve_fit(linDR,dataxdr[3,:],datayf[3,:])

plt.figure(2)
plt.plot(dataxdr[0,:],dataydr[0,:],'.k',label='M1')
plt.plot(np.linspace(dataxdr[0][0],dataxdr[0][-1],101),linDR(np.linspace(dataxdr[0][0],dataxdr[0][-1],101),solvDR1[0][0],solvDR1[0][1]),'k')
plt.plot(dataxdr[1,:],dataydr[1,:],'.r',label='M2')
plt.plot(np.linspace(0,100,101),linDR(np.linspace(0,100,101),solvDR2[0][0],solvDR2[0][1]),'r')
plt.plot(dataxdr[2,:],dataydr[2,:],'.g',label='M3')
plt.plot(np.linspace(0,125,101),linDR(np.linspace(0,125,101),solvDR3[0][0],solvDR3[0][1]),'g')
plt.plot(dataxdr[3,:],dataydr[3,:],'.b',label='M4')
plt.plot(np.linspace(0,125,101),linDR(np.linspace(0,125,101),solvDR4[0][0],solvDR4[0][1]),'b')
plt.xlabel("epsilon^2 (kJ^2/mol^2)")
plt.ylabel("ln(qe)")
plt.legend()

#Perhitungan regresi

#Langmuir
regL = np.zeros(len(dataqe))
solvL = [solvl1[0],solvl2[0],solvl3[0],solvl4[0]]

#Freundlich
regF = np.zeros(len(dataqe))
solvF = [solvf1[0],solvf2[0],solvf3[0],solvf4[0]]

#Dubinin-Radushkevich
regDR = np.zeros(len(dataqe))
solvDR = [solvDR1[0],solvDR2[0],solvDR3[0],solvDR4[0]]

for i in range (len(regL)):
    regL [i] = reg(datayl[i,1:],dataxl[i,1:],solvL[i],0)
    regF [i] = reg(datayf[i,1:],dataxf[i,1:],solvF[i],1)
    regDR [i] = reg(dataydr[i,1:],dataxdr[i,1:],solvDR[i],2)

#Hasil
print("M1 :")
print("Langmuir; Kl = %.3f (L/mg) ,qm = %.3f (mg P/g), regresi = %.4f"%(solvl1[0][0],solvl1[0][1],regL[0]))
print("")
print("Freundlich; Kf = %.3f (L/mg) ,n = %.3f (mg P/g), regresi = %.4f"%(solvf1[0][0],solvf1[0][1],regF[0]))
print("")
print("D-R; Kd = %.3f (L/mg) ,qs = %.3f (mg P/g), regresi = %.4f"%(solvDR1[0][0],solvDR1[0][1],regDR[0]))
print("")

print("M2 :")
print("Langmuir; Kl = %.3f (L/mg) ,qm = %.3f (mg P/g), regresi = %.4f"%(solvl2[0][0],solvl2[0][1],regL[1]))
print("")
print("Freundlich; Kf = %.3f (L/mg) ,n = %.3f (mg P/g), regresi = %.4f"%(solvl2[0][0],solvf2[0][1],regF[1]))
print("")
print("D-R; Kd = %.3f (L/mg) ,qs = %.3f (mg P/g), regresi = %.4f"%(solvDR2[0][0],solvDR2[0][1],regDR[1]))
print("")

print("M3 :")
print("Langmuir; Kl = %.3f (L/mg) ,qm = %.3f (mg P/g), regresi = %.4f"%(solvl3[0][0],solvl3[0][1],regL[2]))
print("")
print("Freundlich; Kf = %.3f (L/mg) ,n = %.3f (mg P/g), regresi = %.4f"%(solvl3[0][0],solvf3[0][1],regF[2]))
print("")
print("D-R; Kd = %.3f (L/mg) ,qs = %.3f (mg P/g), regresi = %.4f"%(solvDR3[0][0],solvDR3[0][1],regDR[2]))
print("")


print("M4 :")
print("Langmuir; Kl = %.3f (L/mg) ,qm = %.3f (mg P/g), regresi = %.4f"%(solvl4[0][0],solvl4[0][1],regL[3]))
print("")
print("Freundlich; Kf = %.3f (L/mg) ,n = %.3f (mg P/g), regresi = %.4f"%(solvl4[0][0],solvf4[0][1],regF[3]))
print("")
print("D-R; Kd = %.3f (L/mg) ,qs = %.3f (mg P/g), regresi = %.4f"%(solvDR4[0][0],solvDR4[0][1],regDR[3]))
print("")
