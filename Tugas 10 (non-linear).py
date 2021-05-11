# -*- coding: utf-8 -*-
"""
Created on Sat Dec  5 20:25:08 2020

@group: Baha, Marsel, Rhesa
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

dataqe =[dataqe1,dataqe2,dataqe3,dataqe4]
dataCe = [dataCe1,dataCe2,dataCe3,dataCe4]


def langmuir (Ce,K,qm):
    return qm*K*Ce/(1+K*Ce)

def freundlich (Ce,K,n):
    return K*Ce**(1/n)

def DR (Ce,K,qs):
    R = 8.314*10**-3
    T = 35 +273.15
    e = R*T*np.log(1+1/Ce)
    return qs*np.exp(-K*e**2)

def reg (data1,data2,arg,i):
    if i == 0:
        calc = langmuir(data2,arg[0],arg[1])
    elif i == 1:
        calc = freundlich(data2,arg[0],arg[1])
    elif i == 2 :
        calc = np.zeros(len(data2))
        for k in range (len(data2)):
            calc[k] = DR(data2[k],arg[0],arg[1])
    avg = np.sum(calc)/len(calc)
    return (np.sum((data1-avg)**2))/(np.sum((data1-avg)**2)+np.sum((data1-calc)**2))
    
#Non-linear
#lagmuir
solvL1 = curve_fit(langmuir,dataCe1,dataqe1)
solvL2 = curve_fit(langmuir,dataCe2,dataqe2)
solvL3 = curve_fit(langmuir,dataCe3,dataqe3)
solvL4 = curve_fit(langmuir,dataCe4,dataqe4)

plt.figure(0,figsize=(9,5.5))
plt.plot(dataCe1,dataqe1,'.k',label="M1")
plt.plot(np.linspace(0,dataCe1[-1]+1,101),langmuir(np.linspace(0,dataCe1[-1]+1,101),solvL1[0][0],solvL1[0][1]),'k')
plt.plot(dataCe2,dataqe2,'.r',label="M2")
plt.plot(np.linspace(0,dataCe2[-1]+1,101),langmuir(np.linspace(0,dataCe2[-1]+1,101),solvL2[0][0],solvL2[0][1]),'r')
plt.plot(dataCe3,dataqe3,'.g',label='M3')
plt.plot(np.linspace(0,dataCe3[-1]+1,101),langmuir(np.linspace(0,dataCe3[-1]+1,101),solvL3[0][0],solvL3[0][1]),'g')
plt.plot(dataCe4,dataqe4,'.b',label='M4')
plt.plot(np.linspace(0,dataCe4[-1]+1,101),langmuir(np.linspace(0,dataCe4[-1]+1,101),solvL4[0][0],solvL4[0][1]),'b')
plt.ylabel("qe (mg P/g)")
plt.xlabel("Ce (mg P/L)")
plt.legend()

#Freundlich
solvf1 = curve_fit(freundlich,dataCe1,dataqe1)
solvf2 = curve_fit(freundlich,dataCe2,dataqe2)
solvf3 = curve_fit(freundlich,dataCe3,dataqe3)
solvf4 = curve_fit(freundlich,dataCe4,dataqe4)

plt.figure(1)
plt.plot(dataCe1,dataqe1,'.k',label="M1")
plt.plot(np.linspace(0,dataCe1[-1]+1,101),freundlich(np.linspace(0,dataCe1[-1]+1,101),solvf1[0][0],solvf1[0][1]),'k')
plt.plot(dataCe2,dataqe2,'.r',label="M2")
plt.plot(np.linspace(0,dataCe2[-1]+1,101),freundlich(np.linspace(0,dataCe2[-1]+1,101),solvf2[0][0],solvf2[0][1]),'r')
plt.plot(dataCe3,dataqe3,'.g',label="M3")
plt.plot(np.linspace(0,dataCe3[-1]+1,101),freundlich(np.linspace(0,dataCe3[-1]+1,101),solvf3[0][0],solvf3[0][1]),'g')
plt.plot(dataCe4,dataqe4,'.b',label="M4")
plt.plot(np.linspace(0,dataCe4[-1]+1,101),freundlich(np.linspace(0,dataCe4[-1]+1,101),solvf4[0][0],solvf4[0][1]),'b')
plt.ylabel("qe (mg P/g)")
plt.xlabel("Ce (mg P/L)")
plt.legend()

#Dubinin-Radushkevich
solvdr1 = curve_fit(DR,dataCe1,dataqe1)
solvdr2 = curve_fit(DR,dataCe2,dataqe2)
solvdr3 = curve_fit(DR,dataCe3,dataqe3)
solvdr4 = curve_fit(DR,dataCe4,dataqe4)

plt.figure(2)
plt.plot(dataCe1,dataqe1,'.k',label="M1")
plt.plot(np.linspace(0.01,dataCe1[-1]+1,101),DR(np.linspace(0.01,dataCe1[-1]+1,101),solvdr1[0][0],solvdr1[0][1]),'k')
plt.plot(dataCe2,dataqe2,'.r',label="M2")
plt.plot(np.linspace(0.01,dataCe2[-1]+1,101),DR(np.linspace(0.01,dataCe2[-1]+1,101),solvdr2[0][0],solvdr2[0][1]),'r')
plt.plot(dataCe3,dataqe3,'.g',label="M3")
plt.plot(np.linspace(0.01,dataCe3[-1]+1,101),DR(np.linspace(0.01,dataCe3[-1]+1,101),solvdr3[0][0],solvdr3[0][1]),'g')
plt.plot(dataCe4,dataqe4,'.b',label="M4")
plt.plot(np.linspace(0.01,dataCe4[-1]+2,101),DR(np.linspace(0.01,dataCe4[-1]+2,101),solvdr4[0][0],solvdr4[0][1]),'b')
plt.legend()

#Perhitungan regresi

#Langmuir
regL = np.zeros(len(dataqe))
solvL = [solvL1[0],solvL2[0],solvL3[0],solvL4[0]]

#Freundlich
regF = np.zeros(len(dataqe))
solvF = [solvf1[0],solvf2[0],solvf3[0],solvf4[0]]

#Dubinin-Radushkevich
regDR = np.zeros(len(dataqe))
solvDR = [solvdr1[0],solvdr2[0],solvdr3[0],solvdr4[0]]

for i in range (len(regL)):
    regL [i] = reg(dataqe[i],dataCe[i],solvL[i],0)
    regF [i] = reg(dataqe[i],dataCe[i],solvF[i],1)
    regDR [i] = reg(dataqe[i],dataCe[i],solvDR[i],2)

#Hasil
print("M1 :")
print("Langmuir; Kl = %.3f (L/mg) ,qm = %.3f (mg P/g), regresi = %.4f"%(solvL1[0][0],solvL1[0][1],regL[0]))
print("")
print("Freundlich; Kf = %.3f (L/mg) ,n = %.3f (mg P/g), regresi = %.4f"%(solvf1[0][0],solvf1[0][1],regF[0]))
print("")
print("D-R; Kd = %.3f (L/mg) ,qs = %.3f (mg P/g), regresi = %.4f"%(solvdr1[0][0],solvdr1[0][1],regDR[0]))
print("")

print("M2 :")
print("Langmuir; Kl = %.3f (L/mg) ,qm = %.3f (mg P/g), regresi = %.4f"%(solvL2[0][0],solvL2[0][1],regL[1]))
print("")
print("Freundlich; Kf = %.3f (L/mg) ,n = %.3f (mg P/g), regresi = %.4f"%(solvf2[0][0],solvf2[0][1],regF[1]))
print("")
print("D-R; Kd = %.3f (L/mg) ,qs = %.3f (mg P/g), regresi = %.4f"%(solvdr2[0][0],solvdr2[0][1],regDR[1]))
print("")

print("M3 :")
print("Langmuir; Kl = %.3f (L/mg) ,qm = %.3f (mg P/g), regresi = %.4f"%(solvL3[0][0],solvL3[0][1],regL[2]))
print("")
print("Freundlich; Kf = %.3f (L/mg) ,n = %.3f (mg P/g), regresi = %.4f"%(solvf3[0][0],solvf3[0][1],regF[2]))
print("")
print("D-R; Kd = %.3f (L/mg) ,qs = %.3f (mg P/g), regresi = %.4f"%(solvdr3[0][0],solvdr3[0][1],regDR[2]))
print("")

print("M4 :")
print("Langmuir; Kl = %.3f (L/mg) ,qm = %.3f (mg P/g), regresi = %.4f"%(solvL4[0][0],solvL4[0][1],regL[3]))
print("")
print("Freundlich; Kf = %.3f (L/mg) ,n = %.3f (mg P/g), regresi = %.4f"%(solvf4[0][0],solvf4[0][1],regF[3]))
print("")
print("D-R; Kd = %.3f (L/mg) ,qs = %.3f (mg P/g), regresi = %.4f"%(solvdr4[0][0],solvdr4[0][1],regDR[3]))
print("")