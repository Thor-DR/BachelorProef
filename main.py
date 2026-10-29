import math

import numpy as np
import matplotlib.pyplot as plt
from classes import Detector, LAL_Detector
from classes  import DetectorPair
from classes import DetectorPairVec
from classes import DetectorNetwork
from time import time

def interpolate_asd(interpolation_frequencies,frequencies,asd):
    interpolated_asd = np.zeros(np.size(interpolation_frequencies))
    reversed_frequencies=np.flip(frequencies)
    for i in range(np.size(interpolation_frequencies)):
        f_lower=frequencies[1]
        index_lower=1
        for j in range(np.size(reversed_frequencies)):
            if reversed_frequencies[j]<interpolation_frequencies[i]:
                f_lower= reversed_frequencies[j]
                index_lower=np.size(reversed_frequencies)+1-j
                break
        f_higher = frequencies[-1]
        index_higher = -1
        for k in range(np.size(frequencies)):
            if frequencies[k] > interpolation_frequencies[i]:
                f_higher = frequencies[k]
                index_higher=k
                break
        f_l_diff=interpolation_frequencies[i]-f_lower
        f_h_diff=f_higher-interpolation_frequencies[i]
        interpolated_asd[i]=asd[index_lower]+f_h_diff*(asd[index_higher]-asd[index_lower])/(f_l_diff+f_h_diff)
    return interpolated_asd

def network_plot(network,freqs,det_respons, x_min=10,x_max=1500, fref=100):
    Omega_0 = network.get_Omega_eff(0)
    Omega_1 = network.get_Omega_eff(1)
    PI_curve, Omega_Bs = network.Power_Integration(betas)
    plt.loglog(freqs, det_respons, color='black')
    plt.loglog(freqs, Omega_0, color='red')
    plt.loglog(freqs, Omega_1, color='green')
    for beta, Omega_B in zip(betas, Omega_Bs):
        plt.loglog(freqs, Omega_B * np.power(freqs / fref, beta), color="gray")
    plt.loglog(freqs, PI_curve, color='blue')
    plt.xlim(x_min, x_max)
    plt.ylim(1e-10, 1e0)
    plt.show()

H0=2.2*10**(-18)
start = time()
starting_frequency=10
ending_frequency=1500
df=0.25
df_inv=int(1/df)
RT_frequencies=np.linspace(starting_frequency, ending_frequency, (ending_frequency - starting_frequency) * df_inv)
# aligo_design ="C:/Users/Thor De Rouck/Documents/0FysicaUAntwerpen/BachelorProef/aligo.txt"
# f_aligo_design, asd_aligo_design = np.transpose(np.loadtxt(aligo_design))
# RT_asd=interpolate_asd(RT_frequencies,f_aligo_design,asd_aligo_design)
# RT_psd=RT_asd**2
# Omega_RT= (2 * np.pi ** 2 * RT_frequencies ** 3 * RT_psd) / (3 * H0 ** 2)
ET_10km="C:/Users/Thor De Rouck/Documents/0FysicaUAntwerpen/BachelorProef/18213_ET10kmcolumns.txt"
ET_15km="C:/Users/Thor De Rouck/Documents/0FysicaUAntwerpen/BachelorProef/18213_ET15kmcolumns.txt"
ET_20km="C:/Users/Thor De Rouck/Documents/0FysicaUAntwerpen/BachelorProef/18213_ET20kmcolumns.txt"
f_ET_10,psd_ET_10,psd_ET_101,psd_ET_102=np.transpose(np.loadtxt(ET_10km))
f_ET_15,psd_ET_15,psd_ET_151,psd_ET_152=np.transpose(np.loadtxt(ET_15km))
f_ET_20,psd_ET_20,psd_ET_201,psd_ET_202=np.transpose(np.loadtxt(ET_20km))
RT_ET_10_psd=interpolate_asd(RT_frequencies,f_ET_10,psd_ET_10)
RT_ET_15_psd=interpolate_asd(RT_frequencies,f_ET_10,psd_ET_15)
RT_ET_20_psd=interpolate_asd(RT_frequencies,f_ET_10,psd_ET_20)
end = time()
print(end-start, "s")
start = time()
# u1=np.array([-0.22389266154, 0.79983062746, 0.55690487831])
# v1=np.array([-0.91397805644, 0.02609403989, -0.40492342125])
# x1=np.array([-2161414.92636, -3834695.17889, 4600350.22664])
# LIGO1=Detector(x1,u1,v1,RT_psd)
#
# u2=np.array([-0.95457412153, -0.14158077340, -0.26218911324])
# v2=np.array([0.29774156894, -0.48791033647, -0.82054461286])
# x2=np.array([-74276.04238, -5496283.71971, 3224257.01744])
# LIGO2=Detector(x2,u2,v2,RT_psd)
# LIGO_detectors=[LIGO1,LIGO2]

#coordinates
ET_S_long=math.radians(9.417)
ET_S_lat=math.radians(40.52)
ET_x_az=math.radians(90)
ET_y_az=math.radians(0)
ET_S_elev=300
ET_MR_long=math.radians(5.921)
ET_MR_lat=math.radians(50.72)
ET_x_az22=math.radians(67.5)
ET_y_az22=math.radians(-337.5)
ET_x_az45=math.radians(45)
ET_y_az45=math.radians(-315)
ET_x_az68=math.radians(22.5)
ET_y_az68=math.radians(-292.5)
ET_x_az90=math.radians(0)
ET_y_az90=math.radians(270)
ET_MR_elev=0
ET_D_long=math.radians(51.1)
ET_D_lat=math.radians(13.8)
ET_D_elev=200
#Detectors
ET_S_10=LAL_Detector(ET_S_long,ET_S_lat,ET_S_elev,ET_x_az,ET_y_az,RT_ET_10_psd)
ET_S_15=LAL_Detector(ET_S_long,ET_S_lat,ET_S_elev,ET_x_az,ET_y_az,RT_ET_15_psd)
ET_S_20=LAL_Detector(ET_S_long,ET_S_lat,ET_S_elev,ET_x_az,ET_y_az,RT_ET_20_psd)

ET_MR_10_0=LAL_Detector(ET_MR_long,ET_MR_lat,ET_MR_elev,ET_x_az,ET_y_az,RT_ET_10_psd)
ET_MR_15_0=LAL_Detector(ET_MR_long,ET_MR_lat,ET_MR_elev,ET_x_az,ET_y_az,RT_ET_15_psd)
ET_MR_20_0=LAL_Detector(ET_MR_long,ET_MR_lat,ET_MR_elev,ET_x_az,ET_y_az,RT_ET_20_psd)

ET_MR_10_22=LAL_Detector(ET_MR_long,ET_MR_lat,ET_MR_elev,ET_x_az22,ET_y_az22,RT_ET_10_psd)
ET_MR_15_22=LAL_Detector(ET_MR_long,ET_MR_lat,ET_MR_elev,ET_x_az22,ET_y_az22,RT_ET_15_psd)
ET_MR_20_22=LAL_Detector(ET_MR_long,ET_MR_lat,ET_MR_elev,ET_x_az22,ET_y_az22,RT_ET_20_psd)

ET_MR_10_45=LAL_Detector(ET_MR_long,ET_MR_lat,ET_MR_elev,ET_x_az45,ET_y_az45,RT_ET_10_psd)
ET_MR_15_45=LAL_Detector(ET_MR_long,ET_MR_lat,ET_MR_elev,ET_x_az45,ET_y_az45,RT_ET_15_psd)
ET_MR_20_45=LAL_Detector(ET_MR_long,ET_MR_lat,ET_MR_elev,ET_x_az45,ET_y_az45,RT_ET_20_psd)

ET_MR_10_68=LAL_Detector(ET_MR_long,ET_MR_lat,ET_MR_elev,ET_x_az68,ET_y_az68,RT_ET_10_psd)
ET_MR_15_68=LAL_Detector(ET_MR_long,ET_MR_lat,ET_MR_elev,ET_x_az68,ET_y_az68,RT_ET_15_psd)
ET_MR_20_68=LAL_Detector(ET_MR_long,ET_MR_lat,ET_MR_elev,ET_x_az68,ET_y_az68,RT_ET_20_psd)

ET_MR_10_90=LAL_Detector(ET_MR_long,ET_MR_lat,ET_MR_elev,ET_x_az90,ET_y_az90,RT_ET_10_psd)
ET_MR_15_90=LAL_Detector(ET_MR_long,ET_MR_lat,ET_MR_elev,ET_x_az90,ET_y_az90,RT_ET_15_psd)
ET_MR_20_90=LAL_Detector(ET_MR_long,ET_MR_lat,ET_MR_elev,ET_x_az90,ET_y_az90,RT_ET_20_psd)

ET_D_10_0=LAL_Detector(ET_D_long,ET_D_lat,ET_D_elev,ET_x_az,ET_y_az,RT_ET_10_psd)
ET_D_15_0=LAL_Detector(ET_D_long,ET_D_lat,ET_D_elev,ET_x_az,ET_y_az,RT_ET_15_psd)
ET_D_20_0=LAL_Detector(ET_D_long,ET_D_lat,ET_D_elev,ET_x_az,ET_y_az,RT_ET_20_psd)

ET_D_10_22=LAL_Detector(ET_D_long,ET_D_lat,ET_D_elev,ET_x_az22,ET_y_az22,RT_ET_10_psd)
ET_D_15_22=LAL_Detector(ET_D_long,ET_D_lat,ET_D_elev,ET_x_az22,ET_y_az22,RT_ET_15_psd)
ET_D_20_22=LAL_Detector(ET_D_long,ET_D_lat,ET_D_elev,ET_x_az22,ET_y_az22,RT_ET_20_psd)

ET_D_10_45=LAL_Detector(ET_D_long,ET_D_lat,ET_D_elev,ET_x_az45,ET_y_az45,RT_ET_10_psd)
ET_D_15_45=LAL_Detector(ET_D_long,ET_D_lat,ET_D_elev,ET_x_az45,ET_y_az45,RT_ET_15_psd)
ET_D_20_45=LAL_Detector(ET_D_long,ET_D_lat,ET_D_elev,ET_x_az45,ET_y_az45,RT_ET_20_psd)

ET_D_10_68=LAL_Detector(ET_D_long,ET_D_lat,ET_D_elev,ET_x_az68,ET_y_az68,RT_ET_10_psd)
ET_D_15_68=LAL_Detector(ET_D_long,ET_D_lat,ET_D_elev,ET_x_az68,ET_y_az68,RT_ET_15_psd)
ET_D_20_68=LAL_Detector(ET_D_long,ET_D_lat,ET_D_elev,ET_x_az68,ET_y_az68,RT_ET_20_psd)

ET_D_10_90=LAL_Detector(ET_D_long,ET_D_lat,ET_D_elev,ET_x_az90,ET_y_az90,RT_ET_10_psd)
ET_D_15_90=LAL_Detector(ET_D_long,ET_D_lat,ET_D_elev,ET_x_az90,ET_y_az90,RT_ET_15_psd)
ET_D_20_90=LAL_Detector(ET_D_long,ET_D_lat,ET_D_elev,ET_x_az90,ET_y_az90,RT_ET_20_psd)
#Networks
nr=np.size(RT_frequencies)
#S-MR
ET_S_MR_10_0=DetectorNetwork([ET_S_10,ET_MR_10_0],2,RT_frequencies,nr)
ET_S_MR_15_0=DetectorNetwork([ET_S_15,ET_MR_15_0],2,RT_frequencies,nr)
ET_S_MR_20_0=DetectorNetwork([ET_S_20,ET_MR_20_0],2,RT_frequencies,nr)

ET_S_MR_10_22=DetectorNetwork([ET_S_10,ET_MR_10_22],2,RT_frequencies,nr)
ET_S_MR_15_22=DetectorNetwork([ET_S_15,ET_MR_15_22],2,RT_frequencies,nr)
ET_S_MR_20_22=DetectorNetwork([ET_S_20,ET_MR_20_22],2,RT_frequencies,nr)

ET_S_MR_10_45=DetectorNetwork([ET_S_10,ET_MR_10_45],2,RT_frequencies,nr)
ET_S_MR_15_45=DetectorNetwork([ET_S_15,ET_MR_15_45],2,RT_frequencies,nr)
ET_S_MR_20_45=DetectorNetwork([ET_S_20,ET_MR_20_45],2,RT_frequencies,nr)

ET_S_MR_10_68=DetectorNetwork([ET_S_10,ET_MR_10_68],2,RT_frequencies,nr)
ET_S_MR_15_68=DetectorNetwork([ET_S_15,ET_MR_15_68],2,RT_frequencies,nr)
ET_S_MR_20_68=DetectorNetwork([ET_S_20,ET_MR_20_68],2,RT_frequencies,nr)

ET_S_MR_10_90=DetectorNetwork([ET_S_10,ET_MR_10_90],2,RT_frequencies,nr)
ET_S_MR_15_90=DetectorNetwork([ET_S_15,ET_MR_15_90],2,RT_frequencies,nr)
ET_S_MR_20_90=DetectorNetwork([ET_S_20,ET_MR_20_90],2,RT_frequencies,nr)
#S-D
ET_S_D_10_0=DetectorNetwork([ET_S_10,ET_D_10_0],2,RT_frequencies,nr)
ET_S_D_15_0=DetectorNetwork([ET_S_15,ET_D_15_0],2,RT_frequencies,nr)
ET_S_D_20_0=DetectorNetwork([ET_S_20,ET_D_20_0],2,RT_frequencies,nr)

ET_S_D_10_22=DetectorNetwork([ET_S_10,ET_D_10_22],2,RT_frequencies,nr)
ET_S_D_15_22=DetectorNetwork([ET_S_15,ET_D_15_22],2,RT_frequencies,nr)
ET_S_D_20_22=DetectorNetwork([ET_S_20,ET_D_20_22],2,RT_frequencies,nr)

ET_S_D_10_45=DetectorNetwork([ET_S_10,ET_D_10_45],2,RT_frequencies,nr)
ET_S_D_15_45=DetectorNetwork([ET_S_15,ET_D_15_45],2,RT_frequencies,nr)
ET_S_D_20_45=DetectorNetwork([ET_S_20,ET_D_20_45],2,RT_frequencies,nr)

ET_S_D_10_68=DetectorNetwork([ET_S_10,ET_D_10_68],2,RT_frequencies,nr)
ET_S_D_15_68=DetectorNetwork([ET_S_15,ET_D_15_68],2,RT_frequencies,nr)
ET_S_D_20_68=DetectorNetwork([ET_S_20,ET_D_20_68],2,RT_frequencies,nr)

ET_S_D_10_90=DetectorNetwork([ET_S_10,ET_D_10_90],2,RT_frequencies,nr)
ET_S_D_15_90=DetectorNetwork([ET_S_15,ET_D_15_90],2,RT_frequencies,nr)
ET_S_D_20_90=DetectorNetwork([ET_S_20,ET_D_20_90],2,RT_frequencies,nr)
#MR-D
ET_MR_D_10_0=DetectorNetwork([ET_MR_10_0,ET_D_10_0],2,RT_frequencies,nr)
ET_MR_D_15_0=DetectorNetwork([ET_MR_15_0,ET_D_15_0],2,RT_frequencies,nr)
ET_MR_D_20_0=DetectorNetwork([ET_MR_20_0,ET_D_20_0],2,RT_frequencies,nr)

ET_MR_D_10_22=DetectorNetwork([ET_MR_10_0,ET_D_10_22],2,RT_frequencies,nr)
ET_MR_D_15_22=DetectorNetwork([ET_MR_15_0,ET_D_15_22],2,RT_frequencies,nr)
ET_MR_D_20_22=DetectorNetwork([ET_MR_20_0,ET_D_20_22],2,RT_frequencies,nr)

ET_MR_D_10_45=DetectorNetwork([ET_MR_10_0,ET_D_10_45],2,RT_frequencies,nr)
ET_MR_D_15_45=DetectorNetwork([ET_MR_15_0,ET_D_15_45],2,RT_frequencies,nr)
ET_MR_D_20_45=DetectorNetwork([ET_MR_20_0,ET_D_20_45],2,RT_frequencies,nr)

ET_MR_D_10_68=DetectorNetwork([ET_MR_10_0,ET_D_10_68],2,RT_frequencies,nr)
ET_MR_D_15_68=DetectorNetwork([ET_MR_15_0,ET_D_15_68],2,RT_frequencies,nr)
ET_MR_D_20_68=DetectorNetwork([ET_MR_20_0,ET_D_20_68],2,RT_frequencies,nr)

ET_MR_D_10_90=DetectorNetwork([ET_MR_10_0,ET_D_10_90],2,RT_frequencies,nr)
ET_MR_D_15_90=DetectorNetwork([ET_MR_15_0,ET_D_15_90],2,RT_frequencies,nr)
ET_MR_D_20_90=DetectorNetwork([ET_MR_20_0,ET_D_20_90],2,RT_frequencies,nr)

beta_max=8
betas=np.linspace(-beta_max, beta_max, 2*beta_max+1)
f_ref=100
#PI-curves
PI_S_MR_10_0,O_S_MR_10_0=ET_S_MR_10_0.Power_Integration(betas)
PI_S_MR_15_0,O_S_MR_15_0=ET_S_MR_15_0.Power_Integration(betas)
PI_S_MR_20_0,O_S_MR_20_0=ET_S_MR_20_0.Power_Integration(betas)

PI_S_MR_10_22,O_S_MR_10_22=ET_S_MR_10_22.Power_Integration(betas)
PI_S_MR_15_22,O_S_MR_15_22=ET_S_MR_15_22.Power_Integration(betas)
PI_S_MR_20_22,O_S_MR_20_22=ET_S_MR_20_22.Power_Integration(betas)

PI_S_MR_10_45,O_S_MR_10_45=ET_S_MR_10_45.Power_Integration(betas)
PI_S_MR_15_45,O_S_MR_15_45=ET_S_MR_15_45.Power_Integration(betas)
PI_S_MR_20_45,O_S_MR_20_45=ET_S_MR_20_45.Power_Integration(betas)

PI_S_MR_10_68,O_S_MR_10_68=ET_S_MR_10_68.Power_Integration(betas)
PI_S_MR_15_68,O_S_MR_15_68=ET_S_MR_15_68.Power_Integration(betas)
PI_S_MR_20_68,O_S_MR_20_68=ET_S_MR_20_68.Power_Integration(betas)

PI_S_MR_10_90,O_S_MR_10_90=ET_S_MR_10_90.Power_Integration(betas)
PI_S_MR_15_90,O_S_MR_15_90=ET_S_MR_15_90.Power_Integration(betas)
PI_S_MR_20_90,O_S_MR_20_90=ET_S_MR_20_90.Power_Integration(betas)

PI_S_D_10_0,O_S_D_10_0=ET_S_D_10_0.Power_Integration(betas)
PI_S_D_15_0,O_S_D_15_0=ET_S_D_15_0.Power_Integration(betas)
PI_S_D_20_0,O_S_D_20_0=ET_S_D_20_0.Power_Integration(betas)

PI_S_D_10_22,O_S_D_10_22=ET_S_D_10_22.Power_Integration(betas)
PI_S_D_15_22,O_S_D_15_22=ET_S_D_15_22.Power_Integration(betas)
PI_S_D_20_22,O_S_D_20_22=ET_S_D_20_22.Power_Integration(betas)

PI_S_D_10_45,O_S_D_10_45=ET_S_D_10_45.Power_Integration(betas)
PI_S_D_15_45,O_S_D_15_45=ET_S_D_15_45.Power_Integration(betas)
PI_S_D_20_45,O_S_D_20_45=ET_S_D_20_45.Power_Integration(betas)

PI_S_D_10_68,O_S_D_10_68=ET_S_D_10_68.Power_Integration(betas)
PI_S_D_15_68,O_S_D_15_68=ET_S_D_15_68.Power_Integration(betas)
PI_S_D_20_68,O_S_D_20_68=ET_S_D_20_68.Power_Integration(betas)

PI_S_D_10_90,O_S_D_10_90=ET_S_D_10_90.Power_Integration(betas)
PI_S_D_15_90,O_S_D_15_90=ET_S_D_15_90.Power_Integration(betas)
PI_S_D_20_90,O_S_D_20_90=ET_S_D_20_90.Power_Integration(betas)

PI_MR_D_10_0,O_MR_D_10_0=ET_MR_D_10_0.Power_Integration(betas)
PI_MR_D_15_0,O_MR_D_15_0=ET_MR_D_15_0.Power_Integration(betas)
PI_MR_D_20_0,O_MR_D_20_0=ET_MR_D_20_0.Power_Integration(betas)

PI_MR_D_10_22,O_MR_D_10_22=ET_MR_D_10_22.Power_Integration(betas)
PI_MR_D_15_22,O_MR_D_15_22=ET_MR_D_15_22.Power_Integration(betas)
PI_MR_D_20_22,O_MR_D_20_22=ET_MR_D_20_22.Power_Integration(betas)

PI_MR_D_10_45,O_MR_D_10_45=ET_MR_D_10_45.Power_Integration(betas)
PI_MR_D_15_45,O_MR_D_15_45=ET_MR_D_15_45.Power_Integration(betas)
PI_MR_D_20_45,O_MR_D_20_45=ET_MR_D_20_45.Power_Integration(betas)

PI_MR_D_10_68,O_MR_D_10_68=ET_MR_D_10_68.Power_Integration(betas)
PI_MR_D_15_68,O_MR_D_15_68=ET_MR_D_15_68.Power_Integration(betas)
PI_MR_D_20_68,O_MR_D_20_68=ET_MR_D_20_68.Power_Integration(betas)

PI_MR_D_10_90,O_MR_D_10_90=ET_MR_D_10_90.Power_Integration(betas)
PI_MR_D_15_90,O_MR_D_15_90=ET_MR_D_15_90.Power_Integration(betas)
PI_MR_D_20_90,O_MR_D_20_90=ET_MR_D_20_90.Power_Integration(betas)
#plots -15_1
# fig,axs = plt.subplots(3,3)
# axs[0,0].set_title('Sardinia-MeuseRhein-10km')
# axs[0,0].set_xscale('log')
# axs[0,0].set_yscale('log')
# axs[0,0].set_xlim(starting_frequency, ending_frequency)
# axs[0,0].set_ylim(1e-15,1e0)
# axs[0,0].plot(RT_frequencies,PI_S_MR_10_0,color='black',label='theta=0')
# axs[0,0].plot(RT_frequencies,PI_S_MR_10_22,color='blue',label='theta=22.5')
# axs[0,0].plot(RT_frequencies,PI_S_MR_10_45,color='red',label='theta=45')
# axs[0,0].plot(RT_frequencies,PI_S_MR_10_68,color='green',label='theta=67.5')
# axs[0,0].plot(RT_frequencies,PI_S_MR_10_90,color='purple',label='theta=90')
# axs[0,0].legend()
#
# axs[0,1].set_title('Sardinia-Dresden-10km')
# axs[0,1].set_xscale('log')
# axs[0,1].set_yscale('log')
# axs[0,1].set_xlim(starting_frequency, ending_frequency)
# axs[0,1].set_ylim(1e-15,1e0)
# axs[0,1].plot(RT_frequencies,PI_S_D_10_0,color='black',label='theta=0')
# axs[0,1].plot(RT_frequencies,PI_S_D_10_22,color='blue',label='theta=22.5')
# axs[0,1].plot(RT_frequencies,PI_S_D_10_45,color='red',label='theta=45')
# axs[0,1].plot(RT_frequencies,PI_S_D_10_68,color='green',label='theta=67.5')
# axs[0,1].plot(RT_frequencies,PI_S_D_10_90,color='purple',label='theta=90')
# axs[0,1].legend()
#
# axs[0,2].set_title('MeuseRhein-Dresden-10km')
# axs[0,2].set_xscale('log')
# axs[0,2].set_yscale('log')
# axs[0,2].set_xlim(starting_frequency, ending_frequency)
# axs[0,2].set_ylim(1e-15,1e0)
# axs[0,2].plot(RT_frequencies,PI_MR_D_10_0,color='black',label='theta=0')
# axs[0,2].plot(RT_frequencies,PI_MR_D_10_22,color='blue',label='theta=22.5')
# axs[0,2].plot(RT_frequencies,PI_MR_D_10_45,color='red',label='theta=45')
# axs[0,2].plot(RT_frequencies,PI_MR_D_10_68,color='green',label='theta=67.5')
# axs[0,2].plot(RT_frequencies,PI_MR_D_10_90,color='purple',label='theta=90')
# axs[0,2].legend()
#
# axs[1,0].set_title('Sardinia-MeuseRhein-15km')
# axs[1,0].set_xscale('log')
# axs[1,0].set_yscale('log')
# axs[1,0].set_xlim(starting_frequency, ending_frequency)
# axs[1,0].set_ylim(1e-15,1e0)
# axs[1,0].plot(RT_frequencies,PI_S_MR_15_0,color='black',label='theta=0')
# axs[1,0].plot(RT_frequencies,PI_S_MR_15_22,color='blue',label='theta=22.5')
# axs[1,0].plot(RT_frequencies,PI_S_MR_15_45,color='red',label='theta=45')
# axs[1,0].plot(RT_frequencies,PI_S_MR_15_68,color='green',label='theta=67.5')
# axs[1,0].plot(RT_frequencies,PI_S_MR_15_90,color='purple',label='theta=90')
# axs[1,0].legend()
#
# axs[1,1].set_title('Sardinia-Dresden-15km')
# axs[1,1].set_xscale('log')
# axs[1,1].set_yscale('log')
# axs[1,1].set_xlim(starting_frequency, ending_frequency)
# axs[1,1].set_ylim(1e-15,1e0)
# axs[1,1].plot(RT_frequencies,PI_S_D_15_0,color='black',label='theta=0')
# axs[1,1].plot(RT_frequencies,PI_S_D_15_22,color='blue',label='theta=22.5')
# axs[1,1].plot(RT_frequencies,PI_S_D_15_45,color='red',label='theta=45')
# axs[1,1].plot(RT_frequencies,PI_S_D_15_68,color='green',label='theta=67.5')
# axs[1,1].plot(RT_frequencies,PI_S_D_15_90,color='purple',label='theta=90')
# axs[1,1].legend()
#
# axs[1,2].set_title('MeuseRhein-Dresden-15km')
# axs[1,2].set_xscale('log')
# axs[1,2].set_yscale('log')
# axs[1,2].set_xlim(starting_frequency, ending_frequency)
# axs[1,2].set_ylim(1e-15,1e0)
# axs[1,2].plot(RT_frequencies,PI_MR_D_15_0,color='black',label='theta=0')
# axs[1,2].plot(RT_frequencies,PI_MR_D_15_22,color='blue',label='theta=22.5')
# axs[1,2].plot(RT_frequencies,PI_MR_D_15_45,color='red',label='theta=45')
# axs[1,2].plot(RT_frequencies,PI_MR_D_15_68,color='green',label='theta=67.5')
# axs[1,2].plot(RT_frequencies,PI_MR_D_15_90,color='purple',label='theta=90')
# axs[1,2].legend()
#
# axs[2,0].set_title('Sardinia-MeuseRhein-20km')
# axs[2,0].set_xscale('log')
# axs[2,0].set_yscale('log')
# axs[2,0].set_xlim(starting_frequency, ending_frequency)
# axs[2,0].set_ylim(1e-15,1e0)
# axs[2,0].plot(RT_frequencies,PI_S_MR_20_0,color='black',label='theta=0')
# axs[2,0].plot(RT_frequencies,PI_S_MR_20_22,color='blue',label='theta=22.5')
# axs[2,0].plot(RT_frequencies,PI_S_MR_20_45,color='red',label='theta=45')
# axs[2,0].plot(RT_frequencies,PI_S_MR_20_68,color='green',label='theta=67.5')
# axs[2,0].plot(RT_frequencies,PI_S_MR_20_90,color='purple',label='theta=90')
# axs[2,0].legend()
#
# axs[2,1].set_title('Sardinia-Dresden-20km')
# axs[2,1].set_xscale('log')
# axs[2,1].set_yscale('log')
# axs[2,1].set_xlim(starting_frequency, ending_frequency)
# axs[2,1].set_ylim(1e-15,1e0)
# axs[2,1].plot(RT_frequencies,PI_S_D_20_0,color='black',label='theta=0')
# axs[2,1].plot(RT_frequencies,PI_S_D_20_22,color='blue',label='theta=22.5')
# axs[2,1].plot(RT_frequencies,PI_S_D_20_45,color='red',label='theta=45')
# axs[2,1].plot(RT_frequencies,PI_S_D_20_68,color='green',label='theta=67.5')
# axs[2,1].plot(RT_frequencies,PI_S_D_20_90,color='purple',label='theta=90')
# axs[2,1].legend()
#
# axs[2,2].set_title('MeuseRhein-Dresden-20km')
# axs[2,2].set_xscale('log')
# axs[2,2].set_yscale('log')
# axs[2,2].set_xlim(starting_frequency, ending_frequency)
# axs[2,2].set_ylim(1e-15,1e0)
# axs[2,2].plot(RT_frequencies,PI_MR_D_20_0,color='black',label='theta=0')
# axs[2,2].plot(RT_frequencies,PI_MR_D_20_22,color='blue',label='theta=22.5')
# axs[2,2].plot(RT_frequencies,PI_MR_D_20_45,color='red',label='theta=45')
# axs[2,2].plot(RT_frequencies,PI_MR_D_20_68,color='green',label='theta=67.5')
# axs[2,2].plot(RT_frequencies,PI_MR_D_20_90,color='purple',label='theta=90')
# axs[2,2].legend()
#plots-15_10
# fig,axs = plt.subplots(3,3)
# axs[0,0].set_title('Sardinia-MeuseRhein-10km')
# axs[0,0].set_xscale('log')
# axs[0,0].set_yscale('log')
# axs[0,0].set_xlim(starting_frequency, ending_frequency)
# axs[0,0].set_ylim(1e-15,1e10)
# axs[0,0].plot(RT_frequencies,PI_S_MR_10_0,color='black',label='theta=0')
# axs[0,0].plot(RT_frequencies,PI_S_MR_10_22,color='blue',label='theta=22.5')
# axs[0,0].plot(RT_frequencies,PI_S_MR_10_45,color='red',label='theta=45')
# axs[0,0].plot(RT_frequencies,PI_S_MR_10_68,color='green',label='theta=67.5')
# axs[0,0].plot(RT_frequencies,PI_S_MR_10_90,color='purple',label='theta=90')
# axs[0,0].legend()
#
# axs[0,1].set_title('Sardinia-Dresden-10km')
# axs[0,1].set_xscale('log')
# axs[0,1].set_yscale('log')
# axs[0,1].set_xlim(starting_frequency, ending_frequency)
# axs[0,1].set_ylim(1e-15,1e10)
# axs[0,1].plot(RT_frequencies,PI_S_D_10_0,color='black',label='theta=0')
# axs[0,1].plot(RT_frequencies,PI_S_D_10_22,color='blue',label='theta=22.5')
# axs[0,1].plot(RT_frequencies,PI_S_D_10_45,color='red',label='theta=45')
# axs[0,1].plot(RT_frequencies,PI_S_D_10_68,color='green',label='theta=67.5')
# axs[0,1].plot(RT_frequencies,PI_S_D_10_90,color='purple',label='theta=90')
# axs[0,1].legend()
#
# axs[0,2].set_title('MeuseRhein-Dresden-10km')
# axs[0,2].set_xscale('log')
# axs[0,2].set_yscale('log')
# axs[0,2].set_xlim(starting_frequency, ending_frequency)
# axs[0,2].set_ylim(1e-15,1e10)
# axs[0,2].plot(RT_frequencies,PI_MR_D_10_0,color='black',label='theta=0')
# axs[0,2].plot(RT_frequencies,PI_MR_D_10_22,color='blue',label='theta=22.5')
# axs[0,2].plot(RT_frequencies,PI_MR_D_10_45,color='red',label='theta=45')
# axs[0,2].plot(RT_frequencies,PI_MR_D_10_68,color='green',label='theta=67.5')
# axs[0,2].plot(RT_frequencies,PI_MR_D_10_90,color='purple',label='theta=90')
# axs[0,2].legend()
#
# axs[1,0].set_title('Sardinia-MeuseRhein-15km')
# axs[1,0].set_xscale('log')
# axs[1,0].set_yscale('log')
# axs[1,0].set_xlim(starting_frequency, ending_frequency)
# axs[1,0].set_ylim(1e-15,1e10)
# axs[1,0].plot(RT_frequencies,PI_S_MR_15_0,color='black',label='theta=0')
# axs[1,0].plot(RT_frequencies,PI_S_MR_15_22,color='blue',label='theta=22.5')
# axs[1,0].plot(RT_frequencies,PI_S_MR_15_45,color='red',label='theta=45')
# axs[1,0].plot(RT_frequencies,PI_S_MR_15_68,color='green',label='theta=67.5')
# axs[1,0].plot(RT_frequencies,PI_S_MR_15_90,color='purple',label='theta=90')
# axs[1,0].legend()
#
# axs[1,1].set_title('Sardinia-Dresden-15km')
# axs[1,1].set_xscale('log')
# axs[1,1].set_yscale('log')
# axs[1,1].set_xlim(starting_frequency, ending_frequency)
# axs[1,1].set_ylim(1e-15,1e10)
# axs[1,1].plot(RT_frequencies,PI_S_D_15_0,color='black',label='theta=0')
# axs[1,1].plot(RT_frequencies,PI_S_D_15_22,color='blue',label='theta=22.5')
# axs[1,1].plot(RT_frequencies,PI_S_D_15_45,color='red',label='theta=45')
# axs[1,1].plot(RT_frequencies,PI_S_D_15_68,color='green',label='theta=67.5')
# axs[1,1].plot(RT_frequencies,PI_S_D_15_90,color='purple',label='theta=90')
# axs[1,1].legend()
#
# axs[1,2].set_title('MeuseRhein-Dresden-15km')
# axs[1,2].set_xscale('log')
# axs[1,2].set_yscale('log')
# axs[1,2].set_xlim(starting_frequency, ending_frequency)
# axs[1,2].set_ylim(1e-15,1e10)
# axs[1,2].plot(RT_frequencies,PI_MR_D_15_0,color='black',label='theta=0')
# axs[1,2].plot(RT_frequencies,PI_MR_D_15_22,color='blue',label='theta=22.5')
# axs[1,2].plot(RT_frequencies,PI_MR_D_15_45,color='red',label='theta=45')
# axs[1,2].plot(RT_frequencies,PI_MR_D_15_68,color='green',label='theta=67.5')
# axs[1,2].plot(RT_frequencies,PI_MR_D_15_90,color='purple',label='theta=90')
# axs[1,2].legend()
#
# axs[2,0].set_title('Sardinia-MeuseRhein-20km')
# axs[2,0].set_xscale('log')
# axs[2,0].set_yscale('log')
# axs[2,0].set_xlim(starting_frequency, ending_frequency)
# axs[2,0].set_ylim(1e-15,1e10)
# axs[2,0].plot(RT_frequencies,PI_S_MR_20_0,color='black',label='theta=0')
# axs[2,0].plot(RT_frequencies,PI_S_MR_20_22,color='blue',label='theta=22.5')
# axs[2,0].plot(RT_frequencies,PI_S_MR_20_45,color='red',label='theta=45')
# axs[2,0].plot(RT_frequencies,PI_S_MR_20_68,color='green',label='theta=67.5')
# axs[2,0].plot(RT_frequencies,PI_S_MR_20_90,color='purple',label='theta=90')
# axs[2,0].legend()
#
# axs[2,1].set_title('Sardinia-Dresden-20km')
# axs[2,1].set_xscale('log')
# axs[2,1].set_yscale('log')
# axs[2,1].set_xlim(starting_frequency, ending_frequency)
# axs[2,1].set_ylim(1e-15,1e10)
# axs[2,1].plot(RT_frequencies,PI_S_D_20_0,color='black',label='theta=0')
# axs[2,1].plot(RT_frequencies,PI_S_D_20_22,color='blue',label='theta=22.5')
# axs[2,1].plot(RT_frequencies,PI_S_D_20_45,color='red',label='theta=45')
# axs[2,1].plot(RT_frequencies,PI_S_D_20_68,color='green',label='theta=67.5')
# axs[2,1].plot(RT_frequencies,PI_S_D_20_90,color='purple',label='theta=90')
# axs[2,1].legend()
#
# axs[2,2].set_title('MeuseRhein-Dresden-20km')
# axs[2,2].set_xscale('log')
# axs[2,2].set_yscale('log')
# axs[2,2].set_xlim(starting_frequency, ending_frequency)
# axs[2,2].set_ylim(1e-15,1e10)
# axs[2,2].plot(RT_frequencies,PI_MR_D_20_0,color='black',label='theta=0')
# axs[2,2].plot(RT_frequencies,PI_MR_D_20_22,color='blue',label='theta=22.5')
# axs[2,2].plot(RT_frequencies,PI_MR_D_20_45,color='red',label='theta=45')
# axs[2,2].plot(RT_frequencies,PI_MR_D_20_68,color='green',label='theta=67.5')
# axs[2,2].plot(RT_frequencies,PI_MR_D_20_90,color='purple',label='theta=90')
# axs[2,2].legend()
#plots-12_-6
fig,axs = plt.subplots(3,3)
axs[0,0].set_title('Sardinia-MeuseRhein-10km')
axs[0,0].set_xscale('log')
axs[0,0].set_yscale('log')
axs[0,0].set_xlim(starting_frequency, ending_frequency)
axs[0,0].set_ylim(1e-12,1e-6)
axs[0,0].plot(RT_frequencies,PI_S_MR_10_0,color='black',label='theta=0')
axs[0,0].plot(RT_frequencies,PI_S_MR_10_22,color='blue',label='theta=22.5')
axs[0,0].plot(RT_frequencies,PI_S_MR_10_45,color='red',label='theta=45')
axs[0,0].plot(RT_frequencies,PI_S_MR_10_68,color='green',label='theta=67.5')
axs[0,0].plot(RT_frequencies,PI_S_MR_10_90,color='purple',label='theta=90')
axs[0,0].legend()

axs[0,1].set_title('Sardinia-Dresden-10km')
axs[0,1].set_xscale('log')
axs[0,1].set_yscale('log')
axs[0,1].set_xlim(starting_frequency, ending_frequency)
axs[0,1].set_ylim(1e-12,1e-6)
axs[0,1].plot(RT_frequencies,PI_S_D_10_0,color='black',label='theta=0')
axs[0,1].plot(RT_frequencies,PI_S_D_10_22,color='blue',label='theta=22.5')
axs[0,1].plot(RT_frequencies,PI_S_D_10_45,color='red',label='theta=45')
axs[0,1].plot(RT_frequencies,PI_S_D_10_68,color='green',label='theta=67.5')
axs[0,1].plot(RT_frequencies,PI_S_D_10_90,color='purple',label='theta=90')
axs[0,1].legend()

axs[0,2].set_title('MeuseRhein-Dresden-10km')
axs[0,2].set_xscale('log')
axs[0,2].set_yscale('log')
axs[0,2].set_xlim(starting_frequency, ending_frequency)
axs[0,2].set_ylim(1e-12,1e-6)
axs[0,2].plot(RT_frequencies,PI_MR_D_10_0,color='black',label='theta=0')
axs[0,2].plot(RT_frequencies,PI_MR_D_10_22,color='blue',label='theta=22.5')
axs[0,2].plot(RT_frequencies,PI_MR_D_10_45,color='red',label='theta=45')
axs[0,2].plot(RT_frequencies,PI_MR_D_10_68,color='green',label='theta=67.5')
axs[0,2].plot(RT_frequencies,PI_MR_D_10_90,color='purple',label='theta=90')
axs[0,2].legend()

axs[1,0].set_title('Sardinia-MeuseRhein-15km')
axs[1,0].set_xscale('log')
axs[1,0].set_yscale('log')
axs[1,0].set_xlim(starting_frequency, ending_frequency)
axs[1,0].set_ylim(1e-12,1e-6)
axs[1,0].plot(RT_frequencies,PI_S_MR_15_0,color='black',label='theta=0')
axs[1,0].plot(RT_frequencies,PI_S_MR_15_22,color='blue',label='theta=22.5')
axs[1,0].plot(RT_frequencies,PI_S_MR_15_45,color='red',label='theta=45')
axs[1,0].plot(RT_frequencies,PI_S_MR_15_68,color='green',label='theta=67.5')
axs[1,0].plot(RT_frequencies,PI_S_MR_15_90,color='purple',label='theta=90')
axs[1,0].legend()

axs[1,1].set_title('Sardinia-Dresden-15km')
axs[1,1].set_xscale('log')
axs[1,1].set_yscale('log')
axs[1,1].set_xlim(starting_frequency, ending_frequency)
axs[1,1].set_ylim(1e-12,1e-6)
axs[1,1].plot(RT_frequencies,PI_S_D_15_0,color='black',label='theta=0')
axs[1,1].plot(RT_frequencies,PI_S_D_15_22,color='blue',label='theta=22.5')
axs[1,1].plot(RT_frequencies,PI_S_D_15_45,color='red',label='theta=45')
axs[1,1].plot(RT_frequencies,PI_S_D_15_68,color='green',label='theta=67.5')
axs[1,1].plot(RT_frequencies,PI_S_D_15_90,color='purple',label='theta=90')
axs[1,1].legend()

axs[1,2].set_title('MeuseRhein-Dresden-15km')
axs[1,2].set_xscale('log')
axs[1,2].set_yscale('log')
axs[1,2].set_xlim(starting_frequency, ending_frequency)
axs[1,2].set_ylim(1e-12,1e-6)
axs[1,2].plot(RT_frequencies,PI_MR_D_15_0,color='black',label='theta=0')
axs[1,2].plot(RT_frequencies,PI_MR_D_15_22,color='blue',label='theta=22.5')
axs[1,2].plot(RT_frequencies,PI_MR_D_15_45,color='red',label='theta=45')
axs[1,2].plot(RT_frequencies,PI_MR_D_15_68,color='green',label='theta=67.5')
axs[1,2].plot(RT_frequencies,PI_MR_D_15_90,color='purple',label='theta=90')
axs[1,2].legend()

axs[2,0].set_title('Sardinia-MeuseRhein-20km')
axs[2,0].set_xscale('log')
axs[2,0].set_yscale('log')
axs[2,0].set_xlim(starting_frequency, ending_frequency)
axs[2,0].set_ylim(1e-12,1e-6)
axs[2,0].plot(RT_frequencies,PI_S_MR_20_0,color='black',label='theta=0')
axs[2,0].plot(RT_frequencies,PI_S_MR_20_22,color='blue',label='theta=22.5')
axs[2,0].plot(RT_frequencies,PI_S_MR_20_45,color='red',label='theta=45')
axs[2,0].plot(RT_frequencies,PI_S_MR_20_68,color='green',label='theta=67.5')
axs[2,0].plot(RT_frequencies,PI_S_MR_20_90,color='purple',label='theta=90')
axs[2,0].legend()

axs[2,1].set_title('Sardinia-Dresden-20km')
axs[2,1].set_xscale('log')
axs[2,1].set_yscale('log')
axs[2,1].set_xlim(starting_frequency, ending_frequency)
axs[2,1].set_ylim(1e-12,1e-6)
axs[2,1].plot(RT_frequencies,PI_S_D_20_0,color='black',label='theta=0')
axs[2,1].plot(RT_frequencies,PI_S_D_20_22,color='blue',label='theta=22.5')
axs[2,1].plot(RT_frequencies,PI_S_D_20_45,color='red',label='theta=45')
axs[2,1].plot(RT_frequencies,PI_S_D_20_68,color='green',label='theta=67.5')
axs[2,1].plot(RT_frequencies,PI_S_D_20_90,color='purple',label='theta=90')
axs[2,1].legend()

axs[2,2].set_title('MeuseRhein-Dresden-20km')
axs[2,2].set_xscale('log')
axs[2,2].set_yscale('log')
axs[2,2].set_xlim(starting_frequency, ending_frequency)
axs[2,2].set_ylim(1e-12,1e-6)
axs[2,2].plot(RT_frequencies,PI_MR_D_20_0,color='black',label='theta=0')
axs[2,2].plot(RT_frequencies,PI_MR_D_20_22,color='blue',label='theta=22.5')
axs[2,2].plot(RT_frequencies,PI_MR_D_20_45,color='red',label='theta=45')
axs[2,2].plot(RT_frequencies,PI_MR_D_20_68,color='green',label='theta=67.5')
axs[2,2].plot(RT_frequencies,PI_MR_D_20_90,color='purple',label='theta=90')
axs[2,2].legend()

end = time()
print(end-start, "s")
plt.show()

ET1_long=math.radians(43.63)
ET1_lat=math.radians(10.50)
ET1_xaz=math.radians(19.43)
ET1_yaz=math.radians(319.43)
ET1_elev=52

ET2_long=math.radians(43.72)
ET2_lat=math.radians(10.55)
ET2_xaz=math.radians(259.43)
ET2_yaz=math.radians(199.43)
ET2_elev=60

ET3_long=math.radians(43.70)
ET3_lat=math.radians(10.42)
ET3_xaz=math.radians(139.43)
ET3_yaz=math.radians(79.43)
ET3_elev=60

ET1=LAL_Detector(ET1_long,ET1_lat,ET1_xaz,ET1_yaz,ET1_elev,RT_ET_10_psd)
ET2=LAL_Detector(ET2_long,ET2_lat,ET2_xaz,ET2_yaz,ET2_elev,RT_ET_10_psd)
ET3=LAL_Detector(ET3_long,ET3_lat,ET3_xaz,ET3_yaz,ET3_elev,RT_ET_10_psd)

ET123=[ET1,ET2,ET3]
ET_Triangle=DetectorNetwork(ET123,3,RT_frequencies,nr)
PI_ET_Triangle,O_ET_Triangle=ET_Triangle.Power_Integration(betas)
plt.loglog(RT_frequencies, PI_ET_Triangle, color='black')
plt.title('ET-Triangle')
plt.show()



# LIGO_Network=DetectorNetwork(LIGO_detectors,len(LIGO_detectors),RT_frequencies,np.size(RT_frequencies))
# Omega_0=LIGO_Network.get_Omega_eff(0)
# Omega_1=LIGO_Network.get_Omega_eff(1)
# PI_curve,Omega_Bs=LIGO_Network.Power_Integration(betas)
# plt.loglog(RT_frequencies, Omega_RT, color='black')
# plt.loglog(RT_frequencies,Omega_0,color='red')
# plt.loglog(RT_frequencies,Omega_1,color='green')
# for beta, Omega_B in zip(betas, Omega_Bs):
#     plt.loglog(RT_frequencies, Omega_B*np.power(RT_frequencies/f_ref, beta), color = "gray")
# plt.loglog(RT_frequencies,PI_curve,color='blue')
# plt.xlim(10,1500)
# plt.ylim(1e-10,1e0)
# plt.show()



# lat=0.81079526383
# long=-2.08405676917
# elev=142.554
# x_az=5.65487724844
# y_az=4.08408092164
# LigoLAL=LAL_Detector(long,lat,elev,x_az,y_az)
# print(LigoLAL.x)
# print(LigoLAL.u)
# print(LigoLAL.v)

# LIGO=DetectorPair(LIGO1,LIGO2)
# frequencies=np.linspace(1, 1000, 1000)
# gammas=np.zeros(len(frequencies))
# start = time()
# for i in range(len(frequencies)):
#     LIGO.calculate_alpha(frequencies[i])
#     gammas[i]=LIGO.overlap_function()
# end = time()
# print(end-start, "s")

# LIGO = DetectorPairVec(LIGO1, LIGO2)
# frequencies = np.linspace(1, 1000, 1000)
#
# start = time()
# gammas = LIGO.overlap_function_vec(frequencies)
# end = time()
# print(end - start, "s")
#
# fig, (ax_lin, ax_log) = plt.subplots(1, 2, figsize=(11, 4.5))
#
# ax_lin.plot(frequencies, gammas)
# ax_lin.set_xlabel("f (Hz)")
# ax_lin.set_ylabel(r"$\gamma(f)$")
# ax_lin.set_title("Linear frequency scale")
# ax_lin.grid(alpha=0.3)
#
# ax_log.plot(frequencies, gammas)
# ax_log.set_xscale("log")
# ax_log.set_xlabel("f (Hz)")
# ax_log.set_ylabel(r"$\gamma(f)$")
# ax_log.set_title("Logarithmic frequency scale")
# ax_log.grid(alpha=0.3, which="both")
#
# fig.suptitle("Overlap function")
# fig.tight_layout()
# plt.show()
