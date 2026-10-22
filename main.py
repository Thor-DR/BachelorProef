import numpy as np
import matplotlib.pyplot as plt
from classes import Detector, LAL_Detector
from classes import DetectorPair
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

H0=2.2*10**(-18)
start = time()
starting_frequencie=10
ending_frequencie=1500
df=0.25
df_inv=int(1/df)
RT_frequencies=np.linspace(starting_frequencie, ending_frequencie, (ending_frequencie-starting_frequencie)*df_inv)
aligo_design ="C:/Users/Thor De Rouck/Documents/0FysicaUAntwerpen/BachelorProef/aligo.txt"
f_aligo_design, asd_aligo_design = np.transpose(np.loadtxt(aligo_design))
RT_asd=interpolate_asd(RT_frequencies,f_aligo_design,asd_aligo_design)
Omega_RT= (2 * np.pi ** 2 * RT_frequencies ** 3 * RT_asd ** 2) / (3 * H0 ** 2)
end = time()
print(end-start, "s")

u1=np.array([-0.22389266154, 0.79983062746, 0.55690487831])
v1=np.array([-0.91397805644, 0.02609403989, -0.40492342125])
x1=np.array([-2161414.92636, -3834695.17889, 4600350.22664])
LIGO1=Detector(x1,u1,v1,RT_asd)

u2=np.array([-0.95457412153, -0.14158077340, -0.26218911324])
v2=np.array([0.29774156894, -0.48791033647, -0.82054461286])
x2=np.array([-74276.04238, -5496283.71971, 3224257.01744])
LIGO2=Detector(x2,u2,v2,RT_asd)
LIGO_detectors=[LIGO1,LIGO2]
beta_max=8
betas=np.linspace(-beta_max, beta_max, 2*beta_max+1)
f_ref=100

LIGO_Network=DetectorNetwork(LIGO_detectors,len(LIGO_detectors),RT_frequencies,np.size(RT_frequencies))
Omega_0=LIGO_Network.get_Omega_eff(0)
Omega_1=LIGO_Network.get_Omega_eff(1)
PI_curve,Omega_Bs=LIGO_Network.Power_Integration(betas)
plt.loglog(RT_frequencies, Omega_RT, color='black')
plt.loglog(RT_frequencies,Omega_0,color='red')
plt.loglog(RT_frequencies,Omega_1,color='green')
for beta, Omega_B in zip(betas, Omega_Bs):
    plt.loglog(RT_frequencies, Omega_B*np.power(RT_frequencies/f_ref, beta), color = "gray")
plt.loglog(RT_frequencies,PI_curve,color='blue')
plt.xlim(10,1500)
plt.ylim(1e-10,1e0)
plt.show()



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
