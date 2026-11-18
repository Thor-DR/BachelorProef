import math

import numpy as np
import matplotlib.pyplot as plt
from classes import Detector
from classes import DetectorPairVec
from time import time
from main import interpolate_asd
#hoimport healpy as hp

theta = np.radians(np.linspace(-90, 90, 500))[None, :]   # (1, 181)
phi   = np.radians(np.linspace(-180, 180, 500))[:, None]   # (360, 1)

ones = np.ones((500, 500))

class Vector3:
    def __init__(self, x, y, z):
        self.x, self.y, self.z = x, y, z   # each (360, 181) array

n = Vector3(np.sin(theta)*np.cos(phi),
            np.sin(theta)*np.sin(phi),
            np.cos(theta)*ones)

l = Vector3(np.cos(theta)*np.cos(phi),
            np.cos(theta)*np.sin(phi),
            -np.sin(theta)*ones)

m = Vector3(-np.sin(phi)*ones,
            np.cos(phi)*ones,
            np.zeros((500, 500)))

# stack components: shape (360, 181, 3)
l_stack = np.stack([l.x, l.y, l.z], axis=-1)
M = np.stack([m.x, m.y, m.z], axis=-1)

# e+_ab = l_a l_b - m_a m_b  -> shape (360, 181, 3, 3)
e_plus = np.einsum('...a,...b->...ab', l_stack, l_stack) - np.einsum('...a,...b->...ab', M, M)
# e_cross_ab = l_a m_b + m_a l_b  -> shape (360, 181, 3, 3)
e_cross = np.einsum('...a,...b->...ab', l_stack, M) + np.einsum('...a,...b->...ab', M, l_stack)

H0=2.2*10**(-18)
RT_frequencies=np.linspace(1, 1000)
aligo_design ="C:/Users/Thor De Rouck/Documents/0FysicaUAntwerpen/BachelorProef/aligo.txt"
f_aligo_design, asd_aligo_design = np.transpose(np.loadtxt(aligo_design))
RT_asd=interpolate_asd(RT_frequencies,f_aligo_design,asd_aligo_design)
RT_psd=RT_asd**2
Omega_RT= (2 * np.pi ** 2 * RT_frequencies ** 3 * RT_psd) / (3 * H0 ** 2)

u1=np.array([-0.22389266154, 0.79983062746, 0.55690487831])
v1=np.array([-0.91397805644, 0.02609403989, -0.40492342125])
x1=np.array([-2161414.92636, -3834695.17889, 4600350.22664])
LIGO1=Detector(x1,u1,v1,RT_psd)
L_detector=4000
c = 299792458.0
freq=0

def T_rt(f, L, ndotu):
    """Eq. (5.27). ndotu has shape (360, 181)."""
    x = f * L / c                       # f L / c  (dimensionless)
    # np.sinc(y) = sin(pi y)/(pi y), so sinc(pi*x*(1+ndotu)) -> np.sinc(x*(1+ndotu))
    term1 = np.exp(-1j*np.pi*x*(1 - ndotu)) * np.sinc(x*(1 + ndotu))
    term2 = np.exp(+1j*np.pi*x*(1 + ndotu)) * np.sinc(x*(1 - ndotu))
    return (L/c) * np.exp(-2j*np.pi*x) * (term1 + term2)

def R_timing(f, L, u, v, N):
    """Eq. (5.30): R^ab_timing, shape (360, 181, 3, 3). u, v are unit 3-vectors."""
    u = np.asarray(u, float); v = np.asarray(v, float)
    Tu = T_rt(f, L, N @ u)              # (360, 181)
    Tv = T_rt(f, L, N @ v)
    uu = np.outer(u, u)                 # (3, 3)
    vv = np.outer(v, v)
    return 0.5 * (Tu[..., None, None]*uu - Tv[..., None, None]*vv)

n_stack = np.stack([n.x, n.y, n.z], axis=-1)
R_ab = R_timing(freq, L_detector, u1, v1, n_stack)

# Polarization response, Eq. (5.10): R^A = R^ab e^A_ab  -> shape (360, 181)
R_plus  = np.einsum('...ab,...ab->...', R_ab, e_plus)
R_cross = np.einsum('...ab,...ab->...', R_ab, e_cross)

# Strain response, Eq. (5.31): divide by 2L/c
R_strain_plus = R_plus / (2*L_detector/c)
R_strain_cross = R_cross / (2*L_detector/c)
R_strain_1 = R_ab / (2*L_detector/c)
R_strain_2 = np.sqrt(R_strain_plus ** 2 + R_strain_cross ** 2)

# Source - https://stackoverflow.com/a/46069458
# Posted by ImportanceOfBeingErnest
# Retrieved 2026-11-18, License - CC BY-SA 3.0
N=500
lon = np.linspace(-np.pi, np.pi,500)
lat = np.linspace(-np.pi/2, np.pi/2,500)
lon= np.tile(lon, N)
lat =np.repeat(lat, N)
print(lon[1:10])
print(lat[1:10])

fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(111, projection='mollweide')

sc = ax.scatter(
    lon,
    lat,
    c=np.real(R_strain_plus),
    cmap='viridis',
    s=10
)
ax.grid(True)
plt.colorbar(sc, label='Value')

ax.set_title('Mollweide Projection-R_plus')

plt.show()

plt.show()
