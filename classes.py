import numpy as np
import warnings
import math
from scipy.special import spherical_jn
# WGS-84 ellipsoid parameters
_A = 6378137.0          # semi-major axis [m]
_B = 6356752.314245     # semi-minor axis [m]

H0=2.2*10**(-18)
years_to_seconds=60*60*24*365.25
def _R(phi):
    """Local Earth radius R(phi)"""
    return _A**2 / np.sqrt(_A**2 * np.cos(phi) ** 2 + _B**2 * np.sin(phi) ** 2)
def _enu_basis(phi, lam):
    """East (lambda_hat), North (phi_hat), Up (h_hat) unit vectors,
    expressed in the geocentric (i, j, k) frame"""
    lambda_hat = np.array([-np.sin(lam), np.cos(lam), 0.0])
    phi_hat = np.array(
        [-np.sin(phi) * np.cos(lam), -np.sin(phi) * np.sin(lam), np.cos(phi)]
    )
    h_hat = np.array(
        [np.cos(phi) * np.cos(lam), np.cos(phi) * np.sin(lam), np.sin(phi)]
    )
    return lambda_hat, phi_hat, h_hat
def _vertex_and_arms(phi, lam, h, psi_u, psi_v):
    """
    Compute the vertex position and arm unit vectors given latitude
    phi, longitude lam (radians), elevation h (meters), and the
    north-of-east azimuths psi_u/psi_v (radians) of the two arms.
    """
    R = _R(phi)

    x = (R + h) * np.cos(phi) * np.cos(lam)
    y = (R + h) * np.cos(phi) * np.sin(lam)
    z = ((_B**2 / _A**2) * R + h) * np.sin(phi)
    x_ijk = np.array([x, y, z])

    lambda_hat, phi_hat, _h_hat = _enu_basis(phi, lam)

    # eps = 0 -> arms lie entirely in the tangent plane
    u = np.cos(psi_u) * lambda_hat + np.sin(psi_u) * phi_hat
    v = np.cos(psi_v) * lambda_hat + np.sin(psi_v) * phi_hat
    return x_ijk, u, v
def _lal_to_east_azimuth(lal_azimuth_rad):
    """
    Azimuth converter
    """
    return (2.5 * np.pi - lal_azimuth_rad) % (2 * np.pi)

class Detector:
    def __init__(self, x_ijk, arm1_ijk, arm2_ijk, psd):
        """
        Parameters
        ----------
        x_ijk : 3-vector
            Detector position expressed in ijk
        arm1_ijk, arm2_ijk : 3-vectors
            Detector arms expressed in ijk system
        """
        self.x=x_ijk
        self.u = arm1_ijk
        self.v = arm2_ijk
        self.psd=psd
        self._validate()
        self.D = self._detector_tensor()
    def _detector_tensor(self):
        """D^ab = 1/2 (u^a u^b - v^a v^b), the symmetric trace-free
        detector tensor in the small-antenna limit."""
        return 0.5 * (np.outer(self.u, self.u) - np.outer(self.v, self.v))
    def _validate(self, tol=1e-6):
        if self.x.shape != (3,) or self.u.shape != (3,) or self.v.shape != (3,):
            raise ValueError("x, u, v must each be 3-vectors.")
        for label in ("u", "v"):
            vec = getattr(self, label)
            norm = np.linalg.norm(vec)
            if abs(norm - 1.0) > tol:
                warnings.warn(
                    f"arm vector '{label}' has norm {norm:.6f}; normalizing."
                )
                setattr(self, label, vec / norm)

class LAL_Detector(Detector):
    """
    Detector initialized with lattitude, longitude, elevation and arm azimuths
    """

    def __init__(
        self,
        longitude_rad,
        latitude_rad,
        elevation_si,
        arm_x_azimuth_rad,
        arm_y_azimuth_rad,
        psd
    ):
        """
        Parameters
        ----------
        longitude_rad, latitude_rad : float
            Vertex longitude / latitude, in radians
        elevation_si : float
            Vertex elevation above the WGS-84 ellipsoid, in meters
        arm_x_azimuth_rad, arm_y_azimuth_rad : float
            LAL-convention azimuths of the x- and y-arms, in radians

        """
        psi_u = _lal_to_east_azimuth(arm_x_azimuth_rad)
        psi_v = _lal_to_east_azimuth(arm_y_azimuth_rad)

        x_ijk, u, v = _vertex_and_arms(
            latitude_rad, longitude_rad, elevation_si, psi_u, psi_v
        )

        super().__init__(x_ijk, u, v, psd)

class DetectorPair:
    def __init__(self, detector1,detector2):
        self.detector1 = detector1
        self.detector2 = detector2
        self.deltaX=self.detector1.x-self.detector2.x
        self.s=np.linalg.norm(self.deltaX)
        self.s_hat=self.deltaX/self.s
        self.alpha=0

    def calculate_alpha(self,f):
        c=3*10**8
        self.alpha=2*np.pi*f*self.s/c

    def get_ABCDE(self):
        J0 = spherical_jn(0, self.alpha)
        J1 = spherical_jn(1, self.alpha)
        J2 = spherical_jn(2, self.alpha)
        A=1/(2*self.alpha**2)*(-5*self.alpha**2*J0+10*self.alpha*J1+5*J2)
        B = 1 / (2 * self.alpha ** 2) * (5 * self.alpha ** 2 * J0 + -10 * self.alpha * J1 + 5 * J2)
        C = 1 / (2 * self.alpha ** 2) * (5 * self.alpha ** 2 * J0 + -10 * self.alpha * J1 + -25 * J2)
        D = 1 / (2 * self.alpha ** 2) * (-5 * self.alpha ** 2 * J0 + 20 * self.alpha * J1 + -25 * J2)
        E = 1 / (2 * self.alpha ** 2) * (5 * self.alpha ** 2 * J0 + -50 * self.alpha * J1 + 175 * J2)

        return A,B,C,D,E

    def overlap_function(self):
        A,B,C,D,E=self.get_ABCDE()

        Tr1=np.trace(self.detector1.D)
        Tr2=np.trace(self.detector2.D)
        term1=A*Tr1*Tr2

        # D1^ab D2_ab  (full double contraction)
        contraction12 = np.sum(self.detector1.D * self.detector2.D)
        term2=2*B*contraction12

        # s^T D1 s  and  s^T D2 s
        s_D1_s = self.s_hat @ self.detector1.D @ self.s_hat
        s_D2_s = self.s_hat @ self.detector2.D @ self.s_hat
        term3=C * (Tr1 * s_D2_s + Tr2 * s_D1_s)

        # D1^ab D2_a^c s_b s_c = s^T (D1 @ D2) s
        s_D1D2_s = self.s_hat @ (self.detector1.D @ self.detector2.D) @ self.s_hat
        term4=4 * D * s_D1D2_s

        term5=E * s_D1_s * s_D2_s
        gamma = term1 + term2 + term3 + term4 + term5
        return gamma

class DetectorPairVec:
    def __init__(self, detector1, detector2):
        self.detector1 = detector1
        self.detector2 = detector2
        self.deltaX = self.detector1.x - self.detector2.x
        self.s = np.linalg.norm(self.deltaX)
        self.s_hat = self.deltaX / self.s

        # --- geometry-only constants, computed ONCE, independent of frequency ---
        self.Tr1 = np.trace(self.detector1.D)
        self.Tr2 = np.trace(self.detector2.D)
        self.contraction12 = np.sum(self.detector1.D * self.detector2.D)
        self.s_D1_s = self.s_hat @ self.detector1.D @ self.s_hat
        self.s_D2_s = self.s_hat @ self.detector2.D @ self.s_hat
        self.s_D1D2_s = self.s_hat @ (self.detector1.D @ self.detector2.D) @ self.s_hat

    def calculate_alpha_vec(self, freqs):
        c = 3e8
        return 2 * np.pi * freqs * self.s / c

    def get_ABCDE_vec(self, alpha):
        # scipy's spherical_jn is elementwise-vectorized over z already
        J0 = spherical_jn(0, alpha)
        J1 = spherical_jn(1, alpha)
        J2 = spherical_jn(2, alpha)

        alpha2 = alpha ** 2
        inv = 1.0 / (2 * alpha2)   # NOTE: singular at alpha=0, see caveat below

        A = inv * (-1 * alpha2 * J0 + 2 * alpha * J1 + 1 * J2)
        B = inv * ( 1 * alpha2 * J0 - 2 * alpha * J1 + 1 * J2)
        C = inv * ( 1 * alpha2 * J0 - 2 * alpha * J1 - 5 * J2)
        D = inv * (-1 * alpha2 * J0 + 2 * alpha * J1 - 5 * J2)
        E = inv * ( 1 * alpha2 * J0 - 10 * alpha * J1 + 35 * J2)
        return A, B, C, D, E

    def overlap_function_vec(self, freqs):
        alpha = self.calculate_alpha_vec(freqs)
        A, B, C, D, E = self.get_ABCDE_vec(alpha)

        term1 = A * self.Tr1 * self.Tr2
        term2 = 2 * B * self.contraction12
        term3 = C * (self.Tr1 * self.s_D2_s + self.Tr2 * self.s_D1_s)
        term4 = 4 * D * self.s_D1D2_s
        term5 = E * self.s_D1_s * self.s_D2_s

        return term1 + term2 + term3 + term4 + term5

class DetectorNetwork:
    def __init__(self, detectors, nr_of_detectors, frequencies, nr_of_frequencies):
        self.detectors=detectors
        self.nr_of_detectors=nr_of_detectors
        self.frequencies=frequencies
        self.nr_of_frequencies=nr_of_frequencies
        self.R_eff=np.zeros(self.nr_of_frequencies)
        self.S_eff = np.zeros(self.nr_of_frequencies)
        self.Omega_eff = np.zeros(self.nr_of_frequencies)
        self.set_R_eff()
        self.set_S_eff()
        self.set_Omega_eff()
    def set_R_eff(self):
        if self.nr_of_detectors==2:
            pair=DetectorPairVec(self.detectors[0],self.detectors[1])
            self.R_eff=np.sqrt(pair.overlap_function_vec(self.frequencies) ** 2)
        elif self.nr_of_detectors>2:
            product = np.ones(self.nr_of_frequencies)
            for i in range(self.nr_of_detectors-1):
                for j in range(i+1,self.nr_of_detectors):
                    pair=DetectorPairVec(self.detectors[i],self.detectors[j])
                    product=product*pair.overlap_function_vec(self.frequencies)**2
            self.R_eff=np.sqrt(product)
        else:
            return
    def set_S_eff(self):
        product=np.ones(self.nr_of_frequencies)
        for detector in self.detectors:
            product=product * detector.psd
        self.S_eff=np.sqrt(product)/self.R_eff
    def set_Omega_eff(self):
        self.Omega_eff= (2 * np.pi ** 2 * self.frequencies ** 3 * self.S_eff) / (3 * H0 ** 2)
    def get_Omega_eff(self,years,df=0.25):
        if years==0:
            return self.Omega_eff
        elif years>0:
            return self.Omega_eff/np.sqrt(2*years*years_to_seconds*df)
        else:
            return np.zeros(self.nr_of_frequencies)
    def Power_Integration(self,betas,years=1,rho=1,fref=100):
        Omega_betas=np.zeros(np.size(betas))
        for i in range(np.size(betas)):
            Omega_betas[i]=(rho / math.sqrt(2*years*years_to_seconds) /
                            np.sqrt(np.trapezoid((self.frequencies / fref) ** (2 * betas[i]) / self.Omega_eff ** 2,
                                               self.frequencies)))
        Omega_PI=np.zeros(self.nr_of_frequencies)
        for j in range(self.nr_of_frequencies):
            Omega_PI[j]=np.max(Omega_betas * (self.frequencies[j] / fref) ** betas)
        return(Omega_PI,Omega_betas)