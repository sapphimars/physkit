# Physkit - A Python toolkit for constants, unit conversions, and equations.
# Copyright (C) 2024 sapphimars
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.


"""
constants.py

Fundamental and derived physical constants, using 2022 CODATA (SI units),
plus selected particle masses from PDG, and astrophysical/cosmological
parameters from external references (Planck data, etc.).
"""

from .config import get_ureg


class Constants:
    ##########################################################################
    #  UNIVERSAL CONSTANTS
    ##########################################################################
    @property
    def c(self):
        # Speed of light in vacuum (exact by definition)
        val = 299792458.0 * get_ureg()("m / s")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    @property
    def G(self):
        # Newtonian constant of gravitation
        val = 6.67430e-11 * get_ureg()("m^3 / (kg * s^2)")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    @property
    def h(self):
        # Planck constant (exact by definition)
        val = 6.62607015e-34 * get_ureg()("J * s")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    @property
    def h_bar(self):
        # Reduced Planck constant ħ = h / (2π)
        val = 1.054571817e-34 * get_ureg()("J * s")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    ##########################################################################
    #  ELECTROMAGNETIC CONSTANTS
    ##########################################################################

    @property
    def mu_0(self):
        # Vacuum magnetic permeability µ0 (no longer exact post-2019)
        val = 1.25663706212e-6 * get_ureg()("N / A^2")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    @property
    def epsilon_0(self):
        # Vacuum electric permittivity ε0 (no longer exact post-2019)
        val = 8.8541878128e-12 * get_ureg()("F / m")  # (F = farad)
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    @property
    def Z_0(self):
        # Characteristic impedance of vacuum Z0
        val = 376.730313412 * get_ureg()("ohm")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    @property
    def e(self):
        # Elementary charge (exact by definition)
        val = 1.602176634e-19 * get_ureg()("C")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    ##########################################################################
    #  PARTICLE MASSES
    ##########################################################################
    @property
    def m_e(self):
        # Electron mass
        val = 9.1093837015e-31 * get_ureg()("kg")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    @property
    def m_p(self):
        # Proton mass
        val = 1.67262192369e-27 * get_ureg()("kg")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    @property
    def m_n(self):
        # Neutron mass
        val = 1.67492749804e-27 * get_ureg()("kg")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    @property
    def m_nu(self):
        # Muon mass
        val = 1.883531627e-28 * get_ureg()("kg")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    @property
    def m_tau(self):
        # Tau mass (approx)
        val = 3.16754e-27 * get_ureg()("kg")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    @property
    def m_u(self):
        # Atomic mass constant (1 u in kg)
        val = 1.66053906660e-27 * get_ureg()("kg")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    @property
    def m_alpha(self):
        # Alpha particle mass
        val = 6.6446573450e-27 * get_ureg()("kg")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    @property
    def m_h(self):
        # Helion mass (helium-3 nucleus)
        val = 5.0064127862e-27 * get_ureg()("kg")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    ##########################################################################
    #  OTHER FUNDAMENTAL CONSTANTS
    ##########################################################################
    @property
    def N_A(self):
        # Avogadro constant (exact by definition)
        val = 6.02214076e23 * get_ureg()("1 / mol")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    @property
    def k_b(self):
        # Boltzmann constant (exact by definition)
        val = 1.380649e-23 * get_ureg()("J / K")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    @property
    def r_e(self):
        # Classical electron radius (m)
        val = 2.8179403262e-15 * get_ureg()("m")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    ##########################################################################
    #  DERIVED FUNDAMENTAL CONSTANTS
    ##########################################################################
    @property
    def alpha(self):
        # Fine-structure constant (dimensionless)
        val = 0.0072973525693
        return val, ""

    @property
    def mu_B(self):
        # Bohr magneton (J/T)
        val = 9.2740100783e-24 * get_ureg()("J/T")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    @property
    def mu_N(self):
        # Nuclear magneton (J/T)
        val = 5.0507837461e-27 * get_ureg()("J/T")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    @property
    def a0(self):
        # Bohr radius (m)
        val = 5.29177210903e-11 * get_ureg()("m")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    @property
    def R_inf(self):
        # Rydberg constant (1/m)
        val = 10973731.568160 * get_ureg()("1/m")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    ##########################################################################
    #  MAGNETIC MOMENTS, G-FACTORS, MOLAR PLANCK CONSTANT
    ##########################################################################
    @property
    def mu_e(self):
        # Electron magnetic moment (magnitude) [J/T] (2022 CODATA)
        val = 9.2847647043e-24 * get_ureg()("J/T")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    @property
    def mu_p(self):
        # Proton magnetic moment [J/T]
        val = 1.41060679736e-26 * get_ureg()("J/T")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    @property
    def mu_n(self):
        # Neutron magnetic moment (magnitude) [J/T]
        val = 9.6623651e-27 * get_ureg()("J/T")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    @property
    def g_e(self):
        # Electron g-factor (dimensionless)
        val = 2.00231930436256
        return val, ""

    @property
    def g_p(self):
        # Proton g-factor (dimensionless)
        val = 5.5856946893
        return val, ""

    @property
    def g_n(self):
        # Neutron g-factor (dimensionless)
        val = -3.8260837  # negative sign is typical convention
        return val, ""

    @property
    def mu_mu(self):
        # Muon magnetic moment (magnitude) [J/T]
        val = 4.49044830e-26 * get_ureg()("J/T")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    @property
    def g_mu(self):
        # Muon g-factor (dimensionless)
        val = 2.0023318418
        return val, ""

    @property
    def molar_planck(self):
        # Molar Planck constant (N_A * h) [J·s/mol]
        val = 3.990312717628e-10 * get_ureg()("J * s / mol")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    ##########################################################################
    #  THERMODYNAMIC
    ##########################################################################
    @property
    def sigma_SB(self):
        # Stefan-Boltzmann constant (sigma) [W/(m^2·K^4)]
        val = 5.670374419e-8 * get_ureg()("W / (m^2 * K^4)")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    @property
    def R_gas(self):
        # Ideal gas constant (R) [J/(mol·K)]
        val = 8.31446261815324 * get_ureg()("J / (mol * K)")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    @property
    def Rydberg_energy(self):
        # Rydberg energy (~13.605693122 eV in joules)
        val = 2.1798723611035e-18 * get_ureg()("J")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    @property
    def atm(self):
        # 1 atm pressure (exact in SI)
        val = 101325.0 * get_ureg()("Pa")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    ##########################################################################
    #  HADRON MASSES
    ##########################################################################
    @property
    def m_pion_charged(self):
        # Charged pion mass (pi±) ~139.57039 MeV/c^2 => ~2.49e-28 kg
        val = 2.48832e-28 * get_ureg()("kg")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    @property
    def m_pion_neutral(self):
        # Neutral pion mass (pi0) ~134.9770 MeV/c^2 => ~2.41e-28 kg
        val = 2.405e-28 * get_ureg()("kg")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    @property
    def m_kaon_charged(self):
        # Charged kaon mass (K±) ~493.677 MeV/c^2 => ~8.80e-28 kg
        val = 8.80e-28 * get_ureg()("kg")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    @property
    def m_kaon_neutral(self):
        # Neutral kaon mass (K0) ~497.611 MeV/c^2 => ~8.87e-28 kg
        val = 8.87e-28 * get_ureg()("kg")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    ##########################################################################
    #  COMPTON WAVELENGTHS
    ##########################################################################
    @property
    def lambda_compton_e(self):
        # Electron Compton wavelength ~2.42631023867e-12 m
        val = 2.42631023867e-12 * get_ureg()("m")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    @property
    def lambda_compton_p(self):
        # Proton Compton wavelength ~1.32140985539e-15 m
        val = 1.32140985539e-15 * get_ureg()("m")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    @property
    def lambda_compton_n(self):
        # Neutron Compton wavelength ~1.31959090581e-15 m
        val = 1.31959090581e-15 * get_ureg()("m")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    ##########################################################################
    # SCATTERING CROSS SECTIONS
    ##########################################################################
    @property
    def sigma_T(self):
        # Thomson cross section (2022 CODATA) [m^2]
        val = 6.6524587321e-29 * get_ureg()("m^2")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    ##########################################################################
    # ASTRONOMICAL / ASTROPHYSICAL CONSTANTS
    ##########################################################################
    @property
    def M_earth(self):
        # Earth mass (approx)
        val = 5.97219e24 * get_ureg()("kg")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    @property
    def M_sun(self):
        # Solar mass (approx)
        val = 1.98847e30 * get_ureg()("kg")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    @property
    def L_sun(self):
        # Solar luminosity (approx)
        val = 3.828e26 * get_ureg()("W")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    @property
    def AU(self):
        # Astronomical Unit (exact, 2012 definition: 1 AU = 149,597,870,700 m)
        val = 1.49597870700e11 * get_ureg()("m")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    ##########################################################################
    # COSMOLOGICAL CONSTANTS
    ##########################################################################
    @property
    def H_0(self):
        # Hubble constant ~67.66 (km/s)/Mpc => ~2.192e-18 1/s
        val = 2.192e-18 * get_ureg()("1/s")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    @property
    def rho_crit(self):
        # Critical density of the Universe (~8.6e-27 kg/m^3)
        val = 8.6e-27 * get_ureg()("kg / m^3")
        converted = val.to_base_units()
        return converted.magnitude, f"{converted.units:~}"

    @property
    def pi(self):
        # pi
        val = 3.1415926535897932384626433
        return val, ""


constants = Constants()
