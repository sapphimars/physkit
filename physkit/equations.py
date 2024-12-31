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


from .config import get_ureg
from .constants import constants


class Equations:

    @staticmethod
    def gravitational_radius(
        mass_value: float, mass_unit: str = "kg", out_unit: str = "m"
    ) -> tuple[float, str]:
        """
        Compute the gravitational radius for a given mass.
        r_g = G * M / c^2

        Args:
            mass_value (float): input mass
            mass_unit (str, optional): unit for input mass. Defaults to "kg".
            out_unit (str, optional): unit for output radius. Defaults to "m".

        Returns:
            (float, str): output radius and its unit
        """

        # Get the shared registry
        ureg = get_ureg()

        # Get pint Quantity for mass
        mass_qty = mass_value * ureg(mass_unit)

        # Retrieve constants from the constants module
        G_val, G_unit = constants.G
        c_val, c_unit = constants.c

        # Get Pint Quantities for constants
        G_qty = G_val * ureg(G_unit)
        c_qty = c_val * ureg(c_unit)

        r_g_qty = G_qty * mass_qty / (c_qty**2)

        # Convert to output unit
        r_g_converted = r_g_qty.to(out_unit)

        return r_g_converted.magnitude, f"{r_g_converted.units:~}"

    @staticmethod
    def schwartzschild_radius(
        mass_value: float, mass_unit: str = "kg", out_unit: str = "m"
    ) -> tuple[float, str]:
        """
        Compute the Schwartzschild radius for a given mass.

        Args:
            mass_value (float): input mass
            mass_unit (str, optional): unit for input mass. Defaults to "kg".
            out_unit (str, optional): unit for output radius. Defaults to "m".

        Returns:
            tuple[float, str]: output radius and its unit
        """
        mag, unit_str = Equations.gravitational_radius(mass_value, mass_unit, out_unit)
        mag *= 2
        return mag, unit_str

    @staticmethod
    def eddington_luminosity(
        mass_value: float, mass_unit: str = "kg", out_unit: str = "W"
    ) -> tuple[float, str]:
        """
        Compute the Eddington luminosity for a given mass.

        Args:
            mass_value (float): input mass
            mass_unit (str, optional): unit for input mass. Defaults to "kg".
            out_unit (str, optional): unit for output luminosity. Defaults to "W".

        Returns:
            tuple[float, str]: output luminosity and its unit
        """
        ureg = get_ureg()
        mass_qty = mass_value * ureg(mass_unit)

        pi, _ = constants.pi
        G_val, G_unit = constants.G  # e.g. (6.6743e-11, "m^3 / kg / s^2")
        c_val, c_unit = constants.c  # e.g. (3.0e8, "m / s")
        m_p_val, m_p_unit = constants.m_p
        sigma_T_val, sigma_T_unit = constants.sigma_T

        G_qty = G_val * ureg(G_unit)
        c_qty = c_val * ureg(c_unit)
        m_p_qty = m_p_val * ureg(m_p_unit)
        sigma_T_qty = sigma_T_val * ureg(sigma_T_unit)

        L_edd_qty = 4 * pi * G_qty * m_p_qty * c_qty * mass_qty / sigma_T_qty

        L_edd_converted = L_edd_qty.to(out_unit)

        return L_edd_converted.magnitude, f"{L_edd_converted.units:~}"


equations = Equations()
