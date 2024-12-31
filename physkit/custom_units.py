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


def define_custom_units(ureg):
    """
    Define custom units that Pint doesn't include by default.
    Call this AFTER creating the UnitRegistry in config.py.

    Args:
        ureg (pint.UnitRegistry): global ureg object from config.py
    """
    # Jansky (Jy)
    ureg.define("jansky = 1e-26 watt / meter^2 / hertz = Jy")

    # Solar luminosity: 1 L_sun = 3.828e26 W
    ureg.define("solar_luminosity = 3.828e26 watt = L_sun")

    # Solar mass: 1 M_sun = 1.98847e30 kg
    ureg.define("solar_mass = 1.98847e30 kg = M_sun")
