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


import pint
from .custom_units import define_custom_units

ureg = pint.UnitRegistry(system="SI")  # Set default to SI
define_custom_units(ureg)  # Define units in SI


def set_default(unit_system: str) -> None:
    """
    Change global default unit system.

    Args:
        unit_system (str): unit system to change to

    Example:
        physkit.set_default("cgs") -> change to cgs unit system
    """
    global ureg
    ureg = pint.UnitRegistry(system=unit_system)
    define_custom_units(ureg)


def get_ureg():
    return ureg
