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
conversions.py

Convenience functions for unit conversions, using the shared Pint
UnitRegistry from config.py.

The user-facing API only takes/returns floats and strings:
    - convert_unit(value, from_unit, to_unit)
    - convert_to_base_units(value, from_unit)
    - convert_unit_system(value, from_unit, system_name="cgs")

Under the hood, we use Pint's UnitRegistry (from config.py) to do
all conversions. The user never interacts with Pint objects.
"""

import pint
from .config import get_ureg  # , set_default


def convert_unit(value: float, from_unit: str, to_unit: str) -> float:
    """
    Convert a numeric value from `from_unit` to `to_unit`.

    Args:
        value (float): The numeric value to convert.
        from_unit (str): The current unit of `value`.
        to_unit (str): The target unit.

    Returns:
        float: The numeric value converted to the `to_unit`.

    Example:
        convert_unit(5, "erg", "J") -> 5e-7
        convert_unit(1, "meter", "foot") -> ~3.28084
    """
    ureg = get_ureg()
    quantity = value * ureg(from_unit)
    converted = quantity.to(to_unit)
    return converted.magnitude


def convert_to_base_units(value: float, from_unit: str) -> tuple[float, str]:
    """
    Convert a numeric value from `from_unit` to the *default base units*
    of the shared UnitRegistry. Typically this is SI base units unless
    you've called set_default_system(...) to change it.

    Args:
        value (float): input value
        from_unit (str): input unit

    Returns:
        tuple[float, str]: output value and its unit

    Example:
        # If default system is SI:
        convert_to_base_units(39.37, "inch") -> ~1.0 (meter)
    """
    ureg = get_ureg()
    quantity = value * ureg(from_unit)
    base_converted = quantity.to_base_units()
    return base_converted.magnitude, f"{base_converted.units:~}"


def convert_unit_system(
    value: float, from_unit: str, system_name: str = "cgs"
) -> tuple[float, str]:
    """
    Convert a numeric value from `from_unit` in the shared default registry
    to the base units of a different unit system (e.g., 'cgs', 'imperial', etc.)
    without changing the global default system.

    Creates a temporary pint.UnitRegistry with the specified system,
    converts `value` to base SI, then interprets that SI quantity in
    the temporary system's base units.

    Args:
        value (float): input value
        from_unit (str): input unit
        system_name (str, optional): output unit system. Defaults to "cgs".

    Returns:
        tuple[float, str]: output value and its unit

    Example:
        # Suppose global default is SI
        # We convert 1 meter (in SI) to the base units of cgs -> centimeter
        convert_unit_system(1, "meter", system_name="cgs") -> 100.0 (centimeter)
    """

    # Convert the input to base SI using the shared registry
    shared_ureg = get_ureg()
    si_quantity = (value * shared_ureg(from_unit)).to_base_units()

    # Create a temporary registry in the chosen system
    temp_ureg = pint.UnitRegistry(system=system_name)

    # Build the quantity in the temporary registry

    magnitude_si = si_quantity.magnitude
    unit_si_str = str(si_quantity.u)

    new_qty_in_temp = magnitude_si * temp_ureg(unit_si_str)

    # Convert to the base units of the new system
    base_in_temp_system = new_qty_in_temp.to_base_units()

    return base_in_temp_system.magnitude, f"{base_in_temp_system.units:~}"


class UnitConverter:
    lengths = {
        "meter": 1.0,
        "inch": 0.0254,
        "foot": 0.3048,
        "yard": 0.9144,
        "mile": 1609.34,
        "eV^-1": 1.97327e-7,
    }
    masses = {
        "gram": 1.0,
        "pound": 453.59,
        "ounce": 28.35,
        "solar_mass": 1.98847e33,
        "eV": 1.78266192e-33,
    }
    times = {
        "second": 1.0,
        "minute": 60.0,
        "hour": 3600.0,
        "eV^-1": 6.582119e-16,
    }
    energies = {
        "joule": 1.0,
        "calorie": 4.184,
        "watt-hour": 3.6e9,
        "BTU": 1055,
        "eV": 1.602176634e-19,
    }

    @classmethod
    def convert_length(cls, value: float, from_unit: str, to_unit: str) -> float:
        base_value = value * cls.lengths[from_unit]
        converted = base_value / cls.lengths[to_unit]
        return converted

    @classmethod
    def convert_mass(cls, value: float, from_unit: str, to_unit: str) -> float:
        base_value = value * cls.masses[from_unit]
        converted = base_value / cls.masses[to_unit]
        return converted

    @classmethod
    def convert_time(cls, value: float, from_unit: str, to_unit: str) -> float:
        base_value = value * cls.times[from_unit]
        converted = base_value / cls.times[to_unit]
        return converted

    @classmethod
    def convert_energy(cls, value: float, from_unit: str, to_unit: str) -> float:
        base_value = value * cls.energy[from_unit]
        converted = base_value / cls.energy[to_unit]
        return converted
