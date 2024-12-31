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
