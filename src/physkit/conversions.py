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

import re
import pint
from .config import get_ureg  # , set_default
from .constants import constants


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
    # Base conversion dictionaries
    lengths = {
        "meter": 1.0,
        "inch": 0.0254,
        "foot": 0.3048,
        "cm": 0.01,
        "mile": 1609.34,
        "yard": 0.9144,
        "angstrom": 1e-10,
        "eV^-1": 1.97327e-7,
    }
    masses = {
        "gram": 1.0,
        "kg": 1000.0,
        "pound": 453.59,
        "ounce": 28.35,
        "solar_mass": 1.98847e33,
        "eV": 1.78266192e-33,
    }
    times = {
        "second": 1.0,
        "minute": 60.0,
        "hour": 3600.0,
        "day": 86400.0,
        "year": 31536000.0,
        "eV^-1": 6.582119e-16,
    }
    energies = {
        "joule": 1.0,
        "calorie": 4.184,
        "watt-hour": 3600.0,
        "BTU": 1055.0,
        "erg": 1e-7,
        "eV": 1.602176634e-19,
    }

    # Custom compound units defined in terms of base units
    compound_units = {
        "hertz": {
            "definition": "1/second",
            "dimension": "frequency",
            "conversion_factor": 1.0,  # In terms of 1/second
        },
        "watt": {
            "definition": "joule/second",
            "dimension": "energy/time",
            "conversion_factor": 1.0,  # In terms of joule/second
        },
    }

    # Unit aliases - maps alternative forms to canonical form
    unit_aliases = {
        # Length
        "m": "meter",
        "meters": "meter",
        "metre": "meter",
        "metres": "meter",
        "in": "inch",
        "inches": "inch",
        "ft": "foot",
        "feet": "foot",
        "yd": "yard",
        "yards": "yard",
        "mi": "mile",
        "miles": "mile",
        "centimeter": "cm",
        "centimeters": "cm",
        "centimetre": "cm",
        "centimetres": "cm",
        "Å": "angstrom",
        "Angstrom": "angstrom",
        "angstroms": "angstrom",
        # Mass
        "g": "gram",
        "grams": "gram",
        "gramme": "gram",
        "grammes": "gram",
        "kilogram": "kg",
        "kilograms": "kg",
        "kilogramme": "kg",
        "kilogrammes": "kg",
        "lb": "pound",
        "pounds": "pound",
        "lbs": "pound",
        "oz": "ounce",
        "ounces": "ounce",
        "solar mass": "solar_mass",
        "M_sun": "solar_mass",
        "M☉": "solar_mass",
        # Time
        "s": "second",
        "seconds": "second",
        "sec": "second",
        "min": "minute",
        "minutes": "minute",
        "h": "hour",
        "hours": "hour",
        "hr": "hour",
        "hrs": "hour",
        "d": "day",
        "days": "day",
        "yr": "year",
        "years": "year",
        "y": "year",
        # Energy
        "J": "joule",
        "joules": "joule",
        "cal": "calorie",
        "calories": "calorie",
        "Wh": "watt-hour",
        "watt hours": "watt-hour",
        "electron volt": "eV",
        "electronvolt": "eV",
        "electron-volt": "eV",
        "electronvolts": "eV",
        "electron volts": "eV",
        "electron-volts": "eV",
        # Compound units
        "Hz": "hertz",
        "hz": "hertz",
        "Hertz": "hertz",
        "W": "watt",
        "watts": "watt",
        "Watt": "watt",
        "Watts": "watt",
        "Jy": "jansky",
        "jy": "jansky",
        "Jansky": "jansky",
        "janskys": "jansky",
        "Janskys": "jansky",
    }

    # SI prefixes and their multipliers - both symbol and full word forms
    si_prefixes = {
        # Sub-units - symbols
        "y": 1e-24,  # yocto
        "z": 1e-21,  # zepto
        "a": 1e-18,  # atto
        "f": 1e-15,  # femto
        "p": 1e-12,  # pico
        "n": 1e-9,  # nano
        "μ": 1e-6,  # micro (also support u as alias for micro)
        "u": 1e-6,  # micro (ascii version)
        "m": 1e-3,  # milli
        "c": 1e-2,  # centi
        "d": 1e-1,  # deci
        # Sub-units - full words
        "yocto": 1e-24,
        "zepto": 1e-21,
        "atto": 1e-18,
        "femto": 1e-15,
        "pico": 1e-12,
        "nano": 1e-9,
        "micro": 1e-6,
        "milli": 1e-3,
        "centi": 1e-2,
        "deci": 1e-1,
        # Super-units - symbols
        "da": 1e1,  # deca
        "h": 1e2,  # hecto
        "k": 1e3,  # kilo
        "M": 1e6,  # mega
        "G": 1e9,  # giga
        "T": 1e12,  # tera
        "P": 1e15,  # peta
        "E": 1e18,  # exa
        "Z": 1e21,  # zetta
        "Y": 1e24,  # yotta
        # Super-units - full words
        "deca": 1e1,
        "hecto": 1e2,
        "kilo": 1e3,
        "mega": 1e6,
        "giga": 1e9,
        "tera": 1e12,
        "peta": 1e15,
        "exa": 1e18,
        "zetta": 1e21,
        "yotta": 1e24,
    }

    # Set of base units that accept SI prefixes
    prefix_compatible_units = {
        # Base units
        "meter",
        "gram",
        "second",
        "joule",
        "watt-hour",
        "eV",
        "m",
        "g",
        "s",
        "J",
        "Wh",
        # Derived units that can take prefixes
        "hertz",
        "Hz",
        "watt",
        "W",
        "jansky",
        "Jy",
    }

    # Regular expressions for parsing units
    unit_regex = re.compile(r"([a-zA-Z_μ]+)(\^-?\d+)?")
    # Updated prefix regex to match full word prefixes too
    prefix_regex = re.compile(
        r"(yocto|zepto|atto|femto|pico|nano|micro|milli|centi|deci|deca|hecto|kilo|mega|giga|tera|peta|exa|zetta|yotta|[yzafpnμumcdhkMGTPEZY])?([a-zA-Z_]+)"
    )
    compound_regex = re.compile(r"([^*/^]+)(?:\^(-?\d+))?")

    dimensions = {
        "energy": {"joule", "calorie", "eV", "erg", "BTU", "watt-hour"},
        "length": {"meter", "inch", "foot", "cm", "mile", "yard", "angstrom", "eV^-1"},
        "time": {"second", "minute", "hour", "day", "year", "eV^-1"},
        "mass": {"gram", "kg", "pound", "ounce", "solar_mass", "eV"},
        "frequency": {"hertz"},
        "energy/time": {"watt"},
        "magnetic_flux_density": {"tesla"},
        "radioactivity": {"becquerel", "curie"},
        "area": {"barn"},
        "exposure": {"roentgen"},
        "pressure": {"pascal"},
        "volume_flow_rate": {"cubic_meter_per_second"},
        "mass_flow_rate": {"kilogram_per_second"},
        "temperature": {"kelvin"},
        "force": {"newton"},
        "velocity": {"speed_of_light"},
    }

    @classmethod
    def add_unit_alias(cls, alias, canonical_form):
        """
        Add a new unit alias to the system.

        Args:
            alias (str): The alias/alternative form
            canonical_form (str): The canonical form this maps to
        """
        cls.unit_aliases[alias] = canonical_form

    @classmethod
    def add_compound_unit(
        cls, name, definition, dimension, conversion_factor=1.0, aliases=None
    ):
        """
        Add a new compound unit definition.

        CURRENTLY BROKEN FOR SOME UNITS WHERE THE CONVERSION FACTOR IS NECESSARY
        AND OTHERS WHERE CONVERSION FACTOR BREAKS THINGS. FIX!!!

        Args:
            name (str): Canonical name of the unit
            definition (str): Definition in terms of other units (e.g., "joule/second")
            dimension (str): Physical dimension description (e.g., "energy/time")
            conversion_factor (float): Conversion factor to base SI units
            aliases (list): Optional list of aliases for this unit
        """
        # Add the compound unit definition
        cls.compound_units[name] = {
            "definition": definition,
            "dimension": dimension,
            "conversion_factor": conversion_factor,
        }

        # Add dimension if it doesn't exist
        if dimension not in cls.dimensions:
            cls.dimensions[dimension] = {name}
        else:
            cls.dimensions[dimension].add(name)

        # Add any aliases
        if aliases is not None:
            for alias in aliases:
                cls.add_unit_alias(alias, name)

    @classmethod
    def normalize_unit(cls, unit_str):
        """
        Convert unit to canonical form, handling aliases and prefixes.

        Args:
            unit_str (str): Unit string like "km" or "microsec"

        Returns:
            tuple: (canonical_unit, multiplier)
        """
        # Check if it's a simple alias
        if unit_str in cls.unit_aliases:
            return cls.unit_aliases[unit_str], 1.0

        # Check if it's a defined compound unit
        if unit_str in cls.compound_units:
            return unit_str, cls.compound_units[unit_str]["conversion_factor"]

        # Check for prefix
        prefix_match = cls.prefix_regex.match(unit_str)
        if prefix_match:
            prefix, base_unit = prefix_match.groups()

            # If base unit is an alias, convert to canonical form
            if base_unit in cls.unit_aliases:
                base_unit = cls.unit_aliases[base_unit]

            # If this unit can take a prefix and we have a prefix
            if base_unit in cls.prefix_compatible_units and prefix:
                if prefix in cls.si_prefixes:
                    return base_unit, cls.si_prefixes[prefix]

            # If it's a valid base unit with no prefix
            if base_unit in cls.unit_aliases:
                return cls.unit_aliases[base_unit], 1.0

        # If no match, return as is (might be a compound unit or unknown)
        return unit_str, 1.0

    @classmethod
    def parse_compound_unit(cls, unit_str):
        # Check if this is a predefined compound unit
        if unit_str in cls.compound_units:
            # Use the definition to expand it
            definition = cls.compound_units[unit_str]["definition"]
            conversion_factor = cls.compound_units[unit_str]["conversion_factor"]

            # Parse the definition - but skip the compound unit check to avoid recursion
            parts, def_multiplier = cls._parse_compound_unit_internal(definition)
            return parts, def_multiplier * conversion_factor

        # Otherwise parse it directly
        return cls._parse_compound_unit_internal(unit_str)

    @classmethod
    def _parse_compound_unit_internal(cls, unit_str):
        parts = {"numerator": {}, "denominator": {}}
        multiplier = 1.0

        # Split on division
        if "/" in unit_str:
            num_str, denom_str = unit_str.split("/", 1)

            # Handle numerator
            for unit_part in num_str.split("*"):
                unit_part = unit_part.strip()
                if not unit_part:
                    continue

                # Check for numeric coefficient at the start
                coef_match = re.match(r"^(\d+\.?\d*[eE]?[-+]?\d*)(.+)$", unit_part)
                if coef_match:
                    coef_str, unit_part = coef_match.groups()
                    multiplier *= float(coef_str)

                power_match = re.search(r"(\^-?\d+)$", unit_part)
                if power_match:
                    base_unit = unit_part[: power_match.start()]
                    power = int(power_match.group(0)[1:])
                else:
                    base_unit = unit_part
                    power = 1

                # Normalize the base unit and get any prefix multiplier
                base_unit, prefix_mult = cls.normalize_unit(base_unit)
                multiplier *= prefix_mult**power

                # Add to numerator
                if base_unit in parts["numerator"]:
                    parts["numerator"][base_unit] += power
                else:
                    parts["numerator"][base_unit] = power

            # Handle multiple denominators (e.g., meter/second/hour)
            for denom_part in denom_str.split("/"):
                denom_part = denom_part.strip()
                if not denom_part:
                    continue

                # Check for numeric coefficient at the start
                coef_match = re.match(r"^(\d+\.?\d*[eE]?[-+]?\d*)(.+)$", denom_part)
                if coef_match:
                    coef_str, denom_part = coef_match.groups()
                    multiplier /= float(coef_str)

                power_match = re.search(r"(\^-?\d+)$", denom_part)
                if power_match:
                    base_unit = denom_part[: power_match.start()]
                    power = int(power_match.group(0)[1:])
                else:
                    base_unit = denom_part
                    power = 1

                # Normalize the base unit and get any prefix multiplier
                base_unit, prefix_mult = cls.normalize_unit(base_unit)
                multiplier /= prefix_mult**power

                # Add to denominator
                if base_unit in parts["denominator"]:
                    parts["denominator"][base_unit] += power
                else:
                    parts["denominator"][base_unit] = power
        else:
            # Handle numerator only expressions (e.g., "joule*meter^2")
            for unit_part in unit_str.split("*"):
                unit_part = unit_part.strip()
                if not unit_part:
                    continue

                # Check for numeric coefficient at the start
                coef_match = re.match(r"^(\d+\.?\d*[eE]?[-+]?\d*)(.+)$", unit_part)
                if coef_match:
                    coef_str, unit_part = coef_match.groups()
                    multiplier *= float(coef_str)

                power_match = re.search(r"(\^-?\d+)$", unit_part)
                if power_match:
                    base_unit = unit_part[: power_match.start()]
                    power = int(power_match.group(0)[1:])
                else:
                    base_unit = unit_part
                    power = 1

                # Normalize the base unit and get any prefix multiplier
                base_unit, prefix_mult = cls.normalize_unit(base_unit)
                multiplier *= prefix_mult**power

                # Add to numerator
                if base_unit in parts["numerator"]:
                    parts["numerator"][base_unit] += power
                else:
                    parts["numerator"][base_unit] = power

        return parts, multiplier

    @classmethod
    def get_unit_dimension(cls, unit):
        # Special cases for complex natural units
        if unit in ["c", "speed_of_light"]:
            return "velocity"

        if unit in ["eV/c^2", "eV/c²", "electronvolt_mass"]:
            return "mass"

        if unit in ["ħc/eV", "ħ*c/eV", "hbar*c/eV", "eV^-1_length"]:
            return "length"

        if unit in ["ħ/eV", "hbar/eV", "eV^-1_time"]:
            return "time"

        # For named units like pascal, newton
        if unit in cls.compound_units:
            return cls.compound_units[unit]["dimension"]

        # Remove any power notation for lookup
        base_unit = unit.split("^")[0].strip()

        # Normalize unit name
        canonical_unit, _ = cls.normalize_unit(base_unit)

        # If it's a compound unit, return its dimension
        if canonical_unit in cls.compound_units:
            return cls.compound_units[canonical_unit]["dimension"]

        # Otherwise check all dimension categories
        for dim, units in cls.dimensions.items():
            if canonical_unit in units:
                return dim

        raise ValueError(
            f"Unknown dimension for unit: {unit} (canonical: {canonical_unit})"
        )

    @classmethod
    def convert_simple_unit(cls, value, from_unit, to_unit):
        # Strip any leading/trailing whitespace
        from_unit = from_unit.strip()
        to_unit = to_unit.strip()

        # Normalize units and get prefix multipliers
        from_canonical, from_multiplier = cls.normalize_unit(from_unit)
        to_canonical, to_multiplier = cls.normalize_unit(to_unit)

        # Get the dimension for both units to ensure compatibility
        from_dim = cls.get_unit_dimension(from_canonical)
        to_dim = cls.get_unit_dimension(to_canonical)

        if from_dim != to_dim:
            raise ValueError(f"Incompatible dimensions: {from_dim} and {to_dim}")

        # Apply prefix multipliers to value
        value = value * from_multiplier / to_multiplier

        # Special case for compound units like jansky
        if from_canonical in cls.compound_units and to_canonical in cls.compound_units:
            # Both are compound units with the same dimension, convert through SI
            from_factor = cls.compound_units[from_canonical]["conversion_factor"]
            to_factor = cls.compound_units[to_canonical]["conversion_factor"]
            return value * from_factor / to_factor

        # Dispatch to the appropriate conversion method
        if from_dim == "energy":
            return cls.convert_energy(value, from_canonical, to_canonical)
        elif from_dim == "length":
            return cls.convert_length(value, from_canonical, to_canonical)
        elif from_dim == "time":
            return cls.convert_time(value, from_canonical, to_canonical)
        elif from_dim == "mass":
            return cls.convert_mass(value, from_canonical, to_canonical)
        elif from_dim == "frequency":
            # Frequency is just 1/time
            return cls.convert_time(value, to_canonical, from_canonical)  # Inversion
        elif from_dim in cls.dimensions:
            # For other dimensions, use compound unit conversion
            if (
                from_canonical in cls.compound_units
                and to_canonical in cls.compound_units
            ):
                # Get the definitions
                from_def = cls.compound_units[from_canonical]["definition"]
                to_def = cls.compound_units[to_canonical]["definition"]
                # Convert through the definitions
                return cls.convert_compound(value, from_def, to_def)
            else:
                raise ValueError(
                    f"Units {from_canonical} and {to_canonical} not directly convertible"
                )
        else:
            raise ValueError(f"Unsupported dimension: {from_dim}")

    @classmethod
    def is_dimensionally_compatible(cls, from_parts, to_parts):
        # DOES NOT WORK FOR NATURAL UNITS - FIX!!!!

        # Check that the same dimensions appear in numerator and denominator
        # from_num_dims = {cls.get_unit_dimension(u) for u in from_parts["numerator"]}
        # to_num_dims = {cls.get_unit_dimension(u) for u in to_parts["numerator"]}

        # from_denom_dims = {cls.get_unit_dimension(u) for u in from_parts["denominator"]}
        # to_denom_dims = {cls.get_unit_dimension(u) for u in to_parts["denominator"]}

        # Check powers for each dimension
        from_dim_powers = {}
        to_dim_powers = {}

        # Count powers in numerator
        for unit, power in from_parts["numerator"].items():
            dim = cls.get_unit_dimension(unit)
            from_dim_powers[dim] = from_dim_powers.get(dim, 0) + power

        for unit, power in to_parts["numerator"].items():
            dim = cls.get_unit_dimension(unit)
            to_dim_powers[dim] = to_dim_powers.get(dim, 0) + power

        # Subtract powers in denominator
        for unit, power in from_parts["denominator"].items():
            dim = cls.get_unit_dimension(unit)
            from_dim_powers[dim] = from_dim_powers.get(dim, 0) - power

        for unit, power in to_parts["denominator"].items():
            dim = cls.get_unit_dimension(unit)
            to_dim_powers[dim] = to_dim_powers.get(dim, 0) - power

        # Compare final dimension powers
        return from_dim_powers == to_dim_powers

    @classmethod
    def convert_compound(cls, value, from_unit_str, to_unit_str):
        # Handle special case of simple/predefined units
        if from_unit_str in cls.compound_units and to_unit_str in cls.compound_units:
            # Both are predefined compound units, check if they're compatible
            from_dim = cls.compound_units[from_unit_str]["dimension"]
            to_dim = cls.compound_units[to_unit_str]["dimension"]

            if from_dim == to_dim:
                # Use simple conversion
                return cls.convert_simple_unit(value, from_unit_str, to_unit_str)

        # Parse the compound units
        from_parts, from_multiplier = cls.parse_compound_unit(from_unit_str)
        to_parts, to_multiplier = cls.parse_compound_unit(to_unit_str)

        # Apply prefix multipliers
        value = value * from_multiplier / to_multiplier

        # Check dimensional compatibility
        if not cls.is_dimensionally_compatible(from_parts, to_parts):
            raise ValueError(
                f"Units {from_unit_str} and {to_unit_str} are not dimensionally compatible"
            )

        # Convert by applying conversions for each component
        result = value

        # Convert numerator components
        for from_unit, from_power in from_parts["numerator"].items():
            # Find matching dimension in target
            from_dim = cls.get_unit_dimension(from_unit)
            to_unit = None

            for candidate in to_parts["numerator"]:
                if cls.get_unit_dimension(candidate) == from_dim:
                    to_unit = candidate
                    break

            if to_unit:
                # Convert and apply the power
                conversion_factor = cls.convert_simple_unit(1.0, from_unit, to_unit)
                result *= conversion_factor**from_power

        # Convert denominator components
        for from_unit, from_power in from_parts["denominator"].items():
            # Find matching dimension in target
            from_dim = cls.get_unit_dimension(from_unit)
            to_unit = None

            for candidate in to_parts["denominator"]:
                if cls.get_unit_dimension(candidate) == from_dim:
                    to_unit = candidate
                    break

            if to_unit:
                # Convert and apply the power (inverse because it's denominator)
                conversion_factor = cls.convert_simple_unit(1.0, from_unit, to_unit)
                result /= conversion_factor**from_power

        return result

    @classmethod
    def handle_natural_units(cls, unit_str):
        """
        Special handling for eV-based units with natural constants

        Args:
            unit_str (str): Unit string like "eV" or "eV/c^2"

        Returns:
            str: The physical dimension this unit represents

        Notes:
            In natural units where ħ=c=1:
            - eV is energy
            - eV/c² is mass (E=mc²): eV -> eV/c² = eV/(3e8)² ≈ 1.78e-36 kg
            - eV⁻¹ as length (ħc/E): ħc/eV ≈ 1.97e-7 m
            - eV⁻¹ as time (ħ/E): ħ/eV ≈ 6.58e-16 s

            Conversion factors are stored in the base dictionaries.
        """
        # Normalize the unit first
        unit, _ = cls.normalize_unit(unit_str)

        if "eV" in unit:
            # Common patterns in high-energy physics
            if unit == "eV":
                # Context-dependent: could be energy or mass (E=mc²)
                return "energy"  # Default to energy
            elif unit in ["eV/c^2", "eV/c²"]:
                return "mass"
            elif unit == "eV^-1":
                # Could be length (ħc/E) or time (ħ/E)
                return "length"  # Default to length

        # Default to normal handling
        return cls.get_unit_dimension(unit)

    @classmethod
    def convert_length(cls, value, from_unit, to_unit):
        # Normalize units
        from_unit = cls.unit_aliases.get(from_unit, from_unit)
        to_unit = cls.unit_aliases.get(to_unit, to_unit)

        base_value = value * cls.lengths[from_unit]
        converted = base_value / cls.lengths[to_unit]
        return converted

    @classmethod
    def convert_mass(cls, value, from_unit, to_unit):
        # Normalize units
        from_unit = cls.unit_aliases.get(from_unit, from_unit)
        to_unit = cls.unit_aliases.get(to_unit, to_unit)

        base_value = value * cls.masses[from_unit]
        converted = base_value / cls.masses[to_unit]
        return converted

    @classmethod
    def convert_time(cls, value, from_unit, to_unit):
        # Normalize units
        from_unit = cls.unit_aliases.get(from_unit, from_unit)
        to_unit = cls.unit_aliases.get(to_unit, to_unit)

        base_value = value * cls.times[from_unit]
        converted = base_value / cls.times[to_unit]
        return converted

    @classmethod
    def convert_energy(cls, value, from_unit, to_unit):
        # Normalize units
        from_unit = cls.unit_aliases.get(from_unit, from_unit)
        to_unit = cls.unit_aliases.get(to_unit, to_unit)

        base_value = value * cls.energies[from_unit]
        converted = base_value / cls.energies[to_unit]
        return converted

    @classmethod
    def convert(cls, value, from_unit, to_unit):
        """
        Main public conversion method that handles all unit types.

        Args:
            value (float): Value to convert
            from_unit (str): Source unit string (simple or compound)
            to_unit (str): Target unit string (simple or compound)

        Returns:
            float: Converted value
        """
        # Check if these are predefined compound units
        if from_unit in cls.compound_units or to_unit in cls.compound_units:
            # Handle predefined compound units
            if from_unit in cls.compound_units and to_unit in cls.compound_units:
                # Both are compound units, check if they have same dimension
                from_dim = cls.compound_units[from_unit]["dimension"]
                to_dim = cls.compound_units[to_unit]["dimension"]

                if from_dim == to_dim:
                    # Direct conversion between compound units of the same dimension
                    from_factor = cls.compound_units[from_unit]["conversion_factor"]
                    to_factor = cls.compound_units[to_unit]["conversion_factor"]
                    return value * from_factor / to_factor
                else:
                    # Need to expand definitions and convert
                    from_def = cls.compound_units[from_unit]["definition"]
                    to_def = cls.compound_units[to_unit]["definition"]
                    return cls.convert_compound(value, from_def, to_def)
            elif from_unit in cls.compound_units:
                # Expand the from_unit definition and convert to to_unit
                from_def = cls.compound_units[from_unit]["definition"]
                from_factor = cls.compound_units[from_unit]["conversion_factor"]
                return cls.convert_compound(value * from_factor, from_def, to_unit)
            else:  # to_unit in cls.compound_units
                # Expand the to_unit definition and convert from from_unit
                to_def = cls.compound_units[to_unit]["definition"]
                to_factor = cls.compound_units[to_unit]["conversion_factor"]
                return cls.convert_compound(value, from_unit, to_def) / to_factor

        # Check if these are simple units or compound units
        if "/" in from_unit or "*" in from_unit or "/" in to_unit or "*" in to_unit:
            return cls.convert_compound(value, from_unit, to_unit)
        else:
            return cls.convert_simple_unit(value, from_unit, to_unit)


# Access natural constants (first element of tuple is the magnitude)
SPEED_OF_LIGHT = constants.c[0]  # c, m/s
PLANCK_CONSTANT = constants.h[0]  # h, J·s
REDUCED_PLANCK = constants.h_bar[0]  # ħ, J·s
BOLTZMANN_CONSTANT = constants.k_b[0]  # k_B, J/K
ELEMENTARY_CHARGE = constants.e[0]  # e, C

unit_converter = UnitConverter()
# Define newton and speed of light first
unit_converter.add_compound_unit(
    name="newton", definition="kg*m/s^2", dimension="force", aliases=["N"]
)

unit_converter.add_compound_unit(
    name="speed_of_light",
    definition=f"{SPEED_OF_LIGHT} m/s",
    dimension="velocity",
    aliases=["c"],
)

# Add natural unit conversions
unit_converter.add_compound_unit(
    name="electronvolt_mass",
    definition="1.78266192e-36 kg",  # eV/c² in kg
    dimension="mass",
    aliases=["eV/c^2", "eV/c²"],
)

unit_converter.add_compound_unit(
    name="electronvolt_length",
    definition="1.97327e-7 m",  # ħc/eV in m
    dimension="length",
    aliases=["eV^-1_length", "hbar*c/eV", "ħc/eV", "ħ*c/eV"],
)

unit_converter.add_compound_unit(
    name="electronvolt_time",
    definition="6.582119e-16 s",  # ħ/eV in s
    dimension="time",
    aliases=["eV^-1_time", "hbar/eV", "ħ/eV"],
)

# Add specialty units after the class definition
unit_converter.add_compound_unit(
    name="jansky",
    definition="1e-26 watt/meter^2/hertz",
    dimension="energy/time/length^2/frequency",
    #    conversion_factor=1e-26,
    aliases=["Jy", "jy", "Jansky", "janskys", "Janskys"],
)

# Magnetic field units
unit_converter.add_compound_unit(
    name="tesla",
    definition="kg/s^2/A",  # kg per second² per ampere
    dimension="magnetic_flux_density",
    aliases=["T"],
)

unit_converter.add_compound_unit(
    name="gauss",
    definition="1e-4 tesla",
    dimension="magnetic_flux_density",
    conversion_factor=1e-4,  # Explicit conversion factor needed
    aliases=["G", "Gauss"],
)

# Astronomical distance units
unit_converter.add_compound_unit(
    name="parsec",
    definition="3.0857e16 meter",
    dimension="length",
    conversion_factor=3.0857e16,
    aliases=["pc", "parsecs"],
)

unit_converter.add_compound_unit(
    name="light_year",
    definition="9.461e15 meter",
    dimension="length",
    conversion_factor=9.461e15,
    aliases=["ly", "lightyear", "light-year"],
)

unit_converter.add_compound_unit(
    name="astronomical_unit",
    definition="1.496e11 meter",
    dimension="length",
    aliases=["au", "AU"],
)

# Nuclear and radiation units
unit_converter.add_compound_unit(
    name="barn", definition="1e-28 meter^2", dimension="area", aliases=["b", "barns"]
)

unit_converter.add_compound_unit(
    name="becquerel",
    definition="1/second",
    dimension="radioactivity",
    conversion_factor=1,
    aliases=["Bq"],
)

unit_converter.add_compound_unit(
    name="curie",
    definition="3.7e10 becquerel",
    dimension="radioactivity",
    conversion_factor=3.7e10,
    aliases=["Ci"],
)


unit_converter.add_compound_unit(
    name="roentgen",
    definition="2.58e-4 coulomb/kg",
    dimension="exposure",
    aliases=["R"],
)

# Temperature units already handled through base conversion, but add Kelvin explicitly
unit_converter.add_compound_unit(
    name="kelvin",
    definition="1 kelvin",  # Base unit for temperature
    dimension="temperature",
    aliases=["K"],
)

# Pressure units
unit_converter.add_compound_unit(
    name="pascal",
    definition="newton/meter^2",  # N/m²
    dimension="pressure",
    aliases=["Pa"],
)

unit_converter.add_compound_unit(
    name="bar",
    definition="1e5 pascal",  # 100,000 Pa
    dimension="pressure",
    conversion_factor=1e5,
    aliases=["bar"],
)

unit_converter.add_compound_unit(
    name="atmosphere",
    definition="101325 pascal",  # Standard atmospheric pressure
    dimension="pressure",
    aliases=["atm"],
)

unit_converter.add_compound_unit(
    name="torr",
    definition="133.322 pascal",  # 1 mmHg
    dimension="pressure",
    aliases=["Torr", "mmHg"],
)

unit_converter.add_compound_unit(
    name="psi",
    definition="6894.76 pascal",  # Pounds per square inch
    dimension="pressure",
    conversion_factor=6894.76,
    aliases=["PSI", "lb/in^2", "pound_per_square_inch"],
)

# Flow rate units - Volume
unit_converter.add_compound_unit(
    name="cubic_meter_per_second",
    definition="meter^3/second",
    dimension="volume_flow_rate",
    aliases=["m³/s", "m^3/s"],
)

unit_converter.add_compound_unit(
    name="liter_per_second",
    definition="0.001 cubic_meter_per_second",
    dimension="volume_flow_rate",
    aliases=["L/s", "l/s"],
)

unit_converter.add_compound_unit(
    name="cubic_foot_per_second",
    definition="0.0283168 cubic_meter_per_second",
    dimension="volume_flow_rate",
    aliases=["ft³/s", "cfs"],
)

unit_converter.add_compound_unit(
    name="gallon_per_minute",
    definition="6.30902e-5 cubic_meter_per_second",  # US gallon
    dimension="volume_flow_rate",
    aliases=["gpm", "GPM"],
)

# Flow rate units - Mass
unit_converter.add_compound_unit(
    name="kilogram_per_second",
    definition="kg/s",
    dimension="mass_flow_rate",
    aliases=["kg/s"],
)

unit_converter.add_compound_unit(
    name="pound_per_second",
    definition="0.453592 kilogram_per_second",
    dimension="mass_flow_rate",
    aliases=["lb/s"],
)

# Additional energy rate units
unit_converter.add_compound_unit(
    name="btu_per_hour",
    definition="0.293071 watt",  # BTU/h
    dimension="energy/time",
    aliases=["BTU/h", "BTUH"],
)

unit_converter.add_compound_unit(
    name="ton_of_refrigeration",
    definition="3516.85 watt",  # 12,000 BTU/h
    dimension="energy/time",
    aliases=["TR", "ton"],
)
