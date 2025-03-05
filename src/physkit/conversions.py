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
from collections import defaultdict
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
    def __init__(self):
        self.base_dimensions = {"length", "mass", "time", "current", "temperature"}
        self.base_dimension_map = {}  # {"length": "meter", "mass": "gram"...}
        self.si_prefixes = self._create_si_prefixes()
        self.unit_registry = {}
        self.aliases = {}
        self._initialize_base_units()
        self._initialize_compound_units()

    def _create_si_prefixes(self):
        prefixes = {
            # SI prefixes (both symbols and names)
            "y": 1e-24,
            "yocto": 1e-24,
            "z": 1e-21,
            "zepto": 1e-21,
            "a": 1e-18,
            "atto": 1e-18,
            "f": 1e-15,
            "femto": 1e-15,
            "p": 1e-12,
            "pico": 1e-12,
            "n": 1e-9,
            "nano": 1e-9,
            "μ": 1e-6,
            "u": 1e-6,
            "micro": 1e-6,
            "m": 1e-3,
            "milli": 1e-3,
            "c": 1e-2,
            "centi": 1e-2,
            "d": 1e-1,
            "deci": 1e-1,
            "da": 1e1,
            "deca": 1e1,
            "h": 1e2,
            "hecto": 1e2,
            "k": 1e3,
            "kilo": 1e3,
            "M": 1e6,
            "mega": 1e6,
            "G": 1e9,
            "giga": 1e9,
            "T": 1e12,
            "tera": 1e12,
            "P": 1e15,
            "peta": 1e15,
            "E": 1e18,
            "exa": 1e18,
            "Z": 1e21,
            "zetta": 1e21,
            "Y": 1e24,
            "yotta": 1e24,
        }
        self._sorted_prefixes = sorted(prefixes.keys(), key=lambda x: -len(x))
        return prefixes

    def _initialize_base_units(self):
        # Base units
        self._add_base_unit("meter", ["length"], 1.0, ["m", "metre", "metres"])
        self._add_base_unit("gram", ["mass"], 1.0, ["g"])
        self._add_base_unit("second", ["time"], 1.0, ["s", "sec"])
        self._add_base_unit("ampere", ["current"], 1.0, ["A", "amp"])
        self._add_base_unit("kelvin", ["temperature"], 1.0, ["K"])

    def _initialize_compound_units(self):
        # Predefined compound units
        self.add_unit("inch", "0.0254 meter", aliases=["in"])
        self.add_unit("mile", "1609.34 meter", aliases=["mi"])
        self.add_unit("foot", "0.3048 meter", aliases=["ft", "feet"])
        self.add_unit("pound", "453.592 g", aliases=["lb", "lbs", "Lb"])
        self.add_unit("hour", "3600 second", aliases=["hr"])

        # Add physics units
        self.add_unit("newton", "kg*m/s^2", aliases=["N"])
        self.add_unit("joule", "newton*meter", aliases=["J"])
        self.add_unit("watt", "joule/second", aliases=["W"])
        self.add_unit("coulomb", "ampere*second", aliases=["C"])
        self.add_unit("hertz", "second^-1", aliases=["Hz"])
        self.add_unit("tesla", "kg/s^2/A", aliases=["T"])
        self.add_unit("sievert", "m^2 s^-2", aliases=["Sv"])
        self.add_unit("abampere", "10 ampere", aliases=["abA"])
        self.add_unit("statampere", "3.336e-10 ampere", aliases=["statA"])
        self.add_unit("biot", "10 ampere", aliases=["Bi"])

        # Astronomy units
        self.add_unit("parsec", "3.0857e16 meter", aliases=["pc"])
        self.add_unit("light_year", "9.461e15 meter", aliases=["ly"])
        self.add_unit("astronomical_unit", "1.496e11 meter", aliases=["au"])

        # Radiation units
        self.add_unit("becquerel", "1/second", aliases=["Bq"])
        self.add_unit("curie", "3.7e10 becquerel", aliases=["Ci"])
        self.add_unit("gray", "joule/kilogram", aliases=["Gy"])
        self.add_unit("rad", "0.01 gray", aliases=["rd"])

        # Pressure units
        self.add_unit("pascal", "newton/meter^2", aliases=["Pa"])
        self.add_unit("bar", "1e5 pascal", aliases=["bar"])
        self.add_unit("atmosphere", "101325 pascal", aliases=["atm"])
        self.add_unit("torr", "133.322 pascal", aliases=["mmHg"])
        self.add_unit("psi", "6894.76 pascal", aliases=["lb/in^2"])

        # Radiation units
        self.add_unit("roentgen", "2.58e-4 coulomb/kg", aliases=["R"])

        # Flow rates
        self.add_unit("cubic_meter_per_second", "m^3/s", aliases=["m³/s"])
        self.add_unit(
            "liter_per_second", "0.001 cubic_meter_per_second", aliases=["L/s"]
        )
        self.add_unit(
            "gallon_per_minute", "6.30902e-5 cubic_meter_per_second", aliases=["gpm"]
        )

        # Energy rates
        self.add_unit("btu_per_hour", "0.293071 watt", aliases=["BTU/h"])
        self.add_unit("ton_of_refrigeration", "3516.85 watt", aliases=["TR"])

        # Special units
        self.add_unit("jansky", "1e-26 watt/meter^2/hertz)", aliases=["Jy"])
        self.add_unit("gauss", "1e-4 tesla", aliases=["G"])
        self.add_unit("calorie", "4.184 joule", aliases=["cal"])

        # Nuclear
        self.add_unit("barn", "1e-28 meter^2", aliases=["b"])
        self.add_unit("rem", "0.01 sievert", aliases=["rem"])  # Dose equivalent

        # Natural units
        self.add_natural_units()

    def _add_base_unit(self, name, dimensions, factor, aliases=None):
        if len(dimensions) != 1:
            raise ValueError("Base units must have exactly one dimension")
        dimension = self._create_dimension_dict(dimensions)

        # base_dim = dimension.split("^")[0] if "^" in dimension else dimension
        # self.base_dimension_map[base_dim] = name
        for dim in dimensions:
            if "^" in dim:
                base_dim = dim.split("^")[0]  # Handle "length^2" cases
            else:
                base_dim = dim
            self.base_dimension_map[base_dim] = name

        self.unit_registry[name] = {
            "factor": factor,
            "dimension": dimension,
            "is_base": True,
        }
        if aliases:
            for alias in aliases:
                self.aliases[alias] = name

    def add_natural_units(self):
        """Add units related to natural constants and physics"""
        # Natural constants (values from scipy.constants)
        c = 299792458.0  # Speed of light (m/s)
        h_bar = 1.054571817e-34  # Reduced Planck constant (J·s)
        eV = 1.602176634e-19  # 1 eV in joules

        # Speed of light
        self.add_unit("speed_of_light", f"{c} m/s", aliases=["c"])

        # Natural units system
        self.add_unit("electronvolt", f"{eV} joule", aliases=["eV"])
        self.add_unit("electronvolt_mass", f"{eV/c**2} kg", aliases=["eV/c^2", "eV/c²"])

        self.add_unit(
            "electronvolt_length", f"{h_bar*c/eV} m", aliases=["eV^-1_length", "ħc/eV"]
        )

        self.add_unit(
            "electronvolt_time", f"{h_bar/eV} s", aliases=["eV^-1_time", "ħ/eV"]
        )

    def add_unit(self, name, definition, aliases=None):
        if name in self.unit_registry:
            raise ValueError(f"Unit {name} already exists")
        if any(a in self.aliases for a in aliases or []):
            conflicts = [a for a in aliases if a in self.aliases]
            raise ValueError(f"Aliases {conflicts} already registered")

        factor, components = self.parse_definition(definition)

        self.unit_registry[name] = {
            "factor": factor,
            "dimension": components,
            "is_base": False,
        }

        if aliases:
            for alias in aliases:
                self.aliases[alias] = name

    def parse_definition(self, definition):
        # Convert implicit multiplication (spaces) to explicit *
        definition = self._resolve_aliases_in_definition(definition)
        definition = definition.replace(" ", "*")

        factor = 1.0
        components = defaultdict(float)

        # Split into tokens by * and /
        tokens = re.split(
            r"(\*|/)", definition.replace("(", "").replace(")", "").strip()
        )

        current_operation = "*"  # initial operation is multiplication

        i = 0
        while i < len(tokens):
            token = tokens[i].strip()
            if not token:
                i += 1
                continue

            if token in ("*", "/"):
                current_operation = token
                i += 1
                continue

            # Handle each term (could have exponent)
            term = token

            # Split into unit and exponent part
            if "^" in term:
                unit_part, exp_part = term.split("^", 1)
                try:
                    exp = float(exp_part)
                except ValueError:
                    raise ValueError(f"Invalid exponent in term '{term}'")
            else:
                unit_part = term
                exp = 1.0

            # Apply current operation to exponent
            if current_operation == "/":
                exp = -exp

            # Reset current operation to * for next term
            current_operation = "*"

            # Extract numerical coefficient from unit_part
            coef = 1.0
            coef_match = re.match(r"^([+-]?\d+\.?\d*([eE][+-]?\d+)?)", unit_part)
            if coef_match:
                coef_str = coef_match.group(1)
                coef = float(coef_str)
                unit_part = unit_part[len(coef_str) :].strip()

            # Handle pure numerical terms (no unit)
            if not unit_part:
                factor *= coef**exp
                i += 1
                continue

            # Split prefix and base unit
            prefix, base = self.split_prefix(unit_part)

            # Check if base unit exists
            if not self.is_valid_unit(base):
                raise ValueError(f"Unknown unit '{base}' in term '{term}'")

            # Get factors
            prefix_factor = self.si_prefixes.get(prefix, 1.0)
            base_factor = self.get_base_unit_factor(base)
            total_unit_factor = prefix_factor * base_factor

            # Update overall factor
            factor *= coef * (total_unit_factor**exp)

            # Update dimensions
            unit_dim = self.get_unit_dimension(base)
            for dimension, power in unit_dim.items():
                components[dimension] += power * exp

            i += 1

        # Remove components that have 0 power
        components = self._clean_dimensions(components)
        return factor, dict(components)

    def _resolve_aliases_in_definition(self, definition):
        # Replace aliases in reverse-length order to handle compound units first
        sorted_aliases = sorted(self.aliases.keys(), key=lambda x: (-len(x), x))

        # Use word boundaries to prevent partial matches
        for alias in sorted_aliases:
            pattern = r"\b" + re.escape(alias) + r"\b"
            definition = re.sub(pattern, self.aliases[alias], definition)

        return definition

    def split_prefix(self, unit):
        # First check if the entire unit is valid
        if self.is_valid_unit(unit):
            return "", unit

        # Then check for prefixes
        for prefix in self._sorted_prefixes:
            if unit.startswith(prefix):
                remaining = unit[len(prefix) :]
                if self.is_valid_unit(remaining):
                    return prefix, remaining
        return "", unit

    def is_valid_unit(self, unit):
        # Check both direct registry and aliases
        return unit in self.unit_registry or unit in self.aliases

    def get_base_unit_factor(self, unit):
        # Follow aliases to canonical name
        canonical = self.aliases.get(unit, unit)
        if canonical not in self.unit_registry:
            raise ValueError(f"Unknown unit: {unit}")
        return self.unit_registry.get(canonical, {}).get("factor", 1.0)

    def get_unit_dimension(self, unit):
        canonical = self.aliases.get(unit, unit)
        if canonical not in self.unit_registry:
            raise ValueError(f"Unknown unit: {unit}")
        return self.unit_registry[canonical]["dimension"]

    def _create_dimension_dict(self, dimensions):
        dim_dict = defaultdict(float)
        for dim in dimensions:
            if "^" in dim:
                base, exp = dim.split("^")
                dim_dict[base] += float(exp)
            else:
                dim_dict[dim] += 1
        return dict(dim_dict)

    def _clean_dimensions(self, dimension_dict):
        """Remove dimensions with zero exponents"""
        return {k: v for k, v in dimension_dict.items() if v != 0}

    def convert(self, value, from_unit, to_unit):
        """
        Convert between compatible units with dimensional analysis

        Args:
            value: Numerical value to convert
            from_unit: Source unit (supports aliases and prefixes)
            to_unit: Target unit (supports aliases and prefixes)

        Returns:
            Converted value in target units

        Raises:
            ValueError: On dimensional mismatch or unknown units
        """
        from_factor, from_dim = self._get_unit_info(from_unit)
        to_factor, to_dim = self._get_unit_info(to_unit)

        # Clean dimensions (redundant safeguard)
        from_dim = self._clean_dimensions(from_dim)
        to_dim = self._clean_dimensions(to_dim)

        if from_dim != to_dim:
            raise ValueError(f"Incompatible dimensions: {from_dim} vs {to_dim}")

        return value * from_factor / to_factor

    def _get_unit_info(self, unit):
        # Get both the parsed factor and components
        parsed_factor, components = self.parse_definition(unit)

        # Start with the parsed factor (numerical coefficients + prefixes)
        total_factor = parsed_factor

        # Multiply by base unit conversion factors
        for dim, exp in components.items():
            base_unit = self._dim_to_base(dim)
            base_factor = self.unit_registry[base_unit]["factor"]
            total_factor *= base_factor**exp

        return total_factor, components

    def _dim_to_base(self, dim):
        try:
            return self.base_dimension_map[dim]
        except KeyError:
            raise ValueError(f"No base unit found for dimension {dim}")


# Access natural constants (first element of tuple is the magnitude)
SPEED_OF_LIGHT = constants.c[0]  # c, m/s
PLANCK_CONSTANT = constants.h[0]  # h, J·s
REDUCED_PLANCK = constants.h_bar[0]  # ħ, J·s
BOLTZMANN_CONSTANT = constants.k_b[0]  # k_B, J/K
ELEMENTARY_CHARGE = constants.e[0]  # e, C

# unit_converter = UnitConverter()
# No need to define joule here - already defined in class-level compound_units
#
## Define newton
# unit_converter.add_compound_unit(
#    name="newton", definition="kg*m/s^2", dimension="force", aliases=["N"]
# )
#
# unit_converter.add_compound_unit(
#    name="speed_of_light",
#    definition=f"{SPEED_OF_LIGHT} m/s",
#    dimension="velocity",
#    aliases=["c"],
# )
#
## Add natural unit conversions
# unit_converter.add_compound_unit(
#    name="electronvolt_mass",
#    definition="1.78266192e-36 kg",  # eV/c² in kg
#    dimension="mass",
#    aliases=["eV/c^2", "eV/c²"],
# )
#
# unit_converter.add_compound_unit(
#    name="electronvolt_length",
#    definition="1.97327e-7 m",  # ħc/eV in m
#    dimension="length",
#    aliases=["eV^-1_length", "hbar*c/eV", "ħc/eV", "ħ*c/eV"],
# )
#
# unit_converter.add_compound_unit(
#    name="electronvolt_time",
#    definition="6.582119e-16 s",  # ħ/eV in s
#    dimension="time",
#    aliases=["eV^-1_time", "hbar/eV", "ħ/eV"],
# )
#
## Add specialty units after the class definition
# unit_converter.add_compound_unit(
#    name="jansky",
#    definition="1e-26 watt/meter^2/hertz",
#    dimension="mass/time^2",  # Simplified from mass*length^2/time^3 * 1/length^2 * time = mass/time^2
#    aliases=["Jy", "jy", "Jansky", "janskys", "Janskys"],
# )
#
## Define gauss as a derived unit of tesla
# unit_converter.add_compound_unit(
#    name="gauss",
#    definition="1e-4 tesla",  # 1 gauss = 1e-4 tesla
#    # No dimension provided - should be inferred from tesla
#    aliases=["G", "Gauss"],
# )
#
## Astronomical distance units
# unit_converter.add_compound_unit(
#    name="parsec",
#    definition="3.0857e16 meter",
#    dimension="length",
#    aliases=["pc", "parsecs"],
# )
#
# unit_converter.add_compound_unit(
#    name="light_year",
#    definition="9.461e15 meter",
#    dimension="length",
#    aliases=["ly", "lightyear", "light-year"],
# )
#
# unit_converter.add_compound_unit(
#    name="astronomical_unit",
#    definition="1.496e11 meter",
#    dimension="length",
#    aliases=["au", "AU"],
# )
#
## Nuclear and radiation units
# unit_converter.add_compound_unit(
#    name="barn",
#    definition="1e-28 meter^2",
#    dimension="length^2",
#    aliases=["b", "barns"],
# )
#
# unit_converter.add_compound_unit(
#    name="becquerel",
#    definition="1/second",
#    dimension="1/time",  # Radioactivity is decays per time unit
#    aliases=["Bq"],
# )
#
# unit_converter.add_compound_unit(
#    name="curie",
#    definition="3.7e10 becquerel",
#    dimension="1/time",  # Same dimension as becquerel
#    aliases=["Ci"],
# )
#
# unit_converter.add_compound_unit(
#    name="roentgen",
#    definition="2.58e-4 coulomb/kg",
#    dimension="current*time/mass",  # Charge/mass = current*time/mass
#    aliases=["R"],
# )
#
## Temperature units already handled through base conversion, but add Kelvin explicitly
# unit_converter.add_compound_unit(
#    name="kelvin",
#    definition="1 kelvin",  # Base unit for temperature
#    dimension="temperature",
#    aliases=["K"],
# )
#
## Pressure units
# unit_converter.add_compound_unit(
#    name="pascal",
#    definition="newton/meter^2",  # N/m²
#    dimension="mass/length/time^2",  # Force/area = mass*acceleration/area
#    aliases=["Pa"],
# )
#
# unit_converter.add_compound_unit(
#    name="bar",
#    definition="1e5 pascal",  # 100,000 Pa
#    dimension="mass/length/time^2",  # Same as pascal
#    aliases=["bar"],
# )
#
# unit_converter.add_compound_unit(
#    name="atmosphere",
#    definition="101325 pascal",  # Standard atmospheric pressure
#    dimension="mass/length/time^2",  # Same as pascal
#    aliases=["atm"],
# )
#
# unit_converter.add_compound_unit(
#    name="torr",
#    definition="133.322 pascal",  # 1 mmHg
#    dimension="mass/length/time^2",  # Same as pascal
#    aliases=["Torr", "mmHg"],
# )
#
# unit_converter.add_compound_unit(
#    name="psi",
#    definition="6894.76 pascal",  # Pounds per square inch
#    dimension="mass/length/time^2",  # Same as pascal
#    aliases=["PSI", "lb/in^2", "pound_per_square_inch"],
# )
#
## Flow rate units - Volume
# unit_converter.add_compound_unit(
#    name="cubic_meter_per_second",
#    definition="meter^3/second",
#    dimension="length^3/time",  # Volume per time
#    aliases=["m³/s", "m^3/s"],
# )
#
# unit_converter.add_compound_unit(
#    name="liter_per_second",
#    definition="0.001 cubic_meter_per_second",
#    dimension="length^3/time",  # Same as cubic_meter_per_second
#    aliases=["L/s", "l/s"],
# )
#
# unit_converter.add_compound_unit(
#    name="cubic_foot_per_second",
#    definition="0.0283168 cubic_meter_per_second",
#    dimension="length^3/time",  # Same as cubic_meter_per_second
#    aliases=["ft³/s", "cfs"],
# )
#
# unit_converter.add_compound_unit(
#    name="gallon_per_minute",
#    definition="6.30902e-5 cubic_meter_per_second",  # US gallon
#    dimension="length^3/time",  # Same as cubic_meter_per_second
#    aliases=["gpm", "GPM"],
# )
#
## Flow rate units - Mass
# unit_converter.add_compound_unit(
#    name="kilogram_per_second",
#    definition="kg/s",
#    dimension="mass/time",  # Mass per time
#    aliases=["kg/s"],
# )
#
# unit_converter.add_compound_unit(
#    name="pound_per_second",
#    definition="0.453592 kilogram_per_second",
#    dimension="mass/time",  # Same as kilogram_per_second
#    aliases=["lb/s"],
# )
#
## Additional energy rate units
# unit_converter.add_compound_unit(
#    name="btu_per_hour",
#    definition="0.293071 watt",  # BTU/h
#    dimension="energy/time",  # Energy per time
#    aliases=["BTU/h", "BTUH"],
# )
#
# unit_converter.add_compound_unit(
#    name="ton_of_refrigeration",
#    definition="3516.85 watt",  # 12,000 BTU/h
#    dimension="energy/time",  # Same as watt
#    aliases=["TR", "ton"],
# )
#
