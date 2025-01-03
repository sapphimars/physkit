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
physkit

A Python toolkit providing physical constants, unit conversions,
and plotting styles for scientific work.
"""

__version__ = "0.1.0"

from .config import get_ureg, set_default
from .constants import constants
from .conversions import convert_unit, convert_to_base_units, convert_unit_system
from .plot_styler import plot_styler
from .equations import equations


__version__ = "0.1.21"

__all__ = [
    "set_default",
    "get_ureg",
    "constants",
    "convert_unit",
    "convert_to_base_units",
    "convert_unit_system",
    "plot_styler",
    "equations",
]
