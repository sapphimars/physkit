# tests/test_constants.py

import math
from physkit.config import set_default
from physkit.constants import constants


def test_speed_of_light_si():
    # Ensure we're in SI
    set_default("SI")
    # Suppose your .c property returns (float_value, "m / s") in SI
    val, unit = constants.c

    # Check magnitude
    assert math.isclose(val, 2.99792458e8, rel_tol=1e-7)
    # Check unit string
    assert unit == "m / s"


def test_speed_of_light_cgs():
    # Switch to cgs
    set_default("cgs")
    val, unit = constants.c

    # c ~ 3.0e10 in cm/s
    assert math.isclose(val, 2.99792458e10, rel_tol=1e-7)
    assert unit == "cm / s"


def test_fine_structure():
    # alpha might be dimensionless
    # Suppose .alpha returns (0.007297..., "")
    set_default("SI")  # or cgs, shouldn't matter
    val, unit = constants.alpha
    assert math.isclose(val, 0.0072973525693, rel_tol=1e-9)
    assert unit == ""
