# tests/test_conversions.py

import math
from physkit.config import set_default
from physkit.conversions import convert_unit, convert_to_base_units


def test_convert_unit_basic():
    set_default("SI")
    # 1 meter = 100 cm
    cm = convert_unit(1, "meter", "cm")
    assert math.isclose(cm, 100.0, rel_tol=1e-9)

    # 1 meter ≈ 39.3701 inches
    inch = convert_unit(1, "meter", "inch")
    assert math.isclose(inch, 39.3700787, rel_tol=1e-6)


def test_convert_unit_cgs():
    set_default("cgs")
    # 1 cm -> base units in cgs is still "cm"
    val, _ = convert_to_base_units(1, "cm")
    assert math.isclose(val, 1.0, rel_tol=1e-9)


def test_convert_energy():
    set_default("SI")
    # 1 erg = 1e-7 J, so 5 erg = 5e-7 J
    joules = convert_unit(5, "erg", "J")
    assert math.isclose(joules, 5e-7, rel_tol=1e-12)


def test_convert_complex():
    gps = 7.1e-17  # g/s
    gps_converted = convert_unit(gps, "g / s", "M_sun / year")

    assert math.isclose(gps_converted, 1.12679e-42, rel_tol=1e-6)
