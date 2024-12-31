# tests/test_equations.py

from physkit.config import set_default
from physkit.equations import equations
import math


def test_schwartzschild_radius():
    set_default("cgs")
    schwartz_radius, _ = equations.schwartzschild_radius(500, "M_sun", "km")
    assert math.isclose(schwartz_radius, 1477, rel_tol=1e-2)


def test_eddington_luminosity():
    set_default("cgs")
    L_edd, _ = equations.eddington_luminosity(500, "M_sun", "erg / s")
    assert math.isclose(L_edd, 6.28750e40, rel_tol=1e-3)
