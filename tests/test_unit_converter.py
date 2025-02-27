"""
Test module for the new UnitConverter class
"""

import pytest
from physkit.conversions import UnitConverter


class TestUnitConverter:
    def test_basic_conversions(self):
        """Test basic unit conversions between different systems"""
        # Length
        assert round(UnitConverter.convert(1.0, "meter", "inch"), 2) == 39.37
        # Mass
        assert round(UnitConverter.convert(1.0, "kg", "pound"), 2) == 2.20
        # Time
        assert UnitConverter.convert(1.0, "hour", "second") == 3600.0
        # Energy
        assert (
            pytest.approx(UnitConverter.convert(1.0, "joule", "eV"), rel=1e-3)
            == 6.242e18
        )

    def test_prefix_handling(self):
        """Test SI prefix handling with both symbols and full words"""
        # Symbol prefixes
        assert UnitConverter.convert(1.0, "km", "meter") == 1000.0
        assert UnitConverter.convert(1.0, "cm", "meter") == 0.01
        assert UnitConverter.convert(1.0, "μs", "second") == 1e-6

        # Full word prefixes
        assert UnitConverter.convert(1.0, "kilometer", "meter") == 1000.0
        assert UnitConverter.convert(1.0, "millimeter", "meter") == 0.001
        assert UnitConverter.convert(1.0, "megajoule", "joule") == 1e6

    def test_compound_units(self):
        """Test compound unit conversions"""
        # Simple compound
        assert UnitConverter.convert(1.0, "joule/second", "watt") == 1.0

        # Different units in compound
        assert (
            pytest.approx(
                UnitConverter.convert(1.0, "meter/second", "km/hour"), rel=1e-5
            )
            == 3.6
        )

        # With powers
        assert UnitConverter.convert(5.0, "newton/meter^2", "pascal") == 5.0

        result = UnitConverter.convert(1.0e-26, "watt/meter^2/hertz", "jansky")
        assert pytest.approx(result, rel=1e-5) == 1.0

    def test_natural_units(self):
        """Test natural unit conversions for high energy physics"""
        # Energy
        assert (
            pytest.approx(UnitConverter.convert(1.0, "eV", "joule"), rel=1e-5)
            == 1.602176634e-19
        )

        # Mass - special case for eV/c^2
        assert (
            pytest.approx(
                UnitConverter.convert(1.0, "electronvolt_mass", "kg"), rel=1e-5
            )
            == 1.78266192e-36
        )

        assert (
            pytest.approx(
                UnitConverter.convert(1.0, "electronvolt_length", "meter"), rel=1e-5
            )
            == 1.97327e-7
        )
        assert (
            pytest.approx(
                UnitConverter.convert(1.0, "electronvolt_time", "second"), rel=1e-5
            )
            == 6.582119e-16
        )

    def test_specialty_units(self):
        """Test specialty unit conversions"""
        # Magnetic field
        assert (
            pytest.approx(UnitConverter.convert(1.0, "tesla", "gauss"), rel=1e-5)
            == 10000.0
        )

        assert (
            pytest.approx(UnitConverter.convert(3.26, "light_year", "parsec"), rel=1e-2)
            == 1.0
        )

        assert pytest.approx(UnitConverter.convert(14.5, "psi", "bar"), rel=1e-2) == 1.0

        assert UnitConverter.convert(1.0, "curie", "becquerel") == 3.7e10

    def test_dimension_compatibility(self):
        """Test dimension compatibility checks"""
        # Should raise ValueError for incompatible dimensions
        with pytest.raises(ValueError):
            UnitConverter.convert(1.0, "meter", "second")
