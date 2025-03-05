import pytest
from physkit.conversions import UnitConverter


@pytest.fixture
def converter():
    """Fixture providing a configured UnitConverter"""
    uc = UnitConverter()

    # Add custom units needed for tests
    uc.add_unit(
        "weber",
        "tesla*meter^2",
        aliases=["Wb"],
    )

    uc.add_unit(
        "maxwell",
        "1e-8 weber",
        aliases=["Mx"],
    )

    uc.add_unit(
        "electronvolt_per_second",
        "eV/second",
        aliases=["eV/s"],
    )

    uc.add_unit("liter", "0.001 meter^3", aliases=["L", "l"])

    return uc


class TestUnitConverter:
    def test_basic_conversions(self, converter):
        """Test basic unit conversions between different systems"""
        # Length
        assert round(converter.convert(1.0, "meter", "inch"), 2) == 39.37
        # Mass
        assert round(converter.convert(1.0, "kg", "pound"), 2) == 2.20
        # Time
        assert converter.convert(1.0, "hour", "second") == 3600.0
        # Energy
        assert (
            pytest.approx(converter.convert(1.0, "joule", "eV"), rel=1e-3) == 6.242e18
        )

    def test_prefix_handling(self, converter):
        """Test SI prefix handling with both symbols and full words"""
        # Symbol prefixes
        assert converter.convert(1.0, "km", "meter") == 1000.0
        assert converter.convert(1.0, "cm", "meter") == 0.01
        assert converter.convert(1.0, "μs", "second") == 1e-6

        # Full word prefixes
        assert converter.convert(1.0, "kilometer", "meter") == 1000.0
        assert converter.convert(1.0, "millimeter", "meter") == 0.001
        assert converter.convert(1.0, "megajoule", "joule") == 1e6

    def test_compound_units(self, converter):
        """Test compound unit conversions"""
        # Simple compound
        assert converter.convert(1.0, "joule/second", "watt") == 1.0

        # Different units in compound
        assert (
            pytest.approx(converter.convert(1.0, "meter/second", "km/hour"), rel=1e-5)
            == 3.6
        )

        # With powers
        assert converter.convert(5.0, "newton/meter^2", "pascal") == 5.0

        result = converter.convert(1.0e-26, "watt/meter^2/hertz", "jansky")
        assert pytest.approx(result, rel=1e-5) == 1.0

    def test_natural_units(self, converter):
        """Test natural unit conversions for high energy physics"""
        # Energy
        assert (
            pytest.approx(converter.convert(1.0, "eV", "joule"), rel=1e-5)
            == 1.602176634e-19
        )

        # Mass - special case for eV/c^2
        assert (
            pytest.approx(converter.convert(1.0, "eV/c^2", "kg"), rel=1e-5)
            == 1.78266192e-36
        )

        assert (
            pytest.approx(converter.convert(1.0, "ħc/eV", "meter"), rel=1e-5)
            == 1.97327e-7
        )
        assert (
            pytest.approx(converter.convert(1.0, "ħ/eV", "second"), rel=1e-5)
            == 6.582119e-16
        )

        assert (
            pytest.approx(converter.nat_convert(1.78266192e-27, "kg"), rel=1e-5)
        ) == 1e9

        assert (
            pytest.approx(converter.nat_convert(5.344286e-19, "kg*m/s"), rel=1e-5)
            == 1e9
        )

        assert (pytest.approx(converter.nat_convert(11604.51812, "K"), rel=1e-5)) == 1.0

    def test_specialty_units(self, converter):
        """Test specialty unit conversions"""
        # Magnetic field
        assert (
            pytest.approx(converter.convert(1.0, "tesla", "gauss"), rel=1e-3) == 10000.0
        )

        assert (
            pytest.approx(converter.convert(3.26, "light_year", "parsec"), rel=1e-2)
            == 1.0
        )
        assert pytest.approx(converter.convert(14.5, "psi", "bar"), rel=1e-2) == 1.0
        assert converter.convert(1.0, "curie", "becquerel") == 3.7e10

    def test_current_units(self, converter):
        """Test electric current unit conversions"""
        # Test SI to CGS conversion
        assert (
            pytest.approx(converter.convert(1.0, "ampere", "abampere"), rel=1e-6) == 0.1
        )
        assert pytest.approx(converter.convert(1.0, "ampere", "biot"), rel=1e-6) == 0.1

        # Test SI to ESU conversion
        statamps_per_amp = 1.0 / 3.336e-10
        assert (
            pytest.approx(converter.convert(1.0, "ampere", "statampere"), rel=1e-6)
            == statamps_per_amp
        )

    def test_dimension_compatibility(self, converter):
        """Test dimension compatibility checks"""
        with pytest.raises(ValueError):
            converter.convert(1.0, "meter", "second")

    def test_ev_aliases(self, converter):
        """Test that eV aliases work properly for all dimensions"""
        assert (
            pytest.approx(converter.convert(1.0, "electronvolt", "joule"), rel=1e-5)
            == 1.602176634e-19
        )
        assert (
            pytest.approx(converter.convert(1.0, "eV", "joule"), rel=1e-5)
            == 1.602176634e-19
        )

        assert (
            pytest.approx(converter.convert(1.0, "eV/c^2", "kg"), rel=1e-5)
            == 1.78266192e-36
        )

    def test_custom_unit_definition(self, converter):
        """Test defining and using new custom units"""
        assert (
            pytest.approx(converter.convert(1.0, "weber", "maxwell"), rel=1e-5) == 1e8
        )
        assert (
            pytest.approx(converter.convert(2.5, "tesla*meter^2", "weber"), rel=1e-5)
            == 2.5
        )
        assert (
            pytest.approx(converter.convert(5.0, "weber", "tesla*meter^2"), rel=1e-5)
            == 5.0
        )

    def test_natural_unit_compound_conversion(self, converter):
        """Test conversions involving natural units in compound expressions"""
        assert (
            pytest.approx(converter.convert(1.0, "eV/s", "joule/second"), rel=1e-5)
            == 1.602176634e-19
        )
        assert (
            pytest.approx(converter.convert(1.0, "joule/second", "eV/s"), rel=1e-5)
            == 6.241509e18
        )

    def test_dimension_analysis_and_compatibility(self, converter):
        """Test dimension analysis for complex unit expressions"""
        assert (
            pytest.approx(converter.convert(1.0, "pascal", "psi"), rel=1e-2) == 0.000145
        )
        assert (
            pytest.approx(converter.convert(1.0, "meter^3", "liter"), rel=1e-2)
            == 1000.0
        )

    def test_add_new_unit_with_extended_aliases(self, converter):
        """Test adding new units with extended aliases"""
        assert pytest.approx(converter.convert(1.0, "gray", "rad"), rel=1e-5) == 100.0
        assert (
            pytest.approx(converter.convert(5.0, "gray", "joule/kg"), rel=1e-5) == 5.0
        )

    def test_base_dimension_conversions(self, converter):
        """Test that units defined with base physical dimensions work correctly"""

        assert (
            pytest.approx(converter.convert(2.0, "hertz", "1/second"), rel=1e-5) == 2.0
        )
        assert (
            pytest.approx(converter.convert(101325, "pascal", "atmosphere"), rel=1e-5)
            == 1.0
        )
