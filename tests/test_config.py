import pytest
from rheoscale.config import RheoscaleConfig, FixedKeysDict


class TestRheoscaleConfigCreation:
    """Test that valid configs can be created and invalid ones raise errors."""

    def test_minimal_config(self):
        config = RheoscaleConfig(protein_name="TestProtein")
        assert config.protein_name == "TestProtein"
        assert config.log_scale is False
        assert config.dead_extremum == "Min"

    def test_protein_name_must_be_string(self):
        with pytest.raises(ValueError, match="MUST be a string"):
            RheoscaleConfig(protein_name=123)

    def test_bins_too_low(self):
        with pytest.raises(ValueError, match="at least 2 bins"):
            RheoscaleConfig(protein_name="Test", number_of_bins=1)

    def test_bins_too_high(self):
        with pytest.raises(ValueError, match="no more than 20"):
            RheoscaleConfig(protein_name="Test", number_of_bins=21)

    def test_valid_bin_range(self):
        config = RheoscaleConfig(protein_name="Test", number_of_bins=10)
        assert config.number_of_bins == 10

    def test_number_of_positions_must_be_positive(self):
        with pytest.raises(ValueError, match="at least one position"):
            RheoscaleConfig(protein_name="Test", number_of_positions=0)

    def test_wt_val_and_error_together(self):
        config = RheoscaleConfig(protein_name="Test", WT_val=1.0, WT_error=0.1)
        assert config.WT_val == 1.0
        assert config.WT_error == 0.1

    def test_wt_val_alone_is_valid(self):
        config = RheoscaleConfig(protein_name="Test", WT_val=1.0)
        assert config.WT_val == 1.0

    def test_wt_error_without_val_raises(self):
        with pytest.raises(ValueError):
            RheoscaleConfig(protein_name="Test", WT_error=0.1)


class TestRheoscaleConfigDefaults:
    """Test default values are set correctly."""

    def test_default_columns(self):
        config = RheoscaleConfig(protein_name="Test")
        assert config.columns["position"] == "Position"
        assert config.columns["substitution"] == "Substitution"
        assert config.columns["value"] == "Value"
        assert config.columns["error"] == "Error"

    def test_default_thresholds(self):
        assert RheoscaleConfig.enhancing_threshold == 0.8
        assert RheoscaleConfig.neutral_threshold == 0.7
        assert RheoscaleConfig.rheostat_threshold == 0.5
        assert RheoscaleConfig.toggle_threshold == 0.64

    def test_default_output_dir(self):
        config = RheoscaleConfig(protein_name="Test")
        assert config.output_dir == "Rheoscale_analysis"

    def test_default_even_bins(self):
        config = RheoscaleConfig(protein_name="Test")
        assert config.even_bins is True


class TestChangeColumns:
    """Test the change_colums method."""

    def test_change_single_column(self):
        config = RheoscaleConfig(protein_name="Test")
        new_config = config.change_colums(position="Pos")
        assert new_config.columns["position"] == "Pos"
        # original unchanged
        assert config.columns["position"] == "Position"

    def test_change_multiple_columns(self):
        config = RheoscaleConfig(protein_name="Test")
        new_config = config.change_colums(position="Pos", value="Score")
        assert new_config.columns["position"] == "Pos"
        assert new_config.columns["value"] == "Score"

    def test_invalid_column_key_raises(self):
        config = RheoscaleConfig(protein_name="Test")
        with pytest.raises(KeyError, match="Invalid column key"):
            config.change_colums(nonexistent="Bad")


class TestNumericOrNoneDict:
    """Test the numeric_or_none_dict method."""

    def test_returns_numeric_fields(self):
        config = RheoscaleConfig(protein_name="Test", WT_val=1.5, number_of_bins=10)
        d = config.numeric_or_none_dict()
        assert "WT_val" in d
        assert d["WT_val"] == 1.5
        assert "number_of_bins" in d
        assert d["number_of_bins"] == 10

    def test_excludes_non_numeric(self):
        config = RheoscaleConfig(protein_name="Test")
        d = config.numeric_or_none_dict()
        assert "protein_name" not in d
        assert "log_scale" not in d
        assert "columns" not in d

    def test_none_values_included(self):
        config = RheoscaleConfig(protein_name="Test")
        d = config.numeric_or_none_dict()
        assert "WT_val" in d
        assert d["WT_val"] is None


class TestFixedKeysDict:
    """Test that FixedKeysDict prevents adding new keys."""

    def test_update_existing_key(self):
        d = FixedKeysDict({"a": 1, "b": 2})
        d["a"] = 10
        assert d["a"] == 10

    def test_add_new_key_raises(self):
        d = FixedKeysDict({"a": 1})
        with pytest.raises(KeyError, match="cannot be added"):
            d["new_key"] = 99


class TestJsonRoundTrip:
    """Test JSON serialization and deserialization."""

    def test_to_and_from_json(self, tmp_path):
        config = RheoscaleConfig(
            protein_name="Test",
            WT_val=1.0,
            WT_error=0.1,
            number_of_bins=10,
        )
        json_path = str(tmp_path / "config.json")
        config.to_json(json_path)

        loaded = RheoscaleConfig.from_json(json_path)
        assert loaded.protein_name == "Test"
        assert loaded.WT_val == 1.0
        assert loaded.WT_error == 0.1
        assert loaded.number_of_bins == 10
