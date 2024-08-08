import pytest
from rdkit import Chem
from mdonatello.mapper import FunctionalGroupHandler, PharmacophoreColorMapper


@pytest.fixture
def ethanol_molecule():
    """Fixture for creating an ethanol molecule for testing."""
    return Chem.MolFromSmiles("CCO")  # Ethanol


@pytest.fixture
def functional_group_handler():
    """Fixture for creating an instance of FunctionalGroupHandler."""
    return FunctionalGroupHandler()


def test_calculate_functional_groups(
    functional_group_handler, ethanol_molecule
):
    """Test functional group detection in the molecule."""
    fg_counts = functional_group_handler.calculate_functional_groups(
        ethanol_molecule
    )

    # Assert the expected functional group counts
    assert "Hydroxyl group (-OH)" in fg_counts
    assert fg_counts["Hydroxyl group (-OH)"] == [2]  # Oxygen atom index

    assert "Primary amine (-NH2)" in fg_counts
    assert fg_counts["Primary amine (-NH2)"] == []  # No amine in ethanol

    assert "Carboxyl group (-COOH)" in fg_counts
    assert fg_counts["Carboxyl group (-COOH)"] == []  # No carboxyl group


def test_get_color_for_functional_group(functional_group_handler):
    """Test color assignment for functional groups."""
    color = functional_group_handler.get_color_for_functional_group(
        "Hydroxyl group (-OH)"
    )
    assert color == (1.0, 0.7, 0.7)  # Expected color for oxygen

    color = functional_group_handler.get_color_for_functional_group(
        "Primary amine (-NH2)"
    )
    assert color == (0.5, 0.5, 1.0)  # Expected color for nitrogen


def test_get_color_for_pharmacophore():
    """Test color assignment for pharmacophore features."""
    color = PharmacophoreColorMapper.get_color_for_pharmacophore("Donor")
    assert color == (0.0, 1.0, 0.0)  # Expected green for donor

    color = PharmacophoreColorMapper.get_color_for_pharmacophore("Acceptor")
    assert color == (1.0, 0.7, 0.7)  # Expected rosa for acceptor


if __name__ == "__main__":
    pytest.main()
