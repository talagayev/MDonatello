import pytest
from unittest.mock import Mock
from rdkit import Chem
from rdkit.Chem import AllChem, ChemicalFeatures
from rdkit.Chem.Scaffolds import MurckoScaffold
from mdonatello.drawer import (
    PharmacophoreHighlighter,
    FunctionalGroupHighlighter,
    RotatableBondsHighlighter,
    PartialChargeHighlighter,
    StereocenterHighlighter,
    MurckoScaffoldHighlighter
)

@pytest.fixture
def molecule():
    # Create a test molecule (e.g., benzene)
    mol = Chem.MolFromSmiles('c1ccccc1')
    AllChem.Compute2DCoords(mol)
    return mol

@pytest.fixture
def mock_checkboxes():
    mock_checkboxes = {
        'Aromatic': Mock(value=True),
        'Alcohol': Mock(value=False),
        'Rotatable': Mock(value=True),
        'PartialCharge': Mock(value=True),
        'Stereocenters': Mock(value=True),
        'MurckoScaffold': Mock(value=True)
    }
    return mock_checkboxes

@pytest.fixture
def mock_factory():
    factory = Mock()
    # Mock the GetFeaturesForMol method
    factory.GetFeaturesForMol = Mock(return_value=[])
    return factory


def test_pharmacophore_highlighter(molecule, mock_checkboxes, mock_factory):
    # Setup mock features
    mock_feature = Mock()
    mock_feature.GetFamily = Mock(return_value='Aromatic')
    mock_feature.GetAtomIds = Mock(return_value=[0, 1, 2])
    mock_factory.GetFeaturesForMol = Mock(return_value=[mock_feature])

    # Initialize the highlighter
    highlighter = PharmacophoreHighlighter(molecule, mock_checkboxes, mock_factory)
    highlights, highlight_colors = highlighter.determine_highlights()

    # Assert that highlights contain atoms and bonds
    assert 'atoms' in highlights
    assert 'bonds' in highlights
    assert isinstance(highlights['atoms'], list)
    assert isinstance(highlights['bonds'], list)

    # Assert that highlights contain colors
    assert isinstance(highlight_colors, dict)
    assert all(isinstance(color, tuple) and len(color) == 3 for color in highlight_colors.values())

    # Assert that if Aromatic checkbox is checked, atoms are highlighted
    if mock_checkboxes['Aromatic'].value:
        assert len(highlights['atoms']) > 0
        assert all(atom in highlight_colors for atom in highlights['atoms'])
