import pytest
from rdkit import Chem
from rdkit.Chem import AllChem, ChemicalFeatures
from unittest.mock import MagicMock
from mdonatello.drawer import (
    PharmacophoreHighlighter,
    FunctionalGroupHighlighter,
    RotatableBondsHighlighter,
    PartialChargeHighlighter,
    StereocenterHighlighter,
    MurckoScaffoldHighlighter,
    MoleculeDrawer,
)

@pytest.fixture
def molecule():
    smiles = "c1ccccc1"  # Benzene
    mol = Chem.MolFromSmiles(smiles)
    AllChem.Compute2DCoords(mol)
    return mol

@pytest.fixture
def checkbox_mock():
    mock = MagicMock()
    mock.value = True
    return mock

@pytest.fixture
def pharmacophore_checkboxes(checkbox_mock):
    return {
        "Aromatic": checkbox_mock,
        "Hydrophobic": checkbox_mock,
        "HBA": checkbox_mock,
        "HBD": checkbox_mock,
        "PosIonizable": checkbox_mock,
        "NegIonizable": checkbox_mock,
    }

@pytest.fixture
def functional_groups_checkboxes(checkbox_mock):
    return {"Alcohol": checkbox_mock}

@pytest.fixture
def rotatable_bonds_checkbox(checkbox_mock):
    return checkbox_mock

@pytest.fixture
def partial_charges_checkbox(checkbox_mock):
    return checkbox_mock

@pytest.fixture
def partial_charges_heatmap_checkbox(checkbox_mock):
    return checkbox_mock

@pytest.fixture
def stereocenters_checkbox(checkbox_mock):
    return checkbox_mock

@pytest.fixture
def murcko_scaffold_checkbox(checkbox_mock):
    return checkbox_mock

@pytest.fixture
def factory():
    mock_factory = MagicMock()
    mock_factory.GetFeaturesForMol = MagicMock(return_value=[])
    return mock_factory

def test_pharmacophore_highlighter(molecule, pharmacophore_checkboxes, factory):
    highlighter = PharmacophoreHighlighter(molecule, pharmacophore_checkboxes, factory)
    highlights, highlight_colors = highlighter.determine_highlights()
    assert "atoms" in highlights
    assert "bonds" in highlights
    assert isinstance(highlights["atoms"], list)
    assert isinstance(highlights["bonds"], list)
    assert isinstance(highlight_colors, dict)

def test_functional_group_highlighter(molecule, functional_groups_checkboxes):
    highlighter = FunctionalGroupHighlighter(molecule, functional_groups_checkboxes)
    highlights, highlight_colors = highlighter.determine_highlights()
    assert "atoms" in highlights
    assert "bonds" in highlights
    assert isinstance(highlights["atoms"], list)
    assert isinstance(highlights["bonds"], list)
    assert isinstance(highlight_colors, dict)

def test_rotatable_bonds_highlighter(molecule, rotatable_bonds_checkbox):
    highlighter = RotatableBondsHighlighter(molecule, rotatable_bonds_checkbox)
    highlights, highlight_colors = highlighter.determine_highlights()
    assert "atoms" in highlights
    assert "bonds" in highlights
    assert isinstance(highlights["atoms"], list)
    assert isinstance(highlights["bonds"], list)
    assert isinstance(highlight_colors, dict)

def test_partial_charge_highlighter(molecule, partial_charges_checkbox):
    highlighter = PartialChargeHighlighter(molecule, partial_charges_checkbox)
    highlights, highlight_colors = highlighter.determine_highlights()
    assert "atoms" in highlights
    assert isinstance(highlights["atoms"], list)
    assert isinstance(highlight_colors, dict)

def test_stereocenter_highlighter(molecule, stereocenters_checkbox):
    highlighter = StereocenterHighlighter(molecule, stereocenters_checkbox)
    highlights, highlight_colors = highlighter.determine_highlights()
    assert "atoms" in highlights
    assert isinstance(highlights["atoms"], list)
    assert isinstance(highlight_colors, dict)

def test_murcko_scaffold_highlighter(molecule, murcko_scaffold_checkbox):
    highlighter = MurckoScaffoldHighlighter(molecule, murcko_scaffold_checkbox)
    highlights, highlight_colors = highlighter.determine_highlights()
    assert "atoms" in highlights
    assert "bonds" in highlights
    assert isinstance(highlights["atoms"], list)
    assert isinstance(highlights["bonds"], list)
    assert isinstance(highlight_colors, dict)

def test_molecule_drawer(
    molecule,
    pharmacophore_checkboxes,
    functional_groups_checkboxes,
    rotatable_bonds_checkbox,
    partial_charges_checkbox,
    partial_charges_heatmap_checkbox,
    stereocenters_checkbox,
    murcko_scaffold_checkbox,
    factory
):
    drawer = MoleculeDrawer(
        molecule,
        pharmacophore_checkboxes,
        functional_groups_checkboxes,
        rotatable_bonds_checkbox,
        partial_charges_checkbox,
        partial_charges_heatmap_checkbox,
        stereocenters_checkbox,
        murcko_scaffold_checkbox,
        factory
    )
    svg = drawer.draw_molecule(show_atom_indices=True, width=300, height=300)
    assert isinstance(svg, str)
