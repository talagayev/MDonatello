"""
Unit and regression test for the mdonatello package.
"""

# Import package, test suite, and other packages as needed
import mdonatello
from mdonatello import MoleculeVisualizer
import pytest
import sys
from ipywidgets import VBox, HTML, Dropdown, Checkbox, Button
from rdkit import Chem
from rdkit.Chem import AllChem
from MDAnalysis import Universe

# Dummy data for testing
mol = Chem.MolFromSmiles("CCO")
AllChem.Compute2DCoords(mol)

def test_mdonatello_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "mdonatello" in sys.modules

def test_mdanalysis_logo_length(mdanalysis_logo_text):
    """Example test using a fixture defined in conftest.py"""
    logo_lines = mdanalysis_logo_text.split("\n")
    assert len(logo_lines) == 46, "Logo file does not have 46 lines!"

class MockAtomGroup:
    def convert_to(self, format):
        return mol

@pytest.fixture
def molecule_visualizer():
    ag = MockAtomGroup()
    return MoleculeVisualizer(ag)
    
def test_widget_initialization(molecule_visualizer):
    """Test Ipywidget initialization"""
    assert isinstance(molecule_visualizer.dropdown, Dropdown)
    assert isinstance(molecule_visualizer.show_atom_indices_checkbox, Checkbox)
    assert isinstance(molecule_visualizer.highlight_aromatic_checkbox, Checkbox)
    assert isinstance(molecule_visualizer.physiochem_props_checkbox, Checkbox)
    assert isinstance(molecule_visualizer.hbond_props_checkbox, Checkbox)
    assert isinstance(molecule_visualizer.save_button, Button)

def test_checkbox_interaction(molecule_visualizer):
    """Test the Ipywidgets checkbox updates after interaction"""
    assert not molecule_visualizer.show_atom_indices_checkbox.value
    molecule_visualizer.show_atom_indices_checkbox.value = True
    assert molecule_visualizer.show_atom_indices_checkbox.value

    assert not molecule_visualizer.highlight_aromatic_checkbox.value
    molecule_visualizer.highlight_aromatic_checkbox.value = True
    assert molecule_visualizer.highlight_aromatic_checkbox.value

    assert not molecule_visualizer.physiochem_props_checkbox.value
    molecule_visualizer.physiochem_props_checkbox.value = True
    assert molecule_visualizer.physiochem_props_checkbox.value
