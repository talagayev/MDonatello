"""
Unit and regression test for the mdonatello package.
"""

# Import package, test suite, and other packages as needed
import mdonatello
from mdonatello import MoleculeVisualizer, MolecularWeight, LogP, TPSA, RotatableBonds, HydrogenBondAcceptors, HydrogenBondDonors
import pytest
import sys
from ipywidgets import VBox, HTML, Dropdown, Checkbox, Button
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from MDAnalysis import Universe

# Dummy data for testing
mol = Chem.MolFromSmiles("CCO")
AllChem.Compute2DCoords(mol)

@pytest.fixture
def molecule_benzene():
    smiles = 'c1ccccc1'
    return Chem.MolFromSmiles(smiles)

def test_molecular_weight(molecule_benzene):
    mw = MolecularWeight(molecule_benzene)
    assert mw.property_value == Descriptors.MolWt(molecule_benzene)
    assert repr(mw) == f"<p style='margin: 0; margin-left: 100px;'>Molecular Weight: <b>{Descriptors.MolWt(benzene):.2f}</b></p>"

def test_logp(molecule_benzene):
    logp = LogP(molecule_benzene)
    assert logp.property_value == Descriptors.MolLogP(molecule_benzene)
    assert repr(logp) == f"<p style='margin: 0; margin-left: 100px;'>LogP: <b>{Descriptors.MolLogP(molecule_benzene):.2f}</b></p>"

def test_tpsa(molecule_benzene):
    tpsa = TPSA(molecule_benzene)
    assert tpsa.property_value == Descriptors.TPSA(molecule_benzene)
    assert repr(tpsa) == f"<p style='margin: 0; margin-left: 100px;'>TPSA: <b>{Descriptors.TPSA(molecule_benzene):.2f}</b></p>"

def test_rotatable_bonds(molecule_benzene):
    rb = RotatableBonds(molecule_benzene)
    assert rb.property_value == Descriptors.NumRotatableBonds(molecule_benzene)
    assert repr(rb) == f"<p style='margin: 0; margin-left: 100px;'>Rotatable Bonds: <b>{Descriptors.NumRotatableBonds(molecule_benzene)}</b></p>"

def test_hydrogen_bond_acceptors(molecule_benzene):
    hba = HydrogenBondAcceptors(molecule_benzene)
    assert hba.property_value == Descriptors.NumHAcceptors(molecule_benzene)
    assert repr(hba) == f"<p style='margin: 0; margin-left: 100px;'>Hydrogen Bond Acceptors: <b>{Descriptors.NumHAcceptors(molecule_benzene)}</b></p>"

def test_hydrogen_bond_donors(molecule_benzene):
    hbd = HydrogenBondDonors(molecule_benzene)
    assert hbd.property_value == Descriptors.NumHDonors(molecule_benzene)
    assert repr(hbd) == f"<p style='margin: 0; margin-left: 100px;'>Hydrogen Bond Donors: <b>{Descriptors.NumHDonors(molecule_benzene)}</b></p>"

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
