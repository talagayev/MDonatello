import pytest
import MDAnalysis as mda
from rdkit import Chem
from ipywidgets import HTML
from mdonatello import MoleculeVisualizer
import os

@pytest.fixture
def ligand_atoms():
    """Creates an MDAnalysis AtomGroup for the residue name UNK from input.pdb."""
    u = mda.Universe("mdonatello/data/input.pdb")
    ag = u.select_atoms("resname UNK")
    return ag

def test_initialization(ligand_atoms):
    visualizer = MoleculeVisualizer(ligand_atoms)
    assert visualizer.mol is not None
    assert visualizer.mol_noh is not None
    assert isinstance(visualizer.molecule_list, list)
    assert isinstance(visualizer.fragments, dict)

def test_widget_initialization(ligand_atoms):
    visualizer = MoleculeVisualizer(ligand_atoms)
    assert visualizer.dropdown is not None
    assert visualizer.show_atom_indices_checkbox is not None
    assert visualizer.partial_charges_checkbox is not None
    assert visualizer.save_button is not None

def test_display_update(ligand_atoms):
    visualizer = MoleculeVisualizer(ligand_atoms)
    visualizer.update_display()
    assert visualizer.current_mol is not None
    assert visualizer.output_molecule.children != []

def test_save_molecule(ligand_atoms):
    visualizer = MoleculeVisualizer(ligand_atoms)
    visualizer.dropdown.value = visualizer.molecule_list[0]
    visualizer.save_selected_molecule(None)
    smiles = visualizer.dropdown.value
    filename = f"{smiles}.png"
    assert os.path.isfile(filename)
    os.remove(filename)

def test_functional_groups(ligand_atoms):
    visualizer = MoleculeVisualizer(ligand_atoms)
    visualizer.functional_groups_checkbox.value = True
    visualizer.update_display()
    assert visualizer.functional_group_checkboxes != {}

def test_physiochem_properties(ligand_atoms):
    visualizer = MoleculeVisualizer(ligand_atoms)
    visualizer.physiochem_props_checkbox.value = True
    visualizer.update_display()
    properties_found = any(
        isinstance(child, HTML) and "Molecular Weight" in child.value
        for child in visualizer.output_molecule.children
    )
    assert properties_found

def test_hbond_properties(ligand_atoms):
    visualizer = MoleculeVisualizer(ligand_atoms)
    visualizer.hbond_props_checkbox.value = True
    visualizer.update_display()
    hbond_found = any(
        isinstance(child, HTML) and "Hydrogen Bond" in child.value
        for child in visualizer.output_molecule.children
    )
    assert hbond_found
