import MDAnalysis as mda
from ipywidgets import (
    interact,
    Layout,
    VBox,
    HTML,
    Dropdown,
    Button,
    Checkbox,
    HBox,
)
from rdkit import Chem, RDConfig
from rdkit.Chem import (
    Draw,
    AllChem,
    ChemicalFeatures,
)
from IPython.display import display, clear_output
from io import BytesIO
import base64
import logging
import os
from functools import cached_property
from mdonatello.mapper import FunctionalGroupHandler
from mdonatello.drawer import MoleculeDrawer
from mdonatello.properties import (
    MolecularWeight,
    LogP,
    TPSA,
    RotatableBonds,
    HydrogenBondAcceptors,
    HydrogenBondDonors,
    Stereocenters,
)


class MoleculeVisualizer:
    """A class for small molecule 2D visualization in jupyter notebook

    Parameters:
    -----------
    ag : MDAnalysis.core.groups.AtomGroup
        An AtomGroup object representing the molecules that need to be visualized.
    show_atom_indices : bool, optional
        Whether to display atom indices of the molecule. Default is False.
    width : int, optional
        The width of the image in pixels. Default is -1.
    height : int, optional
        The height of the image in pixels. Default is -1.

    """

    def __init__(
        self,
        ag: mda.core.groups.AtomGroup,
        show_atom_indices: bool = False,
        width: int = -1,
        height: int = -1,
    ):
        """
        Initializes the MoleculeVisualizer with an AtomGroup and visualization options.
        """
        self.width = width
        self.height = height
        self.show_atom_indices = show_atom_indices

        self.mol: Chem.Mol = ag.convert_to("RDKit")
        self.mol_noh: Chem.Mol = Chem.RemoveHs(self.mol)
        AllChem.Compute2DCoords(self.mol_noh)

        # Get individual fragments
        fragments: list[Chem.Mol] = Chem.GetMolFrags(self.mol_noh, asMols=True)
        self.molecule_list: list[str] = [
            Chem.MolToSmiles(frag) for frag in fragments
        ]
        self.fragments: dict[str, Chem.Mol] = {
            smiles: frag for smiles, frag in zip(self.molecule_list, fragments)
        }

        self.fdefName = os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef")
        self.factory = ChemicalFeatures.BuildFeatureFactory(self.fdefName)

        self.initialize_widgets()
        self.initialize_output()
        self.link_widget_callbacks()

        self.update_display()

    def initialize_widgets(self):
        """
        Initializes the interactive widgets for molecule visualization.
        """
        self.dropdown = Dropdown(
            options=self.molecule_list,
            description="Select molecule:",
            layout=Layout(width="50%"),
        )
        self.show_atom_indices_checkbox = Checkbox(
            value=self.show_atom_indices, description="Show atom indices"
        )
        self.partial_charges_checkbox = Checkbox(
            value=False, description="Show partial charges"
        )
        self.partial_charges_heatmap_checkbox = Checkbox(
            value=False, description="Show partial charge heatmap"
        )
        self.stereocenters_checkbox = Checkbox(
            value=False, description="Show Stereocenters"
        )
        self.murcko_scaffold_checkbox = Checkbox(
            value=False, description="Show Murcko Scaffold"
        )
        self.physiochem_props_checkbox = Checkbox(
            value=False, description="Show Physiochemical Properties"
        )
        self.hbond_props_checkbox = Checkbox(
            value=False, description="Show H-Bond Donors/Acceptors"
        )
        self.rotatable_bonds_checkbox = Checkbox(
            value=False, description="Show Rotatable Bonds"
        )
        self.functional_groups_checkbox = Checkbox(
            value=False, description="Show Functional Groups"
        )
        self.save_button = Button(description="Save as PNG")

        pharmacophore_families = [
            "Donor",
            "Acceptor",
            "Hydrophobe",
            "PosIonizable",
            "NegIonizable",
            "Aromatic",
            "LumpedHydrophobe",
        ]

        self.pharmacophore_checkboxes = {
            family: Checkbox(value=False, description=f"Highlight {family}")
            for family in pharmacophore_families
        }

    def initialize_output(self):
        """
        Initializes the output display and arranges the widgets in the interface.
        """
        properties_header = HTML("<h3>Properties</h3>")
        highlighting_header = HTML("<h3>Atom & Bond Highlighting</h3>")
        pharmacophores_header = HTML("<h3>Pharmacophores</h3>")
        molecule_header = HTML("<h3>Molecule</h3>")

        pharmacophore_checkbox_rows = []
        checkboxes = list(self.pharmacophore_checkboxes.values())
        for i in range(0, len(checkboxes), 4):
            row = checkboxes[i : i + 4]
            pharmacophore_checkbox_rows.append(HBox(row))

        self.output_dropdown = VBox()
        self.output_dropdown.children = (
            [
                HBox([self.dropdown]),
                properties_header,
                HBox(
                    [
                        self.physiochem_props_checkbox,
                        self.partial_charges_checkbox,
                        self.hbond_props_checkbox,
                        self.show_atom_indices_checkbox,
                    ]
                ),
                highlighting_header,
                HBox(
                    [
                        self.rotatable_bonds_checkbox,
                        self.partial_charges_heatmap_checkbox,
                        self.functional_groups_checkbox,
                        self.stereocenters_checkbox,
                    ]
                ),
                HBox(
                    [
                        self.murcko_scaffold_checkbox,
                    ]
                ),
                pharmacophores_header,
            ]
            + pharmacophore_checkbox_rows
            + [molecule_header]
        )

        self.output_molecule = VBox()
        self.output = VBox()
        self.output.children = [self.output_molecule, self.save_button]

        display(self.output_dropdown, self.output)

    def link_widget_callbacks(self):
        """
        Links the interactive widgets to the update_display method to reflect changes.
        """
        self.dropdown.observe(self.update_display, names="value")
        self.show_atom_indices_checkbox.observe(
            self.update_display, names="value"
        )
        self.partial_charges_checkbox.observe(
            self.update_display, names="value"
        )
        self.partial_charges_heatmap_checkbox.observe(
            self.update_display, names="value"
        )
        self.stereocenters_checkbox.observe(self.update_display, names="value")
        self.murcko_scaffold_checkbox.observe(
            self.update_display, names="value"
        )
        self.physiochem_props_checkbox.observe(
            self.update_display, names="value"
        )
        self.hbond_props_checkbox.observe(self.update_display, names="value")
        self.rotatable_bonds_checkbox.observe(
            self.update_display, names="value"
        )
        self.functional_groups_checkbox.observe(
            self.update_display, names="value"
        )
        for checkbox in self.pharmacophore_checkboxes.values():
            checkbox.observe(self.update_display, names="value")
        self.save_button.on_click(self.save_selected_molecule)

    def update_display(self, _=None):
        """
        Updates the molecule display based on the current selections.

        Parameters
        ----------
        _ : any, optional
            A placeholder parameter for widget callback compatibility.
        """
        smiles = self.dropdown.value
        self.current_mol = self.fragments[smiles]

        # Update functional group checkboxes dynamically
        fg_counts = FunctionalGroupHandler.calculate_functional_groups(
            self.current_mol
        )
        self.functional_group_checkboxes = {}
        for fg, atom_indices in fg_counts.items():
            if atom_indices:
                fg_checkbox_name = f"fg_checkbox_{fg}"
                if not hasattr(self, fg_checkbox_name):
                    checkbox = Checkbox(value=False, description=fg)
                    setattr(self, fg_checkbox_name, checkbox)
                    checkbox.observe(self.update_display, names="value")
                self.functional_group_checkboxes[fg] = getattr(
                    self, fg_checkbox_name
                )

        drawer = MoleculeDrawer(
            molecule=self.current_mol,
            pharmacophore_checkboxes=self.pharmacophore_checkboxes,
            functional_groups_checkboxes=self.functional_group_checkboxes,
            rotatable_bonds_checkbox=self.rotatable_bonds_checkbox,
            partial_charges_checkbox=self.partial_charges_checkbox,
            partial_charges_heatmap_checkbox=self.partial_charges_heatmap_checkbox,
            stereocenters_checkbox=self.stereocenters_checkbox,
            murcko_scaffold_checkbox=self.murcko_scaffold_checkbox,
            factory=self.factory,
        )

        children = [
            HTML(
                drawer.draw_molecule(
                    self.show_atom_indices_checkbox.value,
                    self.width,
                    self.height,
                )
            ),
            HTML(f"<h3 style='margin: 0;'>SMILES: {smiles}</h3>"),
        ]

        if self.physiochem_props_checkbox.value:
            physiochem_properties = [
                MolecularWeight(self.current_mol),
                LogP(self.current_mol),
                TPSA(self.current_mol),
                RotatableBonds(self.current_mol),
                Stereocenters(self.current_mol),
            ]
            physiochem_html = [
                HTML(str(prop)) for prop in physiochem_properties
            ]
            children.extend(physiochem_html)

        if self.hbond_props_checkbox.value:
            hbond_properties = [
                HydrogenBondAcceptors(self.current_mol),
                HydrogenBondDonors(self.current_mol),
            ]
            hbond_html = [HTML(str(prop)) for prop in hbond_properties]
            children.extend(hbond_html)

        # Show or hide functional group checkboxes based on the main functional groups checkbox
        if self.functional_groups_checkbox.value:
            functional_groups_header = HTML("<h3>Functional Groups</h3>")
            fg_checkboxes = list(self.functional_group_checkboxes.values())
            if fg_checkboxes:
                fg_hbox = HBox(fg_checkboxes)
                children.append(functional_groups_header)
                children.append(fg_hbox)

        self.output_molecule.children = children

    def save_selected_molecule(self, _):
        """
        Saves the currently selected molecule as a PNG file.

        Parameters
        ----------
        _ : any, optional
            A placeholder parameter for button callback compatibility.
        """
        smiles = self.dropdown.value
        mol = self.fragments[smiles]
        filename = f"{smiles}.png"
        img = Draw.MolToImage(mol)
        img.save(filename)
        logging.info(f"Molecule saved as '{filename}'")
