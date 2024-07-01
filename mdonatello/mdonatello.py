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
from rdkit.Chem import Draw, AllChem, Descriptors, ChemicalFeatures
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import display, clear_output
from io import BytesIO
import base64
import os
from functools import cached_property


class Property:
    name = 'property'
    
    def __init__(self, mol):
        self.mol = mol
    
    @cached_property
    def property_value(self):
        raise NotImplementedError('Subclasses should implement this.')
    
    def __repr__(self):
        return f"<p style='margin: 0; margin-left: 100px;'>Number of {self.name}: {self.property_value:.2f}</p>"

class MolecularWeight(Property):
    name = 'Molecular Weight'

    @cached_property
    def property_value(self):
        return Descriptors.MolWt(self.mol)

class LogP(Property):
    name = 'LogP'

    @cached_property
    def property_value(self):
        return Descriptors.MolLogP(self.mol)

class TPSA(Property):
    name = 'TPSA'

    @cached_property
    def property_value(self):
        return Descriptors.TPSA(self.mol)

class RotatableBonds(Property):
    name = 'Rotatable Bonds'

    @cached_property
    def property_value(self):
        return Descriptors.NumRotatableBonds(self.mol)

class HydrogenBondAcceptors(Property):
    name = 'Hydrogen Bond Acceptors'

    @cached_property
    def property_value(self):
        return Descriptors.NumHAcceptors(self.mol)

class HydrogenBondDonors(Property):
    name = 'Hydrogen Bond Donors'

    @cached_property
    def property_value(self):
        return Descriptors.NumHDonors(self.mol)

class PharmacophoreColorMapper:
    @staticmethod
    def get_color_for_pharmacophore(family):
        color_map = {
            "Donor": (0.0, 1.0, 0.0),    # Green
            "Acceptor": (1.0, 0.7, 0.7),  # Rosa
            "Hydrophobe": (1.0, 1.0, 0.0),  # Yellow
            "PosIonizable": (0.0, 1.0, 1.0),  # Turquoise
            "NegIonizable": (1.0, 0.0, 1.0),  # Pink
            "Aromatic": (0.5, 0.5, 1.0),  # Light Blue
            "LumpedHydrophobe": (1.0, 0.5, 0.0),  # Orange
        }
        return color_map.get(family, (0.5, 0.5, 0.5))  # Default to grey if not specified

class FunctionalGroupHandler:
    @staticmethod
    def calculate_functional_groups(mol):
        functional_groups = {
            'Hydroxyl group (-OH)': '[OX2H]',
            'Primary amine (-NH2)': '[NX3H2]',
            'Primary ammonium (-NH3+)': '[+NX4;H3]',
            'Secondary amine (-NH-)': '[NX3H][#6]',
            'Tertiary amine (-N<)': '[NX3;H0]([#6])[#6]',
            'Carboxyl group (-COOH)': 'C(=O)[OX2H1]',
            'Ester (-COOR)': 'C(=O)[OX2H0][#6]',
            'Amide (-CON-)': 'C(=O)[NX3]',
            'Aldehyde (-CHO)': '[CX3H1](=O)[#6]',
            'Ketone (C=O)': '[CX3](=O)[#6]',
            'Ether (R-O-R)': '[#6][OX2][#6]',
            'Thiocarbonyl group (C=S)': 'C(=S)',
            'Imine group (-C=N-)': '[CX3](=N)',
            'Hydroxylamine group (-N(OH))': '[NX3][OX2H]',
            'Thiol group (-SH)': '[SX2H]',
            'Azide group (-N3)': 'N=[NX1]=[NX1]',
            'Furan ring': 'c1occc1',
            'Guanidine group (-C(=NH)(N)(NH2))': 'C(=N)(N)[NH2]',
            'Isothiocyanate (-N=C=S)': '[NX2]=C=[SX2]',
            'Isocyanate (-N=C=O)': '[NX2]=C=[OX1]',
            'Lactone (C=O-O)': '[CX3](=O)[OX2][CX3](=O)',
            'Lactam (C=O-N)': '[CX3](=O)[NX3][CX3](=O)',
            'Methoxy group (-OCH3)': '[OX2][CH3]',
            'Nitro group (-NO2)': '[NX3](=O)=O',
            'Nitroso group (-NO)': '[NX2]=O',
            'Oxazole ring': 'c1noccc1',
            'Oxime group (-C=N-OH)': '[CX3](=N[OX2H])',
            'Epoxide': 'C1CO1',
            'Nitrile': 'C#N',
            'Sulfone': 'S(=O)(=O)([#6])([#6])',
            'Sulfonamide': 'S(=O)(=O)([#6])N',
            'Sulfide': '[SX2]',
            'Urea': 'C(=O)(N)(N)',
            'Phosphoric Ester': 'P(=O)(O)([OX2H0;R1])',
            'Phosphoric Acid': 'P(=O)(O)(O)'
        }
        
        fg_counts = {}
        for fg, smarts in functional_groups.items():
            substruct_matches = mol.GetSubstructMatches(
                Chem.MolFromSmarts(smarts)
            )
            fg_counts[fg] = [
                atom_idx for match in substruct_matches for atom_idx in match
            ]
        return fg_counts
    
    @staticmethod
    def get_color_for_functional_group(fg):
        parts = fg.split('(')
        smarts_part = parts[1]
        if "P" in smarts_part:
            return (1.0, 0.5, 0.0)  # Orange for phosphore containing groups
        elif "S" in smarts_part:
            return (1.0, 1.0, 0.0)  # Yellow for sulfure containing groups
        elif "N" in smarts_part:
            return (0.5, 0.5, 1.0)  # Light blue for nitrogen containing groups
        elif "O" in smarts_part:
            return (1.0, 0.7, 0.7)  # Red for oxygen containing groups
        else:
            return (1.0, 0.5, 0.0)  # Pink if no specific color assigned

class MoleculeVisualizer:
    """A class for small molecule 2D visualization in jupyter notebook

    Parameters:
    -----------
    ag : MDAnalysis.core.groups.AtomGroup
        An AtomGroup object representing the molecules that need to be visualized.
    show_atom_indices : bool, optional
        Whether to display atom indices of the molecule. Default is False.
    width : int, optional
        The width of the image in pixels. Default is 300.
    height : int, optional
        The height of the image in pixels. Default is 300. 
        
    """
   
    def __init__(self, ag: mda.core.groups.AtomGroup, show_atom_indices: bool = False, width: int = -1, height: int = -1):
        self.mol: Chem.Mol = ag.convert_to("RDKit")
        self.mol_noh: Chem.Mol = Chem.RemoveHs(self.mol)
        AllChem.Compute2DCoords(self.mol_noh)

        # Get individual fragments
        fragments: list[Chem.Mol] = Chem.GetMolFrags(self.mol_noh, asMols=True)
        self.molecule_list: list[str] = [Chem.MolToSmiles(frag) for frag in fragments]
        self.fragments: dict[str, Chem.Mol] = {
            smiles: frag for smiles, frag in zip(self.molecule_list, fragments)
        }

        # Add height and width
        self.width: int = width
        self.height: int = height
        
        # Create the dropdown and other widgets
        self.dropdown = Dropdown(
            options=self.molecule_list,
            description="Select molecule:",
            layout=Layout(width="50%")
        )
        self.show_atom_indices_checkbox = Checkbox(
            value=show_atom_indices, description="Show atom indices"
        )
        self.physiochem_props_checkbox = Checkbox(
            value=False, description="Show Physiochemical Properties"
        )
        self.hbond_props_checkbox = Checkbox(
            value=False, description="Show H-Bond Donors/Acceptors"
        )
        self.functional_groups_checkbox = Checkbox(
            value=False, description="Show Functional Groups"
        )
        self.save_button = Button(description="Save as PNG")

        # Pharmacophore feature detection
        self.fdefName = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
        self.factory = ChemicalFeatures.BuildFeatureFactory(self.fdefName)
        self.pharmacophore_checkboxes = {}

        # Define all possible pharmacophore features
        pharmacophore_families = [
            "Donor",
            "Acceptor",
            "Hydrophobe",
            "PosIonizable",
            "NegIonizable",
            "Aromatic",
            "LumpedHydrophobe",
        ]

        for family in pharmacophore_families:
            self.pharmacophore_checkboxes[family] = Checkbox(
                value=False, description=f"Highlight {family}"
            )
        
        # Save button click event
        self.save_button.on_click(self.save_selected_molecule)
        
        # Display widgets
        pharmacophore_checkbox_rows = []
        checkboxes = list(self.pharmacophore_checkboxes.values())
        for i in range(0, len(checkboxes), 4):
            row = checkboxes[i:i+4]
            pharmacophore_checkbox_rows.append(HBox(row))
        
        properties_header = HTML("<h3>Properties</h3>")
        pharmacophores_header = HTML("<h3>Pharmacophores</h3>")
        
        self.output_dropdown = VBox()
        self.output_dropdown.children = [
            HBox([self.dropdown, self.show_atom_indices_checkbox]),
            properties_header,
            HBox(
                [
                    self.physiochem_props_checkbox,
                    self.hbond_props_checkbox,
                    self.functional_groups_checkbox,
                ]
            ),
            pharmacophores_header,
        ] + pharmacophore_checkbox_rows
        
        self.output_molecule = VBox()
        self.output = VBox()
        self.output.children = [self.output_molecule, self.save_button]
        
        display(self.output_dropdown, self.output)
        
        # Update display when dropdown value changes
        self.update_display()
        
        # Link widgets to display update
        self.dropdown.observe(self.update_display, names="value")
        self.show_atom_indices_checkbox.observe(
            self.update_display, names="value"
        )
        self.physiochem_props_checkbox.observe(
            self.update_display, names="value"
        )
        self.hbond_props_checkbox.observe(self.update_display, names="value")
        self.functional_groups_checkbox.observe(
            self.update_display, names="value"
        )
        for checkbox in self.pharmacophore_checkboxes.values():
            checkbox.observe(self.update_display, names="value")
    
    def draw_molecule(self, mol, show_atom_indices, width, height):
        highlights = {"atoms": [], "bonds": []}
        highlight_colors = {}

        # Update pharmacophore features for the selected molecule
        feats = self.factory.GetFeaturesForMol(mol)

        # Pharmacophore highlighting
        for feat in feats:
            family = feat.GetFamily()
            if self.pharmacophore_checkboxes[family].value:
                atom_ids = feat.GetAtomIds()
                highlights["atoms"].extend(atom_ids)
                color = PharmacophoreColorMapper.get_color_for_pharmacophore(family)
                for atom_id in atom_ids:
                    highlight_colors[atom_id] = color

        # Specific highlighting for Aromatic pharmacophore
        if self.pharmacophore_checkboxes["Aromatic"].value:
            hit_ats = [
                atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()
            ]
            hit_bonds = [
                bond.GetIdx()
                for bond in mol.GetBonds()
                if bond.GetBeginAtom().GetIsAromatic()
                and bond.GetEndAtom().GetIsAromatic()
            ]
            highlights["atoms"].extend(hit_ats)
            highlights["bonds"].extend(hit_bonds)
            color = PharmacophoreColorMapper.get_color_for_pharmacophore("Aromatic")
            for atom_id in hit_ats:
                highlight_colors[atom_id] = color
            for bond_id in hit_bonds:
                highlight_colors[bond_id] = color

        # Functional group highlighting
        if self.functional_groups_checkbox.value:
            fg_counts = FunctionalGroupHandler.calculate_functional_groups(mol)
            for fg, atom_indices in fg_counts.items():
                fg_checkbox_name = f"fg_checkbox_{fg}"
                if (
                    hasattr(self, fg_checkbox_name)
                    and getattr(self, fg_checkbox_name).value
                ):
                    highlights["atoms"].extend(atom_indices)
                    color = FunctionalGroupHandler.get_color_for_functional_group(fg)
                    for atom_id in atom_indices:
                        highlight_colors[atom_id] = color
                    # Highlight bonds between functional group atoms
                    for bond in mol.GetBonds():
                        if (
                            bond.GetBeginAtomIdx() in atom_indices
                            and bond.GetEndAtomIdx() in atom_indices
                        ):
                            highlights["bonds"].append(bond.GetIdx())
                            highlight_colors[bond.GetIdx()] = color
        
        d = rdMolDraw2D.MolDraw2DSVG(width, height)
        d.drawOptions().addAtomIndices = show_atom_indices
        d.drawOptions().addStereoAnnotation = True
        rdMolDraw2D.PrepareAndDrawMolecule(
            d,
            mol,
            highlightAtoms=highlights["atoms"],
            highlightBonds=highlights["bonds"],
            highlightAtomColors=highlight_colors,
            highlightBondColors=highlight_colors,
        )
        d.FinishDrawing()
        svg = d.GetDrawingText()
        return svg
        
    def update_display(self, _=None):
        smiles = self.dropdown.value
        self.current_mol = self.fragments[smiles]

        children = [
            HTML(
                self.draw_molecule(
                    self.current_mol,
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
                RotatableBonds(self.current_mol)
            ]
            physiochem_html = [HTML(repr(prop)) for prop in physiochem_properties]
            children.extend(physiochem_html)
            
        if self.hbond_props_checkbox.value:
            hbond_properties = [
                HydrogenBondAcceptors(self.current_mol),
                HydrogenBondDonors(self.current_mol),
            ]
            hbond_html = [HTML(repr(prop)) for prop in hbond_properties]
            children.extend(hbond_html)

        if self.functional_groups_checkbox.value:
            functional_groups_header = HTML("<h3>Functional Groups</h3>")
            fg_counts = self.FunctionalGroupHandler.calculate_functional_groups(self.current_mol)
            fg_checkboxes = []
            for fg, atom_indices in fg_counts.items():
                if atom_indices:
                    fg_checkbox_name = f"fg_checkbox_{fg}"
                    if not hasattr(self, fg_checkbox_name):
                        checkbox = Checkbox(value=False, description=fg)
                        setattr(self, fg_checkbox_name, checkbox)
                        checkbox.observe(self.update_display, names="value")
                    fg_checkboxes.append(getattr(self, fg_checkbox_name))
            
            if fg_checkboxes:
                fg_hbox = HBox(fg_checkboxes)
                children.append(functional_groups_header)
                children.append(fg_hbox)
        
        self.output_molecule.children = children
        
    def save_selected_molecule(self, _):
        smiles = self.dropdown.value
        mol = self.fragments[smiles]
        filename = f"{smiles}.png"
        img = Draw.MolToImage(mol)
        img.save(filename)
        print(f"Molecule saved as '{filename}'")
