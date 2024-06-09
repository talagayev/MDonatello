import MDAnalysis as mda
from ipywidgets import interact, Layout, VBox, HTML, Dropdown, Button, Checkbox
from rdkit import Chem, RDConfig
from rdkit.Chem import Draw, AllChem, Descriptors, ChemicalFeatures
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import display, clear_output
from io import BytesIO
import base64
import os


class MoleculeVisualizer:
    def __init__(self, ag, show_atom_indices=False):
        self.ag = ag.convert_to("RDKit")
        self.ag_noh = Chem.RemoveHs(self.ag)
        AllChem.Compute2DCoords(self.ag_noh)
        self.molecule_list = ["Molecule"]
        
        # Create the dropdown and other widgets
        self.dropdown = Dropdown(
            options=self.molecule_list,
            description="Select molecule:",
            layout=Layout(width="50%")
        )
        self.show_atom_indices_checkbox = Checkbox(value=show_atom_indices, description="Show atom indices")
        self.physiochem_props_checkbox = Checkbox(value=False, description="Show Physiochemical Properties")
        self.hbond_props_checkbox = Checkbox(value=False, description="Show H-Bond Donors/Acceptors")
        self.save_button = Button(description="Save as PNG")

        # Pharmacophore feature detection
        self.fdefName = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
        self.factory = ChemicalFeatures.BuildFeatureFactory(self.fdefName)
        self.feats = self.factory.GetFeaturesForMol(self.ag_noh)

        # Dynamically create checkboxes for each unique pharmacophore type
        self.pharmacophore_checkboxes = {}
        for feat in self.feats:
            family = feat.GetFamily()
            if family not in self.pharmacophore_checkboxes:
                self.pharmacophore_checkboxes[family] = Checkbox(value=False, description=f"Highlight {family}")
        
        # Save button click event
        self.save_button.on_click(self.save_selected_molecule)
        
        # Display widgets
        self.output_dropdown = VBox()
        self.output_dropdown.children = [
            self.dropdown, self.show_atom_indices_checkbox,
            self.physiochem_props_checkbox, self.hbond_props_checkbox
        ] + list(self.pharmacophore_checkboxes.values())
        self.output_molecule = VBox()
        self.output = VBox()
        self.output.children = [self.output_molecule, self.save_button]
        
        display(self.output_dropdown, self.output)
        
        # Update display when dropdown value changes
        self.update_display()
        
        # Link widgets to display update
        self.dropdown.observe(self.update_display, names="value")
        self.show_atom_indices_checkbox.observe(self.update_display, names="value")
        self.physiochem_props_checkbox.observe(self.update_display, names="value")
        self.hbond_props_checkbox.observe(self.update_display, names="value")
        for checkbox in self.pharmacophore_checkboxes.values():
            checkbox.observe(self.update_display, names="value")

    def display_molecule(self, ag, show_atom_indices):
        highlights = {"atoms": [], "bonds": []}
        highlight_colors = {}

        # Pharmacophore highlighting
        for feat in self.feats:
            family = feat.GetFamily()
            if self.pharmacophore_checkboxes[family].value:
                atom_ids = feat.GetAtomIds()
                highlights["atoms"].extend(atom_ids)
                color = self.get_color_for_pharmacophore(family)
                for atom_id in atom_ids:
                    highlight_colors[atom_id] = color

        if show_atom_indices:
            for atom in ag.GetAtoms():
                atom.SetProp("atomNote", str(atom.GetIdx()))
        else:
            for atom in ag.GetAtoms():
                atom.ClearProp("atomNote")

        d = rdMolDraw2D.MolDraw2DSVG(300, 300)
        d.drawOptions().addStereoAnnotation = True
        rdMolDraw2D.PrepareAndDrawMolecule(
            d, ag, highlightAtoms=highlights["atoms"], highlightBonds=highlights["bonds"],
            highlightAtomColors=highlight_colors, highlightBondColors=highlight_colors
        )
        d.FinishDrawing()
        svg = d.GetDrawingText()
        return HTML(svg)
        
    def get_color_for_pharmacophore(self, family):
        color_map = {
            "Donor": (0.0, 1.0, 0.0),      # Green
            "Acceptor": (1.0, 0.7, 0.7),   # Rosa
            "Hydrophobe": (1.0, 1.0, 0.0),  # Yellow
            "PosIonizable": (0.0, 1.0, 1.0),  # Turquoise
            "NegIonizable": (1.0, 0.0, 1.0),  # Pink
            "Aromatic": (0.5, 0.5, 1.0),  # Light Blue
            "LumpedHydrophobe": (1.0, 0.5, 0.0)  # Orange
        }
        return color_map.get(family, (0.5, 0.5, 0.5))  # Default to grey if not specified
    
    def update_display(self, _=None):
        smiles = Chem.MolToSmiles(self.ag_noh)
        
        children = [
            self.display_molecule(self.ag_noh, self.show_atom_indices_checkbox.value),
            HTML(f"<h3 style='margin: 0;'>SMILES: {smiles}</h3>")
        ]
        
        if self.physiochem_props_checkbox.value:
            children.extend([
                self.display_molecular_weight(self.ag_noh),
                self.display_logp(self.ag_noh),
                self.display_tpsa(self.ag_noh),
                self.display_rotatable_bonds(self.ag_noh)
            ])
            
        if self.hbond_props_checkbox.value:
            children.extend([
                self.display_num_h_donors(self.ag_noh),
                self.display_num_h_acceptors(self.ag_noh)
            ])
        
        self.output_molecule.children = children
        
    def display_molecular_weight(self, ag):
        mw = Descriptors.MolWt(ag)
        return HTML("<p style='margin: 0; margin-left: 100px;'>Molecular Weight: {:.2f} g/mol</p>".format(mw))
        
    def display_logp(self, ag):
        logp = Descriptors.MolLogP(ag)
        return HTML("<p style='margin: 0; margin-left: 100px;'>LogP: {:.2f}</p>".format(logp))

    def display_num_h_donors(self, ag):
        num_h_donors = Descriptors.NumHDonors(ag)
        return HTML("<p style='margin: 0; margin-left: 100px;'>Number of H-Bond Donors: {:.0f}</p>".format(num_h_donors))

    def display_num_h_acceptors(self, ag):
        num_h_acceptors = Descriptors.NumHAcceptors(ag)
        return HTML("<p style='margin: 0; margin-left: 100px;'>Number of H-Bond Acceptors: {:.0f}</p>".format(num_h_acceptors))
    
    def display_tpsa(self, ag):
        tpsa = Descriptors.TPSA(ag)
        return HTML("<p style='margin: 0; margin-left: 100px;'>Topological Polar Surface Area (TPSA): {:.2f} Å²</p>".format(tpsa))
        
    def display_rotatable_bonds(self, ag):
        rotatable_bonds = Descriptors.NumRotatableBonds(ag)
        return HTML("<p style='margin: 0; margin-left: 100px;'>Number of Rotatable Bonds: {:.0f}</p>".format(rotatable_bonds))
        
    def save_selected_molecule(self, _):
        filename = "molecule.png"
        img = Draw.MolToImage(self.ag_noh)
        img.save(filename)
        print(f"Molecule saved as '{filename}'")
