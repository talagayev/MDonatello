import MDAnalysis as mda
from ipywidgets import interact, Layout, VBox, HTML, Dropdown, Button, Checkbox
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import AllChem, Descriptors
from IPython.display import display, clear_output
from io import BytesIO
import base64


class MoleculeVisualizer:
    def __init__(self, ag, show_atom_indices=False, highlight_aromatic=False):
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
        self.highlight_aromatic_checkbox = Checkbox(value=highlight_aromatic, description="Highlight aromatic atoms")
        self.physiochem_props_checkbox = Checkbox(value=False, description="Show Physiochemical Properties")
        self.hbond_props_checkbox = Checkbox(value=False, description=self.get_hbond_description())
        self.save_button = Button(description="Save as PNG")
        
        # Save button click event
        self.save_button.on_click(self.save_selected_molecule)
        
        # Display widgets
        self.output_dropdown = VBox()
        self.output_dropdown.children = [
            self.dropdown, self.show_atom_indices_checkbox, self.highlight_aromatic_checkbox,
            self.physiochem_props_checkbox, self.hbond_props_checkbox
        ]
        self.output_molecule = VBox()
        self.output = VBox()
        self.output.children = [self.output_molecule, self.save_button]
        
        display(self.output_dropdown, self.output)
        
        # Update display when dropdown value changes
        self.update_display()
        
        # Link widgets to display update
        self.dropdown.observe(self.update_display, names="value")
        self.show_atom_indices_checkbox.observe(self.update_display, names="value")
        self.highlight_aromatic_checkbox.observe(self.update_display, names="value")
        self.physiochem_props_checkbox.observe(self.update_display, names="value")
        self.hbond_props_checkbox.observe(self.update_display, names="value")

    def get_hbond_description(self):
        num_h_donors = Descriptors.NumHDonors(self.ag_noh)
        num_h_acceptors = Descriptors.NumHAcceptors(self.ag_noh)
        return f"Show H-Bond Donors/Acceptors"

    def display_molecule(self, ag, show_atom_indices, highlight_aromatic):
        if highlight_aromatic:
            hit_ats = [atom.GetIdx() for atom in ag.GetAtoms() if atom.GetIsAromatic()]
            hit_bonds = [
                bond.GetIdx() for bond in ag.GetBonds() 
                if bond.GetBeginAtom().GetIsAromatic() and bond.GetEndAtom().GetIsAromatic()
            ]
            d = rdMolDraw2D.MolDraw2DSVG(300, 300)
            d.drawOptions().addStereoAnnotation = True
            rdMolDraw2D.PrepareAndDrawMolecule(
                d, ag, highlightAtoms=hit_ats, highlightBonds=hit_bonds,
                highlightAtomColors={idx: (0.5, 0, 0.5) for idx in hit_ats},
                highlightBondColors={idx: (0.5, 0, 0.5) for idx in hit_bonds}
            )
            d.FinishDrawing()
            svg = d.GetDrawingText()
            return HTML(svg)
        
        if show_atom_indices:
            for atom in ag.GetAtoms():
                atom.SetProp("atomNote", str(atom.GetIdx()))
        else:
            for atom in ag.GetAtoms():
                atom.ClearProp("atomNote")
        
        bio = BytesIO()
        img = Draw.MolToImage(ag)
        img.save(bio, format="png")
        html = f'<img src="data:image/png;base64,{base64.b64encode(bio.getvalue()).decode()}" />'
        return HTML(html)
        
    def update_display(self, _=None):
        self.hbond_props_checkbox.description = self.get_hbond_description()
        smiles = Chem.MolToSmiles(self.ag_noh)
        
        children = [
            self.display_molecule(self.ag_noh, self.show_atom_indices_checkbox.value, self.highlight_aromatic_checkbox.value),
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

