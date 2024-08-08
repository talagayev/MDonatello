from typing import Dict, List, Tuple
from rdkit import Chem
from rdkit.Chem import (
    Draw,
    AllChem,
    ChemicalFeatures,
    Lipinski,
    Scaffolds,
)
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Lipinski import RotatableBondSmarts
from mdonatello.mapper import (
    PharmacophoreColorMapper,
    FunctionalGroupHandler,
)


class PharmacophoreHighlighter:
    """A class responsible for highlighting of pharmacophore features in a molecule.

    Parameters:
    -----------
    molecule : RDKit.Chem.Mol
        An RDKit molecule object representing the molecule that should be highlighted.
    pharmacophore_checkboxes : dict
        A dictionary of checkbox widgets indicating which pharmacophore features to highlight.
        Keys are pharmacophore feature names (e.g., "Aromatic"), and values are the corresponding checkbox widgets.
    factory : RDKit.Chem.ChemicalFeatures.MolChemicalFeatureFactory
        A factory object used to generate pharmacophore features available in RDKit for the given molecule.

    Returns:
    --------
    highlights : dict
        A dictionary with keys "atoms" and "bonds", each containing lists of atom indices and bond indices
        that should be highlighted, respectively.
    highlight_colors : dict
        A dictionary mapping atom indices and bond indices to their respective colors that should be used for highlighting.
    """

    def __init__(
        self,
        molecule: Chem.Mol,
        pharmacophore_checkboxes: Dict[str, "Checkbox"],
        factory: ChemicalFeatures.MolChemicalFeatureFactory,
    ):
        """Initialize the PharmacophoreHighlighter with a molecule, pharmacophore checkboxes, and an RDKit feature factory."""
        self.molecule = molecule
        self.pharmacophore_checkboxes = pharmacophore_checkboxes
        self.factory = factory

    def determine_highlights(
        self,
    ) -> Tuple[Dict[str, List[int]], Dict[int, str]]:
        """Determines the atom and bond indices to highlight them based on the selected pharmacophore features.

        Returns:
        --------
        highlights : dict
            A dictionary with keys "atoms" and "bonds", each containing lists of atom indices and bond indices
            that should be highlighted.
        highlight_colors : dict
            A dictionary mapping atom indices and bond indices to their respective highlight colors.
        """  
        mol = self.molecule
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
                color = PharmacophoreColorMapper.get_color_for_pharmacophore(
                    family
                )
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
            color = PharmacophoreColorMapper.get_color_for_pharmacophore(
                "Aromatic"
            )
            for atom_id in hit_ats:
                highlight_colors[atom_id] = color
            for bond_id in hit_bonds:
                highlight_colors[bond_id] = color

        return highlights, highlight_colors


class FunctionalGroupHighlighter:
    """A class responsible for highlighting of the functional groups in a molecule.

    Parameters:
    -----------
    molecule : RDKit.Chem.Mol
        An RDKit molecule object representing the molecule that should be highlighted.
    checkboxes : dict
        A dictionary of checkbox widgets indicating which functional groups to highlight.
        Keys are functional group names (e.g., "Alcohol"), and values are the corresponding checkbox widgets.

    Returns:
    --------
    highlights : dict
        A dictionary with keys "atoms" and "bonds", each containing lists of atom indices and bond indices
        that should be highlighted, respectively.
    highlight_colors : dict
        A dictionary mapping atom indices and bond indices to their respective highlight colors.
    """

    def __init__(self, molecule: Chem.Mol, checkboxes: Dict[str, "Checkbox"]):
        """Initialize the FunctionalGroupHighlighter with a molecule and functional group checkboxes."""
        self.molecule = molecule
        self.checkboxes = checkboxes

    def determine_highlights(
        self,
    ) -> Tuple[Dict[str, List[int]], Dict[int, str]]:
        """Determine the atom and bond indices to highlight based on the selected functional groups.

        Returns:
        --------
        highlights : dict
            A dictionary with keys "atoms" and "bonds", each containing lists of atom indices and bond indices
            that should be highlighted.
        highlight_colors : dict
            A dictionary mapping atom indices and bond indices to their respective highlight colors.
        """
        highlights = {"atoms": [], "bonds": []}
        highlight_colors = {}

        fg_counts = FunctionalGroupHandler.calculate_functional_groups(
            self.molecule
        )
        for fg, atom_indices in fg_counts.items():
            if atom_indices and self.checkboxes[fg].value:
                highlights["atoms"].extend(atom_indices)

                atom_index_pairs = [
                    (atom_indices[i], atom_indices[j])
                    for i in range(len(atom_indices))
                    for j in range(i + 1, len(atom_indices))
                ]

                for idx1, idx2 in atom_index_pairs:
                    bond = self.molecule.GetBondBetweenAtoms(idx1, idx2)
                    if bond is not None:
                        highlights["bonds"].append(bond.GetIdx())

                highlight_color = (
                    FunctionalGroupHandler.get_color_for_functional_group(fg)
                )
                for atom_idx in atom_indices:
                    highlight_colors[atom_idx] = highlight_color

        return highlights, highlight_colors


class RotatableBondsHighlighter:
    """A class responsible for highlighting of rotatable bonds in a molecule.

    Parameters:
    -----------
    molecule : RDKit.Chem.rdchem.Mol
        An RDKit molecule object representing the molecule that should be highlighted.
    checkbox : Checkbox
        A checkbox widget indicating whether to highlight rotatable bonds.

    Returns:
    --------
    highlights : dict
        A dictionary with keys "atoms" and "bonds", each containing lists of atom indices and bond indices
        that should be highlighted, respectively.
    highlight_colors : dict
        A dictionary mapping atom indices and bond indices to their respective highlight colors.
    """

    def __init__(self, molecule: Chem.Mol, checkbox: Dict[str, "Checkbox"]):
        """Initialize the RotatableBondsHighlighter with a molecule and a checkbox for rotatable bonds."""
        self.molecule = molecule
        self.checkbox = checkbox

    def determine_highlights(
        self,
    ) -> Tuple[Dict[str, List[int]], Dict[int, str]]:
        """Determines the atom and bond indices to highlight them based on the presence of rotatable bonds.

        Returns:
        --------
        highlights : dict
            A dictionary with keys "atoms" and "bonds", each containing lists of atom indices and bond indices
            that should be highlighted.
        highlight_colors : dict
            A dictionary mapping atom indices and bond indices to their respective highlight colors.
        """
        highlights = {"atoms": [], "bonds": []}
        highlight_colors = {}

        if self.checkbox.value:
            rot_atom_pairs = self.molecule.GetSubstructMatches(
                RotatableBondSmarts
            )
            for atom_pair in rot_atom_pairs:
                bond = self.molecule.GetBondBetweenAtoms(*atom_pair)
                if bond is not None:
                    highlights["bonds"].append(bond.GetIdx())
                    for atom_id in atom_pair:
                        highlight_colors[atom_id] = (
                            1.0,
                            0.6,
                            0.2,
                        )  # Orange for rotatable bonds

        return highlights, highlight_colors


class PartialChargeHighlighter:
    """A class responsible for highlighting of partial charges in a molecule.

    Parameters:
    -----------
    molecule : RDKit.Chem.rdchem.Mol
        An RDKit molecule object representing the molecule that should be highlighted.
    checkbox : Checkbox
        A checkbox widget indicating whether to highlight the partial charges of the molecule.

    Returns:
    --------
    highlights : dict
        A dictionary with keys "atoms" and "bonds", each containing lists of atom indices and bond indices
        that should be highlighted, respectively.
    highlight_colors : dict
        A dictionary mapping atom indices and bond indices to their respective highlight colors.
    """

    def __init__(self, molecule: Chem.Mol, checkbox: Dict[str, "Checkbox"]):
        """Initialize the PartialChargeHighlighter with a molecule and a checkbox for partial charges."""
        self.molecule = molecule
        self.checkbox = checkbox

    def determine_highlights(
        self,
    ) -> Tuple[Dict[str, List[int]], Dict[int, str]]:
        """Determines the atom and bond indices to highlight them based on the partial charges.

        Returns:
        --------
        highlights : dict
            A dictionary with keys "atoms" and "bonds", each containing lists of atom indices and bond indices
            that should be highlighted.
        highlight_colors : dict
            A dictionary mapping atom indices and bond indices to their respective highlight colors.
        """
        highlights = {"atoms": [], "bonds": []}
        highlight_colors = {}

        if self.checkbox.value:
            AllChem.ComputeGasteigerCharges(self.molecule)
            for atom in self.molecule.GetAtoms():
                atom_idx = atom.GetIdx()
                partial_charge = atom.GetProp("_GasteigerCharge")
                partial_charge = float(partial_charge)

                scaled_charge = (partial_charge + 1.0) / 2.0
                red_component = 0.5 + 0.5 * scaled_charge
                blue_component = 1.0 - 0.5 * scaled_charge
                green_component = 0.5 + 0.5 * (1 - scaled_charge)

                color = (
                    red_component,
                    green_component,
                    blue_component,
                )  # Blue to Red spectrum

                highlight_colors[atom_idx] = color
                highlights["atoms"].append(atom_idx)

        return highlights, highlight_colors


class StereocenterHighlighter:
    """A class responsible for highlighting of stereocenters in a molecule.

    Parameters:
    -----------
    molecule : RDKit.Chem.rdchem.Mol
        An RDKit molecule object representing the molecule that should be highlighted.
    checkbox : Checkbox
        A checkbox widget indicating whether to highlight the stereocenters of the molecule.

    Returns:
    --------
    highlights : dict
        A dictionary with keys "atoms" and "bonds", each containing lists of atom indices and bond indices
        that should be highlighted, respectively.
    highlight_colors : dict
        A dictionary mapping atom indices and bond indices to their respective highlight colors.
    """

    def __init__(self, molecule: Chem.Mol, checkbox: Dict[str, "Checkbox"]):
        """Initialize the StereocenterHighlighter with a molecule and a checkbox for stereocenters."""
        self.molecule = molecule
        self.checkbox = checkbox

    def determine_highlights(
        self,
    ) -> Tuple[Dict[str, List[int]], Dict[int, str]]:
        """Determine the atom and bond indices to highlight based on their stereocenter configuration.

        Returns:
        --------
        highlights : dict
            A dictionary with keys "atoms" and "bonds", each containing lists of atom indices and bond indices
            that should be highlighted.
        highlight_colors : dict
            A dictionary mapping atom indices and bond indices to their respective highlight colors.
        """
        highlights = {"atoms": [], "bonds": []}
        highlight_colors = {}

        if self.checkbox.value:
            chiral_centers = Chem.FindMolChiralCenters(
                self.molecule, includeUnassigned=True
            )
            for chiral_atom in chiral_centers:
                atom_idx = chiral_atom[0]
                chirality = chiral_atom[1]
                highlights["atoms"].append(atom_idx)
                if chirality == "S":
                    highlight_colors[atom_idx] = (
                        1.0,
                        0.0,
                        1.0,
                    )  # Magenta for S stereocenters
                elif chirality == "R":
                    highlight_colors[atom_idx] = (
                        0.25,
                        0.88,
                        0.82,
                    )  # Red for R stereocenters
                else:
                    highlight_colors[atom_idx] = (
                        0.5,
                        0.5,
                        0.5,
                    )  # Grey for unassigned stereocenters

        return highlights, highlight_colors


class MurckoScaffoldHighlighter:
    """A class responsible for highlighting of the Murcko scaffold in a molecule.

    Parameters:
    -----------
    molecule : RDKit.Chem.rdchem.Mol
        An RDKit molecule object representing the molecule that should be highlighted.
    checkbox : Checkbox
        A checkbox widget indicating whether to highlight the Murcko scaffold of the molecule.

    Returns:
    --------
    highlights : dict
        A dictionary with keys "atoms" and "bonds", each containing lists of atom indices and bond indices
        that should be highlighted, respectively.
    highlight_colors : dict
        A dictionary mapping atom indices and bond indices to their respective highlight colors.
    """

    def __init__(self, molecule: Chem.Mol, checkbox: Dict[str, "Checkbox"]):
        """Initialize the MurckoScaffoldHighlighter with a molecule and a checkbox for Murcko scaffolds."""
        self.molecule = molecule
        self.checkbox = checkbox

    def determine_highlights(
        self,
    ) -> Tuple[Dict[str, List[int]], Dict[int, str]]:
        """Determines the atom and bond indices to highlight them based on the presence of Murcko scaffolds.

        Returns:
        --------
        highlights : dict
            A dictionary with keys "atoms" and "bonds", each containing lists of atom indices and bond indices
            that should be highlighted.
        highlight_colors : dict
            A dictionary mapping atom indices and bond indices to their respective highlight colors.
        """
        highlights = {"atoms": [], "bonds": []}
        highlight_colors = {}

        if self.checkbox.value:
            scaffold = MurckoScaffold.GetScaffoldForMol(self.molecule)
            scaffold_atoms = scaffold.GetAtoms()
            scaffold_atom_indices = [atom.GetIdx() for atom in scaffold_atoms]
            scaffold_bonds = scaffold.GetBonds()
            scaffold_bond_indices = [bond.GetIdx() for bond in scaffold_bonds]

            highlights["atoms"].extend(scaffold_atom_indices)
            highlights["bonds"].extend(scaffold_bond_indices)

            scaffold_color = (1, 0.8, 0.8)  # Light lila for scaffold
            for atom_idx in scaffold_atom_indices:
                highlight_colors[atom_idx] = scaffold_color
            for bond_idx in scaffold_bond_indices:
                highlight_colors[bond_idx] = scaffold_color

        return highlights, highlight_colors


class MoleculeDrawer:
    """A class for drawing a molecule with various highlighted features.

    Parameters:
    -----------
    molecule : RDKit.Chem.Mol
        An RDKit molecule object representing the molecule to be drawn.
    pharmacophore_checkboxes : dict
        A dictionary of checkbox widgets indicating which pharmacophore features to highlight.
        Keys are pharmacophore feature names (e.g., "Aromatic"), and values are the corresponding checkbox widgets.
    functional_groups_checkboxes : dict
        A dictionary of checkbox widgets indicating which functional groups to highlight.
        Keys are functional group names (e.g., "Alcohol"), and values are the corresponding checkbox widgets.
    rotatable_bonds_checkbox : Checkbox
        A checkbox widget indicating whether to highlight rotatable bonds of the molecule.
    partial_charges_checkbox : Checkbox
        A checkbox widget indicating whether to display partial charges of the molecule.
    partial_charges_heatmap_checkbox : Checkbox
        A checkbox widget indicating whether to use a heatmap for partial charges of the molecule.
    stereocenters_checkbox : Checkbox
        A checkbox widget indicating whether to highlight the stereocenters of the molecule.
    murcko_scaffold_checkbox : Checkbox
        A checkbox widget indicating whether to highlight the Murcko scaffold of the molecule.
    factory : ChemicalFeatures.MolChemicalFeatureFactory
        A factory object used to generate pharmacophore features for the molecule.

    Returns:
    --------
    svg : str
        An SVG string representing the drawn molecule with highlighted features.
    """

    def __init__(
        self,
        molecule: Chem.Mol,
        pharmacophore_checkboxes: Dict[str, "Checkbox"],
        functional_groups_checkboxes: Dict[str, "Checkbox"],
        rotatable_bonds_checkbox: "Checkbox",
        partial_charges_checkbox: "Checkbox",
        partial_charges_heatmap_checkbox: "Checkbox",
        stereocenters_checkbox: "Checkbox",
        murcko_scaffold_checkbox: "Checkbox",
        factory: ChemicalFeatures.MolChemicalFeatureFactory,
    ):
        """Initialize the MoleculeDrawer with a molecule and various highlighting checkboxes and features."""
        self.molecule = molecule
        self.pharmacophore_checkboxes = pharmacophore_checkboxes
        self.functional_groups_checkboxes = functional_groups_checkboxes
        self.rotatable_bonds_checkbox = rotatable_bonds_checkbox
        self.partial_charges_checkbox = partial_charges_checkbox
        self.partial_charges_heatmap_checkbox = partial_charges_heatmap_checkbox
        self.stereocenters_checkbox = stereocenters_checkbox
        self.murcko_scaffold_checkbox = murcko_scaffold_checkbox
        self.factory = factory

    def draw_molecule(
        self, show_atom_indices: bool, width: int, height: int
    ) -> str:
        """Draws the molecule with the highlighted features.

        Parameters:
        -----------
        show_atom_indices : bool
            Whether to show atom indices on the drawn molecule.
        width : int
            Width of the drawing canvas.
        height : int
            Height of the drawing canvas.

        Returns:
        --------
        svg : str
            An SVG string representing the drawn molecule with highlighted features.
        """
        highlights: Dict[str, List[int]] = {"atoms": [], "bonds": []}
        highlight_colors: Dict[int, str] = {}

        # Pharmacophore Highlighter
        pharmacophore_highlighter = PharmacophoreHighlighter(
            self.molecule, self.pharmacophore_checkboxes, self.factory
        )
        pharma_highlights, pharma_colors = (
            pharmacophore_highlighter.determine_highlights()
        )
        highlights["atoms"].extend(pharma_highlights["atoms"])
        highlights["bonds"].extend(pharma_highlights["bonds"])
        highlight_colors.update(pharma_colors)

        # Functional Groups
        fg_highlighter = FunctionalGroupHighlighter(
            self.molecule, self.functional_groups_checkboxes
        )
        fg_highlights, fg_colors = fg_highlighter.determine_highlights()
        highlights["atoms"].extend(fg_highlights["atoms"])
        highlights["bonds"].extend(fg_highlights["bonds"])
        highlight_colors.update(fg_colors)

        # Rotatable Bonds
        rot_bonds_highlighter = RotatableBondsHighlighter(
            self.molecule, self.rotatable_bonds_checkbox
        )
        rot_bonds_highlights, rot_bonds_colors = (
            rot_bonds_highlighter.determine_highlights()
        )
        highlights["atoms"].extend(rot_bonds_highlights["atoms"])
        highlights["bonds"].extend(rot_bonds_highlights["bonds"])
        highlight_colors.update(rot_bonds_colors)

        # Partial Charges
        partial_charge_highlighter = PartialChargeHighlighter(
            self.molecule, self.partial_charges_heatmap_checkbox
        )
        pc_highlights, pc_colors = (
            partial_charge_highlighter.determine_highlights()
        )
        highlights["atoms"].extend(pc_highlights["atoms"])
        highlight_colors.update(pc_colors)

        # Stereocenters
        stereo_highlighter = StereocenterHighlighter(
            self.molecule, self.stereocenters_checkbox
        )
        stereo_highlights, stereo_colors = (
            stereo_highlighter.determine_highlights()
        )
        highlights["atoms"].extend(stereo_highlights["atoms"])
        highlight_colors.update(stereo_colors)

        # Murcko Scaffold
        scaffold_highlighter = MurckoScaffoldHighlighter(
            self.molecule, self.murcko_scaffold_checkbox
        )
        scaffold_highlights, scaffold_colors = (
            scaffold_highlighter.determine_highlights()
        )
        highlights["atoms"].extend(scaffold_highlights["atoms"])
        highlights["bonds"].extend(scaffold_highlights["bonds"])
        highlight_colors.update(scaffold_colors)

        d = rdMolDraw2D.MolDraw2DSVG(width, height)
        draw_options = d.drawOptions()

        if self.partial_charges_checkbox.value:
            AllChem.ComputeGasteigerCharges(self.molecule)
            for atom in self.molecule.GetAtoms():
                atom_idx = atom.GetIdx()
                partial_charge = atom.GetProp("_GasteigerCharge")
                draw_options.atomLabels[atom_idx] = (
                    f"{float(partial_charge):.2f}"
                )

        d.drawOptions().addAtomIndices = show_atom_indices
        d.drawOptions().addStereoAnnotation = True
        rdMolDraw2D.PrepareAndDrawMolecule(
            d,
            self.molecule,
            highlightAtoms=highlights["atoms"],
            highlightBonds=highlights["bonds"],
            highlightAtomColors=highlight_colors,
            highlightBondColors=highlight_colors,
        )
        d.FinishDrawing()
        svg = d.GetDrawingText()
        return svg
