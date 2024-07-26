from typing import Dict, Tuple
from rdkit import Chem


class FunctionalGroupHandler:
    """A class that is responsible for handling functional groups in molecules.

    Parameters:
    -----------
    mol : RDKit.Chem.Mol
        An RDKit molecule object that represents the molecule in which functional groups are identified.
    fg : str
        The functional group for which the color is requested.

    Returns:
    --------
    fg_counts : dict
        A dictionary with functional group names as keys and lists of atom indices as values,
        representing the atoms involved in each functional group.
    color : Tuple[float, float, float]
        An RGB tuple representing the color assigned to the specified functional group.
    """

    @staticmethod
    def calculate_functional_groups(mol):
        """
        Identifies, calculates and counts the functional groups present in the given RDKit molecule.

        Parameters:
        -----------
        mol : RDKit.Chem.Mol
            An RDKit molecule object representing the molecule in which functional groups are identified.

        Returns:
        --------
        fg_counts : dict
            A dictionary with functional group names as keys and lists of atom indices as values,
            representing the atoms involved in each functional group.
        """
        functional_groups = {
            "Hydroxyl group (-OH)": "[OX2H]",
            "Primary amine (-NH2)": "[NX3H2]",
            "Primary ammonium (-NH3+)": "[+NX4;H3]",
            "Secondary amine (-NH-)": "[NX3H][#6]",
            "Tertiary amine (-N<)": "[NX3;H0]([#6])[#6]",
            "Carboxyl group (-COOH)": "C(=O)[OX2H1]",
            "Ester (-COOR)": "C(=O)[OX2H0][#6]",
            "Amide (-CON-)": "C(=O)[NX3]",
            "Aldehyde (-CHO)": "[CX3H1](=O)[#6]",
            "Ketone (C=O)": "[CX3](=O)[#6]",
            "Ether (R-O-R)": "[#6][OX2][#6]",
            "Thiocarbonyl group (C=S)": "C(=S)",
            "Imine group (-C=N-)": "[CX3](=N)",
            "Hydroxylamine group (-N(OH))": "[NX3][OX2H]",
            "Thiol group (-SH)": "[SX2H]",
            "Azide group (-N3)": "N=[NX1]=[NX1]",
            "Furan ring": "c1occc1",
            "Guanidine group (-C(=NH)(N)(NH2))": "C(=N)(N)[NH2]",
            "Isothiocyanate (-N=C=S)": "[NX2]=C=[SX2]",
            "Isocyanate (-N=C=O)": "[NX2]=C=[OX1]",
            "Lactone (C=O-O)": "[CX3](=O)[OX2][CX3](=O)",
            "Lactam (C=O-N)": "[CX3](=O)[NX3][CX3](=O)",
            "Methoxy group (-OCH3)": "[OX2][CH3]",
            "Nitro group (-NO2)": "[NX3](=O)=O",
            "Nitroso group (-NO)": "[NX2]=O",
            "Oxazole ring": "c1noccc1",
            "Oxime group (-C=N-OH)": "[CX3](=N[OX2H])",
            "Epoxide": "C1CO1",
            "Nitrile": "C#N",
            "Sulfone": "S(=O)(=O)([#6])([#6])",
            "Sulfonamide": "S(=O)(=O)([#6])N",
            "Sulfide": "[SX2]",
            "Urea": "C(=O)(N)(N)",
            "Phosphoric Ester": "P(=O)(O)([OX2H0;R1])",
            "Phosphoric Acid": "P(=O)(O)(O)",
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
        """
        Determines the color associated with the given functional group based on its chemical structure.

        Parameters:
        -----------
        fg : str
            The functional group for which the color is required.

        Returns:
        --------
        color : Tuple[float, float, float]
            An RGB tuple representing the color assigned to the specified functional group.
        """
        parts = fg.split("(")
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


class PharmacophoreColorMapper:
    """A class for mapping pharmacophore features according to their respective colors.

    Parameters:
    -----------
    family : str
        The pharmacophore feature family for which the color is requested.

    Returns:
    --------
    color : Tuple[float, float, float]
        An RGB tuple representing the color assigned to the specified pharmacophore feature family.
    """

    @staticmethod
    def get_color_for_pharmacophore(family: str) -> Tuple[float, float, float]:
        """
        Determines the color associated with the given pharmacophore feature family.

        Parameters:
        -----------
        family : str
            The pharmacophore feature family for which the color is required.

        Returns:
        --------
        color : Tuple[float, float, float]
            An RGB tuple representing the color assigned to the specified pharmacophore feature.
        """
        color_map = {
            "Donor": (0.0, 1.0, 0.0),  # Green
            "Acceptor": (1.0, 0.7, 0.7),  # Rosa
            "Hydrophobe": (1.0, 1.0, 0.0),  # Yellow
            "PosIonizable": (0.0, 1.0, 1.0),  # Turquoise
            "NegIonizable": (1.0, 0.0, 1.0),  # Pink
            "Aromatic": (0.5, 0.5, 1.0),  # Light Blue
            "LumpedHydrophobe": (1.0, 0.5, 0.0),  # Orange
        }
        return color_map.get(
            family, (0.5, 0.5, 0.5)
        )  # Default to grey if not specified
