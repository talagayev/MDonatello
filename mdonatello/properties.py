from typing import Union
from functools import cached_property
from rdkit import Chem, RDConfig
from rdkit.Chem import (
    Descriptors,
    Lipinski,
)
from rdkit.Chem.Lipinski import RotatableBondSmarts
from mdonatello.mapper import FunctionalGroupHandler


class Property:
    """Base class for representing a property of a molecule.

    Parameters:
    -----------
    mol : Chem.Mol
        RDKit Mol object for the molecule for which the property is being calculated.

    Returns:
    --------
    repr_str : str
        HTML-formatted string representation of the property name and value.
    """

    name: str = "property"
    values_format: bool = True

    def __init__(self, mol: Chem.Mol):
        self.mol = mol

    @cached_property
    def property_value(self) -> Union[float, str]:
        raise NotImplementedError("Subclasses should implement this.")

    def __repr__(self) -> str:
        value = (
            f"{self.property_value:.2f}"
            if self.values_format
            else str(self.property_value)
        )
        repr_str = repr(f"{self.name}: <b>{value}</b>")
        return repr_str.replace("'", "")


class MolecularWeight(Property):
    """A class for calculating the molecular weight of a molecule.

    Parameters:
    -----------
    mol : Chem.Mol
        The RDKit Mol object of the molecule for which the molecular weight is being calculated.

    Returns:
    --------
    float
        The molecular weight of the molecule.
    """

    name = "Molecular Weight"
    values_format = True

    @cached_property
    def property_value(self) -> float:
        return Descriptors.MolWt(self.mol)


class LogP(Property):
    """A class for calculating the LogP value of a molecule.

    Parameters:
    -----------
    mol : Chem.Mol
        The RDKit Mol object of the molecule for which the LogP value is being calculated.

    Returns:
    --------
    float
        The LogP value of the molecule.
    """

    name = "LogP"
    values_format = True

    @cached_property
    def property_value(self) -> float:
        return Descriptors.MolLogP(self.mol)


class TPSA(Property):
    """A class for calculating the TPSA value of a molecule.

    Parameters:
    -----------
    mol : Chem.Mol
        The RDKit Mol object of the molecule for which the TPSA value is being calculated.

    Returns:
    --------
    float
        The TPSA value of the molecule.
    """

    name = "TPSA"
    values_format = True

    @cached_property
    def property_value(self):
        return Descriptors.TPSA(self.mol)


class RotatableBonds(Property):
    """A class for calculating the number of the rotatable bonds of a molecule.

    Parameters:
    -----------
    mol : Chem.Mol
        The RDKit Mol object of the molecule for which the number of rotatable bonds is being calculated.

    Returns:
    --------
    float
        The number of rotatable bonds of the molecule.
    """

    name = "Rotatable Bonds"
    values_format = False

    @cached_property
    def property_value(self) -> float:
        return Descriptors.NumRotatableBonds(self.mol)


class HydrogenBondAcceptors(Property):
    """A class for calculating the number of the hydrogen bond acceptors of a molecule.

    Parameters:
    -----------
    mol : Chem.Mol
        The RDKit Mol object of the molecule for which the number of hydrogen bond acceptors is being calculated.

    Returns:
    --------
    float
        The number of hydrogen bond acceptors of the molecule.
    """

    name = "Hydrogen Bond Acceptors"
    values_format = False

    @cached_property
    def property_value(self) -> float:
        return Lipinski.NumHAcceptors(self.mol)


class HydrogenBondDonors(Property):
    """A class for calculating the number of the hydrogen bond donors of a molecule.

    Parameters:
    -----------
    mol : Chem.Mol
        The RDKit Mol object of the molecule for which the number of hydrogen bond donors is being calculated.

    Returns:
    --------
    float
        The number of hydrogen bond acceptors of the molecule.
    """

    name = "Hydrogen Bond Donors"
    values_format = False

    @cached_property
    def property_value(self) -> float:
        return Lipinski.NumHDonors(self.mol)


class Stereocenters(Property):
    """A class for calculating the number of stereocenters of a molecule.

    Parameters:
    -----------
    mol : Chem.Mol
        The RDKit Mol object of the molecule for which the number of stereocenters is being calculated.

    Returns:
    --------
    float
        The number of stereocenters of the molecule.
    """

    name = "Stereocenters"
    values_format = False

    @cached_property
    def property_value(self) -> float:
        # Using the provided mol directly assuming it's an RDKit molecule object
        return len(Chem.FindMolChiralCenters(self.mol))
