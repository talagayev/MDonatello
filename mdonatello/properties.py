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
        """
        Initializes the Property class with a molecule.

        Parameters:
        -----------
        mol : Chem.Mol
            The RDKit Mol object of the molecule.
        """
        self.mol = mol

    @cached_property
    def property_value(self) -> Union[float, str]:
        """
        Abstract method to be implemented by subclasses to calculate the property value.

        Returns:
        --------
        Union[float, str]
            The calculated property value.
        """
        raise NotImplementedError("Subclasses should implement this.")

    def __repr__(self) -> str:
        """
        Returns the HTML-formatted string representation of the property name and value.

        Returns:
        --------
        repr_str : str
            HTML-formatted string representation of the property name and value.
        """
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
        """
        Calculates the molecular weight of the molecule.

        Returns:
        --------
        float
            The molecular weight of the molecule.
        """
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
        """
        Calculates the LogP value of the molecule.

        Returns:
        --------
        float
            The LogP value of the molecule.
        """
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
        """
        Calculates the TPSA (Topological Polar Surface Area) value of the molecule.

        Returns:
        --------
        float
            The TPSA value of the molecule.
        """
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
        """
        Calculates the number of rotatable bonds of the molecule.

        Returns:
        --------
        float
            The number of rotatable bonds of the molecule.
        """
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
        """
        Calculates the number of hydrogen bond acceptors of the molecule.

        Returns:
        --------
        float
            The number of hydrogen bond acceptors of the molecule.
        """
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
        """
        Calculates the number of hydrogen bond donors of the molecule.

        Returns:
        --------
        float
            The number of hydrogen bond donors of the molecule.
        """
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
        """
        Calculates the number of stereocenters of the molecule.

        Returns:
        --------
        float
            The number of stereocenters of the molecule.
        """
        return len(Chem.FindMolChiralCenters(self.mol))
