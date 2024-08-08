"""
MDonatello
2D small molecule visualization for MDAnalysis
"""

# Add imports here
from importlib.metadata import version
from .mdonatello import MoleculeVisualizer, MolecularWeight, LogP, TPSA, RotatableBonds, HydrogenBondAcceptors, HydrogenBondDonors

__version__ = version("mdonatello")
