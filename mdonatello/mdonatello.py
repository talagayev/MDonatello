import MDAnalysis as mda
from ipywidgets import interact, Layout, VBox, HTML, Dropdown, Button, Checkbox
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import AllChem, Descriptors
from IPython.display import display, clear_output
from io import BytesIO
import base64


