Display Configuration
===============

This page details the features present in the display of ``mdonatello``. 

The ``mdonatello`` consists of an IpyWidget, which allows to modify the figure
accoridng to the selection of the molecule in the dropdown menu in addition
to the checking of the separate checkboxes.

Here is an example of the display presented by ``mdonatello``:

.. figure:: /_static/display_example.png
    :figwidth: 700px
    :align: center

The checkboxes present in this display are divided into separate sections, each representing
an important aspect of the molecule visualization.


Properties
---------------------

This section is mainly responsible for the highlighting of different features of the molecule.
These are the checkboxes present in the section:

**Show Phyisochemical Properties**: Displays the following Physiochemical properties of the molecule:

.. code-block:: text

    Molecular Weight: The Molecular Weight of the displayed molecule.
    LogP: The LogP value of the displayed molecule.
    TPSA: The TPSA value of the dispalyed molecule.
    Rotatable Bonds: The number of rotatable bonds present in the displayed molecule.
    Stereocenters: The number of stereocenters present in the displayed molecule.
    
**Show partial charges**: Displays the partial charges of each atom in the displayed molecule.

**Show H-Bond Donors/Acceptors**: Displays the number of Hydrogen bond donors and acceptors present in the molecule.

**Show atom indices**: Displays the atom indices of the displayed molecule.


Atom & Bond Highlighting
---------------------

This section is responsible for highlighting the Atoms and bonds of the displayed molecule according to the selected checkbox.

**Show Rotatable Bonds**: Highlights the rotatable bonds present in the molecule.

**Show partial charge heatmap**: Highlights the atoms of the molecule according to their partial charges.

**Show Functional Groups**: Displays additional functional group checkboxes, which when selected highlight the functional group in the molecule.

**Show Stereocenters**: Highlights the stereocenters present in the molecule according to their configuration.

**Show Murcko Scaffold**: Highlights the murcko scaffold present in the molecule.


Pharmacophores
---------------------

This section is mainly responsible for the highlighting of pharmacophore features of the molecule.
These are the checkboxes present in the section:

**Highlight Donor**: Highlights the hydrogen bond donors present in the molecule.

**Highlight Acceptor**: Highlights the hydrogen bond acceptors present in the molecule.

**Highlight Hydrophobe**: Highlights the hydrophobic features present in the molecule.

**Highlight PosIonizable**: Highlights the positive ionizable features present in the molecule.

**Highlight NegIonizable**: Highlights the negative ionizable features present in the molecule.

**Highlight Aromatic**: Highlights the aromatic features present in the molecule.

**Highlight LumpedHydrophobe**: Highlights the lymphed hydrophobic features present in the molecule.
