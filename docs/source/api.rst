API Documentation
==============================

Here is the API Documentation.

The following API represents the MoleculeVisualizer, which is responsible for the display of the
molecule in IPywidgets.

.. toctree::
   :maxdepth: 2
   :caption: API MoleculeVisualizer:

   modules/mdonatello

The following API represents the Drawer, which is responsible for the highlighting and drawing of the
molecule with RDKit, which then serves as input for the display in MoleculeVisualizer.

.. toctree::
   :maxdepth: 2
   :caption: API Modules:

   modules/drawer
   
The following API represents the Mapper, which is responsible for the mapping of the
molecule according to the pharmacophores and functional groups, including the selection of the color code
that should be used for the highlighting, which are then given to the Drawer.
   
.. toctree::
   :maxdepth: 2
   :caption: API Modules:

   modules/mapper
   
The following API represents the Properties, which is responsible for the calculation of certain properties
of the molecule, which are then returend to the display of MoleculeVisualizer
   
.. toctree::
   :maxdepth: 2
   :caption: API Modules:

   modules/properties
