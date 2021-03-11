# pythonpolarisation

Python-based model to calculate electronic polarisation in arrays of organic molecules. 

Developed by Sheridan Few based upon a framework built by James Kirkpatrick.

As used in:

 Few, S., Chia, C., Teo, D., Kirkpatrick, J. & Nelson, J. The impact of chemical structure and molecular packing on the electronic polarisation of fullerene arrays. Phys. Chem. Chem. Phys. 19, (2017). http://doi.org/10.1039/c7cp00317j

 Few, S. Theoretical studies of charge transfer excitations, absorption, and polarisation in organic photovoltaic materials, 2015 https://spiral.imperial.ac.uk/handle/10044/1/33284 (my thesis)

 A draft paper calculating electronic reorganisation energies and transfer integrals in a set of molecules first reported by Minder et al. (https://onlinelibrary.wiley.com/doi/epdf/10.1002/adma.201103960), yet to see the light of day...

### Overview ###

This model has taken much longer than intended to upload. Unforunately I've lost familarity with details of the process, but am uploading the model now in the hope that it might be useful to others considering similar questions.

The following may be helpful:

BasicElements/ contains scripts defining key attributes of the model

BasicElements/tests/ and tests/ - Directories containing a number of small examples for initiating single atoms/ions and calculating interactions between pairs of atoms/ions/molecules which may be helpful in getting to grips with the model

Crystal_Building_Instructions_0515.odt (or txt file of same name) - an overview of how to build a crystal using this model (some creativity may be required)

Crystals/Morpurgo_all_sp_Reorgs/Jobs/Pc/Pc_anion_neut_inner1_outer0/Pc_anion_neut_inner1_outer0.py - a script for calculating the electronic interaction energy of a single pentacene anion with its neighbouring molecules (one unit cell in each direction). Atoms are kept in equilibrium crystal positions and rearrangement of atomic positions in response to the charge is not considered. "Electronic reorganisation energy" associated with the change in electronic interaction when the anion moves without changes to electronic polarisation of surroundings is also calculated in this script, but this may be a less useful concept. 

Crystals/video[...].gnp - I remember having fun using scripts like these to create pngs which I turned into animated gifs showing how atomic dipoles change as you separate a pair of molecules. Perhaps something here could be useful for visualisations.

There is a lot that could be tidied amongst these directories - I have left a lot in place, but much could be safely ignored/deleted