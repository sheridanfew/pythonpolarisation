# pythonpolarisation
Python-based model to calculate electronic polarisation in arrays of organic molecules.

As used in:

 Few, S., Chia, C., Teo, D., Kirkpatrick, J. & Nelson, J. The impact of chemical structure and molecular packing on the electronic polarisation of fullerene arrays. Phys. Chem. Chem. Phys. 19, (2017). http://doi.org/10.1039/c7cp00317j

 Few, S. Theoretical studies of charge transfer excitations, absorption, and polarisation in organic photovoltaic materials, 2015 https://spiral.imperial.ac.uk/handle/10044/1/33284 (my thesis)

 A draft paper calculating electronic reorganisation energies and transfer integrals in a set of molecules first reported by Minder et al. (https://onlinelibrary.wiley.com/doi/epdf/10.1002/adma.201103960) 

### Overview ###

I am a bit embarrassed how long it has taken to upload this model, but hopefully better late than never! Unfortunately, much of the workings have been lost in the mists of time, but the following may be helpful:

tests/ - Contains a number of small examples for initiating single atoms/ions and calculating interactions between pairs of atoms/ions/molecules which may be helpful in getting to grips with the model

Crystal_Building_Instructions_0515.odt (or txt file of same name) - a detailed overview of how to build a crystal using this model (I seem to remember there were some additional steps I wanted to add, some creativity may be required)

Crystals/Morpurgo_all_sp_Reorgs_qsplit_Molscreen_largepiatoms_sq/Jobs/C60/C60_anion_neut_inner1_outer0/C60_anion_neut_inner1_outer0.py - a script for calculating the electronic interaction energy of a single C60 anion with its neighbouring molecules (one unit cell in each direction - NB. atoms are kept in equilibrium crystal positions rearrangement of atomic positions in response to the charge is not considered)

Crystals/video[...].gnp - I remember having fun using scripts like these to create pngs which I turned into animated gifs showing how atomic dipoles change as you separate a pair of molecules. Perhaps something in these could be useful for visualisations.

There is a lot that could be tidied in these directories, but I have left most of it in place as I may no longer be the best judge of what is and isn't useful here.
