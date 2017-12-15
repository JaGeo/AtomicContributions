[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://github.com/JaGeo/AtomicContributions/blob/master/LICENSE) [![DOI](https://zenodo.org/badge/110969681.svg)](https://zenodo.org/badge/latestdoi/110969681)

# AtomicContributions
This python package can visualize the contribution of each atom to the phonon modes at the Gamma point. To do so, you need [```Phonopy```](https://github.com/atztogo/phonopy) and [```VASP```](https://www.vasp.at/). 
<hr></hr>

What to cite
------------
Please cite this repository: 
Janine George, & Richard Dronskowski. (2017, December 12). AtomicContributions (Version 1.1). Zenodo. [http://doi.org/10.5281/zenodo.1112391](http://doi.org/10.5281/zenodo.1112391) ([Bibtex](https://zenodo.org/record/1112391/export/hx)). 

Of course, also [```VASP```](https://www.vasp.at/) and [```Phonopy```](https://github.com/atztogo/phonopy).

Intallation
-----------
To use this package you need to install [```Phonopy```](https://github.com/atztogo/phonopy) correctly. Furthermore, ```numpy``` and ```matplotlib``` are required. Also, the python path should be exported correctly. 

How to
--------
1. Perform a phonon calculation with Phonopy and VASP (finite displacements or DFPT) ([More information on this procedure](https://atztogo.github.io/phonopy/procedure.html))
2. Generate the ```FORCE_SETS``` or ```FORCE_CONSTANTS``` file
3. If neccessary calculate the BORN charges ([More information on this procedure](https://atztogo.github.io/phonopy/procedure.html)) and the ```BORN``` file
4. Download this repository, export the Python path correctly
5. Copy an example script, adapt the names of the files, the supercell size (the one you used for the phonon calculation!) and include the rotational matrix to arrive at the primitive cell if necessary.  
6. Run the script

Result
------

![all text](https://github.com/JaGeo/AtomicContributions/blob/master/Doc/allmodes.png)


Todo
--------
1. Other functionalities
2. Include error handling. Not included so far. Thus, be careful.

Information about the Author
--------

- J. George (RWTH Aachen University)
- PI during the development of the code: [R. Dronskowski, RWTH Aachen University](http://www.ssc.rwth-aachen.de/)

