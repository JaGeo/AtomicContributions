from ContributionsOfAtomsToModes  import AtomicContributionToModes 
from phonopy.units import VaspToCm


#Initialize AtomicContributionToModes 
test=AtomicContributionToModes(PoscarName='POSCAR',ForceConstants=False,ForceFileName='FORCE_SETS',supercell=[[3, 0, 0],[0, 3, 0], [0, 0, 4]],primitive=[[1, 0, 0],[0, 1, 0], [0, 0, 1]],)


#write a file with atomic contributions to each mode
test.write_file();



#atom numbers start at 1, not at 0
atomgroups=[[1,2],[3,4],[5,6,7,8],[9,10,11,12,13,14,15,16]] #this list groups the atoms for the plots
colorofgroups=['black','red','green','grey'] #this list gives the colors for the grouped atoms
legendforgroups=['C','O','N','H'] #this list gives the legend for the grouped atoms

#Plot with all modes
test.plot(atomgroups=atomgroups,colorofgroups=colorofgroups,legendforgroups=legendforgroups,filename='allmodes.eps',transmodes=False)

#Plot with modes that have the irreducible representation B2 and E (leaves out translational modes by default)
#Have a look at https://github.com/atztogo/phonopy/blob/master/phonopy/phonon/irreps.py for the strings you need to represent the irreducible representations
test.plot_irred(atomgroups=atomgroups,colorofgroups=colorofgroups,legendforgroups=legendforgroups,irreps=['B2','E'],filename='IRactivemodes.eps')


