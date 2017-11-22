#Copyright (C) 2017 Janine George

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from phonopy import Phonopy
from phonopy.interface.vasp import read_vasp
from phonopy.file_IO import parse_FORCE_SETS, parse_BORN, write_FORCE_CONSTANTS, parse_FORCE_CONSTANTS
from phonopy.structure.atoms import PhonopyAtoms
from phonopy.units import VaspToCm, VaspToTHz,VaspToEv
from phonopy.phonon.irreps import IrReps
from phonopy.phonon.degeneracy import degenerate_sets as get_degenerate_sets

import os



class AtomicContributionToModes:

	def __init__(self,PoscarName='POSCAR',ForceConstants=False,ForceFileName='FORCE_SETS',BornFileName='BORN',supercell=[[1, 0, 0],[0, 1, 0], [0, 0, 1]],nac=False,symprec=1e-5,masses=[],primitive=[[1, 0, 0],[0, 1, 0], [0, 0, 1]],degeneracy_tolerance=1e-4,factor=VaspToCm,q=[0,0,0]):
		"""Class that calculates contributions of each atom to the phonon modes at Gamma
			Args:
			PoscarNamse (str): name of the POSCAR that was used for the phonon calculation
			BornFileName (str): name of the file with BORN charges (formatted with outcar-born)
			ForceConstants (boolean): If True, ForceConstants are read in. If False, forces are read in.
			ForceFileName (str): name of the file including force constants or forces
			supercell (list of lists): reads in supercell
			nac (boolean): If true, NAC is applied.
			symprec (float): contains symprec tag as used in Phonopy
			masses (list): Masses in this list are used instead of the ones prepared in Phonopy. Useful for isotopes.
			primitive (list of lists): contains rotational matrix to arrive at primitive cell
			factor (float): VaspToCm or VaspToTHz or VaspToEv
			q (list of int): q point for the plot. So far only Gamma works
			
		"""
		
		self.__unitcell =read_vasp(PoscarName)
		self.__supercell=supercell
		self.__phonon= Phonopy(self.__unitcell,supercell_matrix=self.__supercell,primitive_matrix=primitive,factor=factor,symprec=symprec)
	        self.__natoms=self.__phonon.get_primitive().get_number_of_atoms()	
		self.__symbols= self.__phonon.get_primitive().get_chemical_symbols()
		self.__factor=factor
		#If different masses are supplied
		if masses: 
			self.__phonon.set_masses(masses)
                self.__masses=self.__phonon.get_primitive().get_masses() 
	
		#Forces or Force Constants
		if not ForceConstants:
			self.__set_ForcesSets(filename=ForceFileName,phonon=self.__phonon)
		
		if ForceConstants:
			self.__ForceConstants(filename=ForceFileName,phonon=self.__phonon)
	

		#Apply NAC Correction
                if nac:
                   	BORN_file = parse_BORN(self.__phonon.get_primitive(),filename=BornFileName)
                	self.__BORN_CHARGES=BORN_file['born']
			self.__phonon.set_nac_params(BORN_file)
		

		#frequencies and eigenvectors at Gamma
		self.__frequencies,self.__eigvecs=self.__phonon.get_frequencies_with_eigenvectors(q)
	


		self.__NumberOfBands=len(self.__frequencies)

		#Nicer format of the eigenvector file
		self.__FormatEigenvectors()

		#Get dipole approximation of the intensitiess
		self.__set_Contributions()
		
		#irrepsobject
		self.__set_IRLabels(phonon=self.__phonon,degeneracy_tolerance=degeneracy_tolerance,factor=factor,q=q)
		
		
		
		
	def __set_ForcesSets(self,filename,phonon):
		"""
		sets forces

		"""

		force_sets = parse_FORCE_SETS(filename=filename)
		phonon.set_displacement_dataset(force_sets)
		phonon.produce_force_constants()
		


	def __set_ForceConstants(self,filename,phonon):
		"""
		sets force constants
		"""
		force_constants = parse_FORCE_CONSTANTS(filename=filename)
		phonon.set_force_constants(force_constants)	

	def __set_IRLabels(self,phonon,degeneracy_tolerance,factor,q):
		"""
		sets list of irreducible labels and list of frequencies without degeneracy 
		"""
		phonon.set_dynamical_matrix()
		self.__Irrep=IrReps(dynamical_matrix=phonon._dynamical_matrix,q=q,is_little_cogroup=False,nac_q_direction=None,degeneracy_tolerance=degeneracy_tolerance,factor=factor)
		self.__Irrep.run()
		self.__IRLabels=self.__Irrep._get_ir_labels()
		
		self.__ListOfModesWithDegeneracy=self.__Irrep._get_degenerate_sets()
		self.__freqlist={}
		for band in range(len(self.__ListOfModesWithDegeneracy)):
			self.__freqlist[band]=self.__ListOfModesWithDegeneracy[band][0]


	def __FormatEigenvectors(self):
		"""
		Formats eigenvectors to a dictionary: the first argument is the number of bands, the second the number of atoms, the third the Cartesian coordinate
		"""

		self.__EigFormat = {}
		for alpha in  range(self.__NumberOfBands):
			laufer=0
    	 	   	for beta in range(self.__natoms):
        			for xyz in range(0,3):
      	    	     	 		self.__EigFormat[beta,alpha,xyz]=self.__eigvecs[laufer][alpha]
                       			laufer=laufer+1

	def __Eigenvector(self, atom, band, xoryorz ):
		"""
		Gives a certain eigenvector corresponding to one specific atom, band and Cartesian coordinate

		args: 
			atom (int) : number of the atoms (same order as in POSCAR)
			band (int) : number of the frequency (ordered by energy)
			xoryorz (int): Cartesian coordinate of the eigenvector			
			
			
		"""

		return np.real(self.__EigFormat[atom,band,xoryorz])

	def __massEig(self,atom,band,xoryorz):
		"""
		Gives a certain eigenvector divided by sqrt(mass of the atom) corresponding to one specific atom, band and Cartesian coordinate

		args: 
			atom (int) : number of the atoms (same order as in POSCAR)
			band (int) : number of the frequency (ordered by energy)
			xoryorz (int): Cartesian coordinate of the eigenvector	
			

		"""

		return self.__Eigenvector(atom,band,xoryorz)/np.sqrt(self.__masses[atom])

	def __set_Contributions(self):
		"""
		Calculate contribution of each atom to modes"
		"""
		self.__PercentageAtom = {}
		for freq in range(len(self.__frequencies)):
			for atom in range(self.__natoms):
				sum=0;
				for alpha in range(3):
					sum=sum+abs(self.__Eigenvector(atom,freq,alpha)*self.__Eigenvector(atom,freq,alpha))
				self.__PercentageAtom[freq,atom]=sum
	def __get_Contributions(self,band,atom):
		"""
		Gives contribution of specific atom to modes with certain frequency 
		args:
			band (int): number of the frequency (ordered by energy)
		"""
		return self.__PercentageAtom[band,atom]
		
	
	

	
	
	
	def write_file(self,filename="Contributions.txt"):
		"""
		Writes contributions of each atom in file

		args:
		
			filename (string): filename 
		"""
		file  = open(filename, 'w')
		file.write('Frequency Contributions \n')		
		for freq in range(len(self.__frequencies)):
			file.write('%s ' % (self.__frequencies[freq]))
			for atom in range(self.__natoms):
				file.write('%s ' % (self.__get_Contributions(freq,atom)))
			file.write('\n ')

		file.close()



	def plot(self,atomgroups,colorofgroups,legendforgroups,freqstart=[],freqend=[],freqlist=[],labelsforfreq=[],filename="Plot.eps"):
		"""
		Plots contributions of atoms/several atoms to modes with certain frequencies (freqlist starts at 1 here)	
		
		args:
			atomgroups (list of list of ints): list that groups atoms, atom numbers start at 1
			colorofgroups (list of strings): list that matches a color to each group of atoms
			legendforgroups (list of strings): list that gives a legend for each group of atoms
			freqstart (float): min frequency of plot in cm-1
			freqend (float): max frequency of plot in cm-1
			freqlist (list of int): list of frequencies that will be plotted; if no list is given all frequencies in the range from freqstart to freqend are plotted, list begins at 1
			labelsforfreq (list of strings): list of labels (string) for each frequency
			filename (string): filename for the plot
		"""
		
		p={}
		summe={}
		if labelsforfreq==[]:		
			labelsforfreq=self.__IRLabels



		if freqlist==[]:		
			freqlist=self.__freqlist
			
		else: 	
				
			for freq in range(len(freqlist)):
				freqlist[freq]=freqlist[freq]-1
		
		self._plot(atomgroups=atomgroups,colorofgroups=colorofgroups,legendforgroups=legendforgroups,freqstart=freqstart,freqend=freqend,freqlist=freqlist,labelsforfreq=labelsforfreq,filename=filename)


	def _plot(self,atomgroups,colorofgroups,legendforgroups,freqstart=[],freqend=[],freqlist=[],labelsforfreq=[],filename="Plot.eps"):
		"""
		Plots contributions of atoms/several atoms to modes with certain frequencies (freqlist starts at 0 here)		
		
		args:
			atomgroups (list of list of ints): list that groups atoms, atom numbers start at 1
			colorofgroups (list of strings): list that matches a color to each group of atoms
			legendforgroups (list of strings): list that gives a legend for each group of atoms	
			freqstart (float): min frequency of plot in cm-1
			freqend (float): max frequency of plot in cm-1
			freqlist (list of int): list of frequencies that will be plotted; this freqlist starts at 0
			labelsforfreq (list of strings): list of labels (string) for each frequency
			filename (string): filename for the plot
		"""
		#setting of some parameters in matplotlib: http://matplotlib.org/users/customizing.html
		mpl.rcParams["savefig.directory"] = os.chdir(os.getcwd())
		mpl.rcParams["savefig.format"]='eps'
		
		fig, ax1 = plt.subplots()		
		p={}
		summe={}
				
		for group in range(len(atomgroups)):	
			color1=colorofgroups[group]
			Entry={}		
			for freq in range(len(freqlist)):
				Entry[freq]= 0
			for number in atomgroups[group]:	
				#set the first atom to 0
				atom=int(number)-1
				for freq in range(len(freqlist)):
					
         	       			Entry[freq]= Entry[freq]+ self.__get_Contributions(freqlist[freq],atom)
					if group==0:
						summe[freq]=0	
			
			#plot bar chart
			p[group]=ax1.barh(np.arange(len(freqlist)),Entry.values(),left=summe.values(),color=color1,height=1,label=legendforgroups[group] ) 
			#needed for "left" in the bar chart plot
			for freq in range(len(freqlist)):
				if group==0:
					summe[freq]=Entry[freq]
				else:
					summe[freq]=summe[freq]+Entry[freq]			
		labeling={}
		for freq in range(len(freqlist)):
			labeling[freq]=round(self.__frequencies[freqlist[freq]],1)
		#details for the plot
		plt.rc("font", size=8)
		ax1.set_yticklabels(labeling.values())
		ax1.set_yticks(np.arange(0.5,len(self.__frequencies)+0.5))
		ax2 = ax1.twinx()
		ax2.set_yticklabels(labelsforfreq)
		ax2.set_yticks(np.arange(0.5,len(self.__frequencies)+0.5))
		#start and end of the yrange
		start,end=self.__get_freqbordersforplot(freqstart,freqend,freqlist)
		ax1.set_ylim(start,end)
		ax2.set_ylim(start,end)
		ax1.set_xlim(0.0, 1.0)
		ax1.set_xlabel('Contribution of Atoms to Modes')
		if self.__factor==VaspToCm:
			ax1.set_ylabel('Wavenumber (cm-1)')
		elif self.__factor==VaspToTHz:
			ax1.set_ylabel('Frequency (THz)')
		elif self.__factor==VaspToEv:
			ax1.set_ylabel('Frequency (eV)')
		else:
			ax1.set_ylabel('Frequency')
		ax1.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",mode="expand", borderaxespad=0, ncol=len(atomgroups))
		
		plt.savefig(filename, bbox_inches="tight")
	
		plt.show()
		
	def __get_freqbordersforplot(self,freqstart,freqend,freqlist):
        	if freqstart==[]:
                        start=0.0
                else:
                        for freq in range(len(freqlist)):
                                if self.__frequencies[freqlist[freq]] > freqstart:
                                        start=freq
                                        break
                                else:
                                        start=len(freqlist)

                if freqend==[]:
                        end=len(freqlist)
                else:
                        for freq in range(len(freqlist)-1,0,-1):
                                if self.__frequencies[freqlist[freq]] < freqend:
                                        end=freq+1
                                        break
                                else:
                                        end=len(freqlist)

		return start,end


	def plot_irred(self,atomgroups,colorofgroups,legendforgroups,transmodes=False,irreps=[],filename="Plot.eps",freqstart=[],freqend=[]):
		"""
		Plots contributions of atoms/several atoms to modes with certain irreducible representations (selected by Mulliken symbol)
		args:
			atomgroups (list of list of ints): list that groups atoms, atom numbers start at 1
			colorofgroups (list of strings): list that matches a color to each group of atoms
			legendforgroups (list of strings): list that gives a legend for each group of atoms
			transmodes (boolean): translational modes are included if true
			irreps (list of strings): list that includes the irreducible modes that are plotted
			filename (string): filename for the plot
		"""
		
		
		freqlist=[]
		labelsforfreq=[]
		for band in range(len(self.__freqlist)):
			if self.__IRLabels[band] in irreps:
				if not transmodes:
					if not self.__freqlist[band] in [0,1,2]:				
						freqlist.append(self.__freqlist[band])
						labelsforfreq.append(self.__IRLabels[band])
				else:
					freqlist.append(self.__freqlist[band])
					labelsforfreq.append(self.__IRLabels[band])

		
		self._plot(atomgroups=atomgroups,colorofgroups=colorofgroups,legendforgroups=legendforgroups,filename=filename,freqlist=freqlist,labelsforfreq=labelsforfreq,freqstart=freqstart,freqend=freqend)
 


