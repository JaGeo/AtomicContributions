#Copyright (C) 2017 Janine George

import numpy as np
import matplotlib.pyplot as plt
from phonopy import Phonopy
from phonopy.interface.vasp import read_vasp
from phonopy.file_IO import parse_FORCE_SETS, parse_BORN, write_FORCE_CONSTANTS, parse_FORCE_CONSTANTS
from phonopy.structure.atoms import PhonopyAtoms
from phonopy.units import VaspToCm, VaspToTHz
import os



class AtomicContributionToModes:

	def __init__(self,PoscarName='POSCAR',ForceConstants=False,ForceFileName='FORCE_SETS',BornFileName='BORN',supercell=[[1, 0, 0],[0, 1, 0], [0, 0, 1]],nac=False,symprec=1e-5,masses=[],primitive=[[1, 0, 0],[0, 1, 0], [0, 0, 1]]):
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
			
		"""
	
		self.__unitcell =read_vasp(PoscarName)
		self.__supercell=supercell
		self.__phonon= Phonopy(self.__unitcell,supercell_matrix=self.__supercell,primitive_matrix=primitive,factor=VaspToCm,symprec=symprec)
	        self.__natoms=self.__phonon.get_primitive().get_number_of_atoms()	
		self.__symbols= self.__phonon.get_primitive().get_chemical_symbols()

		#If different masses are supplied
		if masses: 
			self.__phonon.set_masses(masses)
                self.__masses=self.__phonon.get_primitive().get_masses() 
	
		#Forces or Force Constants
		if not ForceConstants:
			self.__force_sets = parse_FORCE_SETS(filename=ForceFileName)
			self.__phonon.set_displacement_dataset(self.__force_sets)
			self.__phonon.produce_force_constants()
		
		if ForceConstants:
			force_constants = parse_FORCE_CONSTANTS(filename=ForceFileName)
			self.__phonon.set_force_constants(force_constants)
	



		#Apply NAC Correction
                if nac:
                   	BORN_file = parse_BORN(self.__phonon.get_primitive(),filename=BornFileName)
                	self.__BORN_CHARGES=BORN_file['born']
			self.__phonon.set_nac_params(BORN_file)
		#frequencies and eigenvectors at Gamma
		self.__frequencies,self.__eigvecs=self.__phonon.get_frequencies_with_eigenvectors([0, 0, 0])
	


		self.__NumberOfBands=len(self.__frequencies)

		#Nicer format of the eigenvector file
		self.__FormatEigenvectors()

		#Get dipole approximation of the intensitiess
		self.__set_Contributions()
		
		
	


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
		band (int): number of the frequency (ordered by energy)
		"""
		return self.__PercentageAtom[band,atom]
		
	def write_file(self,filename="Contributions.txt"):
		file  = open(filename, 'w')
		file.write('Frequency Atomare Beitraege \n')		
		for freq in range(len(self.__frequencies)):
			file.write('%s ' % (self.__frequencies[freq]))
			for atom in range(self.__natoms):
				file.write('%s ' % (self.__get_Contributions(freq,atom)))
			file.write('\n ')

		file.close()

	def plot(self,grouping,freqstart=[],freqend=[],freqlist=[],labelsforfreq=[],filename="Plot.eps"):
		"""
		freqstart min frequency of plot in cm-1
		freqend max frequency of plot in cm-1
		freqlist list of frequencies that will be plotted
		"""
		fig, ax1 = plt.subplots()		
		p={}
		summe={}
		if freqlist==[]:		
			freqlist=range(len(self.__frequencies))
			
		else: 	
			#Todo: Frequenzliste: von jedem Element 1 abziehen	
			for freq in range(len(freqlist)):
				freqlist[freq]=freqlist[freq]-1
		
		for group in range(len(grouping['GroupedAtoms'])):	
			color1=grouping['ColorsOfGroupedAtoms'][group]
			Entry={}		
			for freq in range(len(freqlist)):
				Entry[freq]= 0
			for number in grouping['GroupedAtoms'][group]:	
				#set the first atom to 0
				atom=int(number)-1
				for freq in range(len(freqlist)):
					
         	       			Entry[freq]= Entry[freq]+ self.__get_Contributions(freqlist[freq],atom)
					if group==0:
						summe[freq]=0	
			
			#plot bar chart
			p[group]=ax1.barh(np.arange(len(freqlist)),Entry.values(),left=summe.values(),color=color1,height=1,label=grouping['Legend'][group] ) 
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
		ax1.set_ylabel('Wavenumber (cm-1)')
		ax1.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",mode="expand", borderaxespad=0, ncol=len(grouping['GroupedAtoms']))
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
 
#more functions planned: only show IR or Raman active modes automatically
#listen uebergeben, die frequenznummern enthalten und nur diese plotten
