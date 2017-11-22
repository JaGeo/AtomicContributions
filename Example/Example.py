from ContributionsOfAtomsToModes  import AtomicContributionToModes 


test=AtomicContributionToModes(PoscarName='POSCAR',ForceConstants=False,ForceFileName='FORCE_SETS',supercell=[[3, 0, 0],[0, 3, 0], [0, 0, 4]])


#write a file with atomic contributions to each mode
test.write_file();


#now plots: 
#grouping of atoms for plots
grouping={}
#atom numbers start at 1, not at 0





atomgroups=[[1,2],[3,4],[5,6,7,8],[9,10,11,12,13,14,15,16]]
colorofgroups=['black','red','green','grey']
legendforgroups=['C','O','N','H']

test.plot(atomgroups=atomgroups,colorofgroups=colorofgroups,legendforgroups=legendforgroups)
test.plot_irred(atomgroups=atomgroups,colorofgroups=colorofgroups,legendforgroups=legendforgroups,irreps=['B2','E'])

#lists all numbers of frequencies that are plotted, frequency numbers start at 1, not at 0
#freqlist=[1,2,3,4]

#Uebergebe Mulliken Symbole IR-aktiver Moden

#IRactiveModes= [1,3,6,9,11,14,17,19,23,25,27,29,31,33,36,37,40,41,44,45,48]			

#LabelsForIRactiveModes=['E','B2','E','E','E','E','E','B2','E','E','B2','E','B2','E','B2','E','B2','E','B2','E','B2']



#Symmetrien der einzelnen Moden abfragen und eintragen?
#Auswahl verschiedener Frequenzen? Einbau ueber Freqrange und eine neue Frequenzliste, die nur noch Frequenzen innerhalb des Bereiches enthaelt 
#Entartete Moden rauslassen
#test.plot(grouping=grouping,labelsforfreq=LabelsForIRactiveModes,freqlist=IRactiveModes) #,freqstart=0.0,freqend=150,freqlist=freqlist)
