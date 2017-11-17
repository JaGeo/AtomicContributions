from ContributionsOfAtomsToModes  import AtomicContributionToModes 


test=AtomicContributionToModes(PoscarName='POSCAR',ForceConstants=False,ForceFileName='FORCE_SETS',supercell=[[3, 0, 0],[0, 3, 0], [0, 0, 4]])

test.write_file();
grouping={}
#atom numbers start at 1, not at 0

#grouping of atoms for plots
grouping['GroupedAtoms']=[[1,2],[3,4],[5,6,7,8],[9,10,11,12,13,14,15,16]]
grouping['ColorsOfGroupedAtoms']=['black','red','green','grey']
grouping['Legend']=['C','O','N','H']
freqlist=[1,2,3,4,5,6,7,8,9]

#Symmetrien der einzelnen Moden abfragen und eintragen?
#Auswahl verschiedener Frequenzen? Einbau ueber Freqrange und eine neue Frequenzliste, die nur noch Frequenzen innerhalb des Bereiches enthaelt 
#Entartete Moden rauslassen
test.plot(grouping=grouping,freqstart=0.0,freqend=150,freqlist=freqlist)# ,freqlist=freqlist) #,freqend=4000.0) # ,freqstart=0.0,freqend=1000.0);
