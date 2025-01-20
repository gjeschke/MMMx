EnergyUnit kJ/mol
# Choose a nice orientation, get rid of water
NiceOri
DelRes Water
# Create a duplicate of the initial object for comparison
DuplicateObj 1
RemoveObj 1
# Clean for simulation, create cell
Clean
Cell Auto,Extension=20
# Enable automatic correction of cis-peptide bonds and wrong isomers
CorrectCis On
CorrectIso On
# Use YASARA2 force field
ForceField YASARA2,SetPar=Yes
# Get initial energies and quality Z-scores
startenergy = Energy
solvenergy = SolvEnergy
startenergy=startenergy+solvenergy
startzscore=0.
checklist='Dihedrals','Packing1D','Packing3D'
for check in checklist
  zscore = Check (check)
  startzscore=startzscore+zscore
startzscore=startzscore/count checklist
# If there is a protein, optimize side-chain rotamers
protresidues = CountRes Protein
if protresidues
  OptimizeRes Protein,Method=SCWALL
# Optimize the hydrogen-bonding network
OptHyd Method=YASARA
# Create a water shell (force field is not suitable for use in vacuo)
Experiment Neutralization
  NaCl 0
Experiment On
Wait ExpEnd
DelRes Water with distance>6 from not Water
# Perform energy minimization
Experiment Minimization
  Convergence 0.05
Experiment On
Wait ExpEnd
# Delete counter ions and water that moved away
DelRes CIP CIM with distance>6 from not HOH CIP CIM
DelRes Water with distance>6 from !Water
ShowAll
TransferObj 2 4,2,Local=Fix

