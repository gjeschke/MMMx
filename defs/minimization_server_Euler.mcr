# set number of parallel threads, adapt to the number of threads that you want to use
SetThreads 60
# set energy unit
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
RemoveObj Water
# Get final energies and quality Z-scores
endenergy = Energy
solvenergy = SolvEnergy
endenergy=endenergy+solvenergy
endzscore=0.
for check in checklist
  zscore = Check (check)
  endzscore=endzscore+zscore
endzscore=endzscore/count checklist
# Add water again
AddObj Water
# Determine radius of water shell
r=RadiusObj Water
len=15
DelObj SimCell
# Add starting structure again
AddObj 1
# Place the two objects sufficiently far apart
MoveObj !1,X=(r+len)
MoveObj 1,X=(-r-len)
# Link them with an arrow
ShowArrow Start=Point,X=(-len),Y=0,Z=50,End=Point,X=(len),Y=0,Z=50,Radius=1,Color=Yellow
# Fit scene to screen
ZoomAll Steps=0
Move Z=20
# Label start and end structures with energies and Z-scores
namelist='Start','End'
for i=1 to 2
  name=namelist(i)
  NameObj (i),(name)
  LabelObj (i),(name),3.5,Color=Yellow,Y=5,Z=-10
  LabelObj (i),'Energy: (0.0+(name)energy) (EnergyUnit)',2.5,Color=Yellow,Y=0,Z=-10
  LabelObj (i),'Score: (0.00+(name)zscore)',2.5,Color=Yellow,Y=-5,Z=-10
Style Ribbon,BallStick
StickRes Water
# Make sure that the water shell is in the same coordinate system as the solute
TransferObj Water,End,Local=Fix

