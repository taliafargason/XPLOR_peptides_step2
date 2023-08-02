#(opts,args) = xplor.parseArguments(["quick"]) # check for command-line typos

#quick=False
#for opt in opts:
 #   if opt[0]=="quick":  #specify -quick to just test that the script runs
  #      quick=True
   #     pass
   # pass


numberOfStructures=100
numTags=3 
startingStructure=0
ensembleSize=10
RRM1 = "resid 16:90"
RRM2 = "resid 121:196"

#if quick:
 #   numberOfStructures=3
  #  pass

# protocol module has many high-level helper functions.
import protocol

protocol.initRandomSeed(3421)   #explicitly set random seed

startPDB='best_single_struct+RS8b.pdb'
protocol.loadPDB(startPDB)

protocol.addUnknownAtoms()

#duplicate tag atoms
xplorSel=" or name ".join("N O C CA HA HN".split())
for i in range(1,numTags+1):
    xplor.command(
  f'''duplicate 
        sele=(segid "A" and resname CYSP and not (name {xplorSel}) ) 
        segid=ZP{i} 
    end''')
    pass
xplor.simulation.deleteAtoms(
    f'segid "A" and not (name {xplorSel}) and resname CYSP'
                            )

from avePot import AvePot
from xplorPot import XplorPot
from ensembleSimulation import EnsembleSimulation
esim = EnsembleSimulation("ensemble",ensembleSize)
esim.setAveType('sum')

#
# annealing settings
#

command = xplor.command

protocol.fixupCovalentGeom(maxIters=1000,useVDW=1)


# a PotList contains a list of potential terms. This is used to specify which
# terms are active during refinement.
#
from potList import PotList
potList = PotList()
# parameters to ramp up during the simulated annealing protocol
crossTerms = PotList()
from simulationTools import *

rampedParams=[]
highTempParams=[]

#measure of difference between ensemble members - reported, but not used in
# energy calculation
from posRMSDPotTools import RAPPot
posRMSD = RAPPot("posRMSD",f"name C CA N O")
crossTerms.append(posRMSD)


#rap term to keep members of ensemble in similar orientations
#from posRMSDPotTools import RAPPot
#rap = RAPPot("rap","name C CA N")
#rap.setScale( 100.0 )
#rap.setPotType( "square" )
#rap.setTol( 0.5 )
#potList.append(rap)
#crossTerms.append(rap)

#compare atomic Cartesian rmsd with a reference structure
#  backbone and heavy atom RMSDs will be printed in the output
#  structure files
#
#from posDiffPotTools import create_PosDiffPot
#refRMSD = create_PosDiffPot("refRMSD","name CA C N",
#                            pdbFile=startPDB,
#                            cmpSel=None)

# orientation Tensor - used with the dipolar coupling term
#  one for each medium
#   For each medium, specify a name, and initial values of Da, Rh.
#
from varTensorTools import *
media={}
#                        medium  Da   rhombicity
for (medium,Da,Rh) in [ ('t',   -20.36, 0.62),
                        ('b',   -20.36, 0.62)]:
    oTensor = create_VarTensor(medium)
    oTensor.setDa(Da)
    oTensor.setRh(Rh)
    media[medium] = oTensor
    pass
    

# dipolar coupling restraints for protein amide NH.  
#
# collect all RDCs in the rdcs PotList
#
# RDC scaling. Three possible contributions.
#   1) gamma_A * gamma_B / r_AB^3 prefactor. So that the same Da can be used
#      for different expts. in the same medium. Sometimes the data is
#      prescaled so that this is not needed. scale_toNH() is used for this.
#      Note that if the expt. data has been prescaled, the values for rdc rmsd
#      reported in the output will relative to the scaled values- not the expt.
#      values.
#   2) expt. error scaling. Used here. A scale factor equal to 1/err^2
#      (relative to that for NH) is used.
#   3) sometimes the reciprocal of the Da^2 is used if there is a large
#      spread in Da values. Not used here.
#
from rdcPotTools import *
rdcs = PotList('rdc') 
for (medium,expt,file,                 scale) in \
    [('t','NH' ,'RDC_unphos_sup_conf.tbl'       ,1)]:
    rdc = create_RDCPot("%s_%s"%(medium,expt),file,media[medium])

    #1) scale prefactor relative to NH
    #   see python/rdcPotTools.py for exact calculation
    # scale_toNH(rdc) - not needed for these datasets -
    #                        but non-NH reported rmsd values will be wrong.

    #3) Da rescaling factor (separate multiplicative factor)
    # scale *= ( 1. / rdc.oTensor.Da(0) )**2
    rdc.setScale(scale)
    rdc.setShowAllRestraints(1) #all restraints are printed during analysis
    rdc.setThreshold(1.5)       # in Hz
    rdcs.append(rdc)
    pass
potList.append(rdcs)
rampedParams.append( MultRamp(0.05,5.0, "rdcs.setScale( VALUE )") )

# calc. initial tensor orientation
# and setup tensor calculation during simulated annealing
#
for medium in list(media.keys()):
    calcTensorOrientation(media[medium])
    rampedParams.append( StaticRamp("calcTensor(media['%s'])" % medium) )
    pass

from prePotTools import create_PREPot

# Setup of backbone 1H PRE restraints  
#pre = create_PREPot("PRE","RS8PRE.tbl","normal")

pre = PotList('pre')
a = create_PREPot("RS8","RS8PRE.tbl","normal",fixTau=True)
a.setFunType("correlation")
#a.setTaucAtoms()
pre.append(a)
#b = create_PREPot("SR8_RRM2","RS8_RRM2.tbl","normal",fixTau=True)
#b.setFunType("correlation")
#b.setTaucAtoms()
#pre.append(b)

def runSBmode():
    print ('configuring SB mode')
    from prePotTools import setupSBmode
    from simulationTools import flattenPotList
    for p in flattenPotList(pre):
        setupSBmode(p)
        p.setTcType("fix")
        p.setEquType("sb")
        p.setTauC(16.8)
        p.setRlxType("r2dd")
        p.setSqn( 0.5 )
        p.setGfac( 2  )
        p.setFreqI( 850 )
        #p.setFunType("square")
        p.setSclType("sigma")
        p.calcEnergy() 
        pass
    return

def runSBMFmode():
    print ('configuring SBMF mode')
    from prePotTools import setupSBMFmode
    from simulationTools import flattenPotList
    for p in flattenPotList(pre):
        setupSBMFmode(p)
        p.calcEnergy()
        p.setTcType("fix")
        p.setTiType("fix")
        p.setTtType("fix")
        p.setEquType("sbmf")
        p.setSbmfType("taut")
        p.setTauC(16.8)
        p.setTauT(0.000)
        p.setTauI(0.000)
        p.setRlxType("r2dd")
        p.setSqn( 0.5 )
        p.setGfac( 2  )
        p.setFreqI( 850 )
        #p.setFunType("square")
        p.setSclType("sigma")
        p.calcEnergy()
        pass
    return

potList.append(pre)
rampedParams.append( MultRamp(100,1000, "pre.setScale( VALUE )") )
from prePotTools import fitTauc
rampedParams.append( StaticRamp("[fitTauc(term) for term in pre]" ) ) #to optimize TauC
#rampedParams.append( StaticRamp("runSBMFmode()") )
highTempParams.append( StaticRamp("runSBmode()") )




# set up NOE potential (CSPs)
noe=PotList('noe')
potList.append(noe)
from noePotTools import create_NOEPot
for (name,scale,file) in [('all',1,"CSPdata.tbl"),
                          #add entries for additional tables
                          ]:
    pot = AvePot(create_NOEPot(name,file,esim=esim.member()))
    pot.setPotType("soft") # if you think there may be bad NOEs
    pot.setScale(scale)
    noe.append(pot)
rampedParams.append( MultRamp(2,30, "noe.setScale( VALUE )") )

# Set up dihedral angles
protocol.initDihedrals("talosn.tbl",
                       #useDefaults=False  # by default, symmetric sidechain
                                           # restraints are included
                       )
potList.append( AvePot(XplorPot,'CDIH') ) 
highTempParams.append( StaticRamp("potList['CDIH'].setScale(10)") )
rampedParams.append( StaticRamp("potList['CDIH'].setScale(200)") )


# gyration volume term 
#
#from gyrPotTools import create_GyrPot
#gyr = create_GyrPot("Vgyr",
#                    #"resid 1:56" # selection should exclude disordered tails
#                   )
#potList.append(gyr)
#rampedParams.append( MultRamp(.002,1,"gyr.setScale(VALUE)") )


# HBPot - knowledge-based hydrogen bond term
#
from hbPotTools import create_HBPot
hb = AvePot(create_HBPot('hb',selection=AtomSel("not pseudo",esim.member())))
hb.setScale(2.5)
potList.append( hb )

#New torsion angle database potential
#
from torsionDBPotTools import create_TorsionDBPot
torsionDB = AvePot(create_TorsionDBPot('torsionDB',
                                       selection=AtomSel("all",esim.member())
                                       ))
potList.append( torsionDB )
rampedParams.append( MultRamp(.002,2,"torsionDB.setScale(VALUE)") )

#
# setup parameters for atom-atom repulsive term. (van der Waals-like term)
#
from repelPotTools import create_RepelPot,initRepel
repel = AvePot(create_RepelPot('repel',
                        selPairs=[("segid A","segid A"),
                                  ("segid A","segid ZP*")] +
                        [("segid A",f"segid {segid}") for
                         segid in [atom.segmentName() for atom in
                                   AtomSel("tag and segid ZP*")]],
                               selection=AtomSel("not PSEUDO",
                                                 esim.member())
                               ))
potList.append(repel)
rampedParams.append( StaticRamp("initRepel(repel,use14=False)") )
rampedParams.append( MultRamp(.004,4,  "repel.setScale( VALUE)") )
# nonbonded interaction only between CA atoms
highTempParams.append( StaticRamp("""initRepel(repel,
                                               use14=True,
                                               scale=0.004,
                                               repel=1.2,
                                               moveTol=45,
                                               interactingAtoms='name CA'
                                               )""") )

# Selected 1-4 interactions.
from  torsionDBPotTools import create_Terminal14Pot
repel14 = AvePot(create_Terminal14Pot('repel14',
                                      selection=AtomSel("not PSEUDO",
                                                        esim.member())))
potList.append(repel14)
highTempParams.append(StaticRamp("repel14.setScale(0)"))
rampedParams.append(MultRamp(0.004, 4, "repel14.setScale(VALUE)"))


potList.append( AvePot(XplorPot,"BOND") )
potList.append( AvePot(XplorPot,"ANGL") )
potList['ANGL'].setThreshold( 5 )
rampedParams.append( MultRamp(0.4,1,"potList['ANGL'].setScale(VALUE)") )
potList.append( AvePot(XplorPot,"IMPR") )


# set custom values of threshold values for violation calculation
#
potList['CDIH'].setThreshold( 5 ) #5 degrees is the default value, though

potList['IMPR'].setThreshold( 5 )
rampedParams.append( MultRamp(0.1,1,"potList['IMPR'].setScale(VALUE)") )
      


# Give atoms uniform weights, except for the anisotropy axis
#
protocol.massSetup()


for m in list(media.values()):
#    m.setFreedom("fixDa, fixRh")        #fix tensor Rh, Da, vary orientation
    m.setFreedom("varyDa, varyRh")      #vary tensor Rh, Da, vary orientation

# IVM setup
#   the IVM is used for performing dynamics and minimization in torsion-angle
#   space, and in Cartesian space.
#
from ivm import IVM
dyn  = IVM() 

# initialize ivm topology for torsion-angle dynamics

dyn.group(f"({RRM2}) and name C CA C N")
dyn.fix(f"({RRM1}) and name C CA C N")

protocol.torsionTopology(dyn)

# minc used for final cartesian minimization
#
minc = IVM()
protocol.initMinimize(minc)

minc.group(f"({RRM2}) and name C CA C N")
minc.fix(f"({RRM1}) and name C CA C N")

for m in list(media.values()):
    m.setFreedom("varyDa, varyRh")    #allow all tensor parameters float here
    pass
protocol.cartesianTopology(minc)



# object which performs simulated annealing
#
from simulationTools import AnnealIVM
init_t  = 3000.     # Need high temp and slow annealing to converge
cool = AnnealIVM(initTemp =init_t,
                 finalTemp=25,
                 tempStep =12.5,
                 ivm=dyn,
                 rampedParams = rampedParams)

def accept(potList):
    """
    return True if current structure meets acceptance criteria
    """
    #if potList['noe'].violations()>0:
      #  return False
    #if potList['rdc'].rms()>1.2: #this might be tightened some
     #   return False
    #if potList['CDIH'].violations()>0:
     #   return False
    if potList['BOND'].violations()>0:
        return False
    if potList['ANGL'].violations()>0:
        return False
    if potList['IMPR'].violations()>1:
        return False
    
    return True

def calcOneStructure(loopInfo):
    """ this function calculates a single structure, performs analysis on the
    structure, and then writes out a pdb file, with remarks.
    """

    #write one PSF each time this script is run
    if loopInfo.structNum==0:
        if esim.singleThread():
            from simulationTools import genFilename
            xplor.fastCommand(f"write PSF output={genFilename('SCRIPT.psf')} end")
            pass
        esim.multiThread()
        pass

    from monteCarlo import randomizeTorsions
    randomizeTorsions(dyn)

    # initialize parameters for high temp dynamics.
    InitialParams( rampedParams )
    # high-temp dynamics setup - only need to specify parameters which
    #   differfrom initial values in rampedParams
    InitialParams( highTempParams )

    # high temp dynamics
    #
    protocol.initDynamics(dyn,
                          potList=potList, # potential terms to use
                          bathTemp=init_t,
                          initVelocities=1,
                          finalTime=10,    # stops at 10ps or 5000 steps
                          numSteps=5000,   # whichever comes first
                          printInterval=100)

    dyn.setETolerance( init_t/100 )  #used to det. stepsize. default: t/1000 
    dyn.run()

    # initialize parameters for cooling loop
    InitialParams( rampedParams )


    # initialize integrator for simulated annealing
    #
    protocol.initDynamics(dyn,
                          potList=potList,
                          numSteps=200,       #at each temp: 200 steps or
                          finalTime=.4 ,       # .4ps, whichever is less
                          printInterval=100)

    # perform simulated annealing
    #
    cool.run()
              
              
    # final torsion angle minimization
    #
    protocol.initMinimize(dyn,
                          printInterval=50)
    dyn.run()

    # final all- atom minimization
    #
    protocol.initMinimize(minc,
                          potList=potList,
                          dEPred=10)
    minc.run()

    #do analysis and write structure when this function returns
    pass
    #enePre = pre.calcEnergy()
    #qPre = pre.Qfactor()
    #rPre = pre.rms()
    #tcPre = pre.tc() * 1.0e9
    #print (pre.info())
    #print (pre.showRestraints(0))

    #loopInfo.writeStructure(potList)

    #pass



from simulationTools import StructureLoop, FinalParams
StructureLoop(numStructures=numberOfStructures,
              structLoopAction=calcOneStructure,
              pdbTemplate="ENSEMBLE_%d_STRUCTURE_MEMBER.sa" %ensembleSize,
              calcMissingStructs=True, #calculate only missing structures
              doWriteStructures=True,  #analyze and write coords after calc
              genViolationStats=True,
              averagePotList=potList,
              averageSortPots=[potList['BOND'],potList['ANGL'],potList['IMPR'],
                               noe,rdcs,potList['CDIH']],
              averageCrossTerms=crossTerms,
              averageTopFraction=0.5, #report only on best 50% of structs
              averageAccept=accept,   #only use structures which pass accept()
              averageContext=FinalParams(rampedParams),
              averageFilename="SCRIPT_ave.pdb",    #generate regularized ave structure
              averageFitSel=f"({RRM1}) and name CA",
              averageCompSel="not resname ANI and not name H*"     ).run()