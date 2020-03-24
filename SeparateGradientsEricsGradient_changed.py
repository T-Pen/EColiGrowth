
from opentrons import robot, labware, instruments, protocol_api
#from opentrons import 
#import numpy as np
#from copy import deepcopy






# metadata
metadata = {
    'protocolName': 'SeparateGradientsErisGradient',
    'author': 'Eric Behle eric.behle@hhu.de',
    'description': '\
        Protocol designed to pipette gradients of acetate and glucose (or any two substances) respectively\
        on an 48 well plate.\
        Rows 1-3: triplicate glucose gradients up to last column;\
        Rows 4-6: triplicate acetate gradients up to last column\
        Last column: both glucose and acetate are varied at fixed overall carbon concentration.\
        CAREFUL: the Biozym-tiprack used in this version of the script is too small to fit securely\
        into the robot. Use caution when placing it. Currently aligned so that it works best when moved into upper left corner',
}

# =============================== Define custom functions =================================

def ericsVolumeGradient(vol1, vol2, wellNumber, pipette, source, destination, sourceHeight, destHeight):
    """ Creates a volume gradient with a choosable height of pipette above source and destination;
        This was written due to a flaw in the api command for gradient pipetting
    kwargs:
         vol1 --> volume at beginning of gradient (double)       
         vol2 --> volume at end of gradient (double)
         wellNumber --> number of wells included in gradient, i.e. number of wells in destination (int)
         pipette --> pipette used (opentrons.legacy_api.instruments.pipette.Pipette)
         source --> Source from which to pipette (opentrons.legacy_api.containers.placeable.Well)
         destination --> list of wells in which to pipette (opentrons.legacy_api.containers.placeable.WellSeries)
         sourceHeight --> Distance of pipette tip from BOTTOM of source well (int)
         destHeight --> Distance of pipette tip from BOTTOM of destination well (int)
    """
    minVol = 1
    if wellNumber == 1:
        pipette.transfer(abs(vol1-vol2),
                         source.bottom(sourceHeight),
                         destination.bottom(destHeight),
                         new_tip='never'
                         )
        return
    
    dV = -1 * (vol1-vol2)/(wellNumber-1)
    for i in range(wellNumber):
        if vol1+i*dV >= minVol:
            pipette.transfer(vol1 + i * dV,
                     source.bottom(sourceHeight),
                     destination(i).bottom(destHeight),
                     new_tip='never'
                     )
    return 


#================================ Load labware ==================================

plate_name = 'biolector-flower-plate-48-well'
if plate_name not in labware.list():
    custom_plate = labware.create(
        plate_name,                    # name of labware
        grid=(8, 6),                    # amount of (columns, rows)
        spacing=(13, 13),               # distances (mm) between each (column, row)
        diameter=10.85,                     # diameter (mm) of each well on the plate
        depth=33,                       # depth (mm) of each well on the plate
        volume=3200)

tip_name = 'Biozym-tiprack-200ul' # Add Biozym 200 ul tiprack to labware library
if tip_name not in labware.list():
    custom_tip = labware.create(
        tip_name,                    # name of you labware
        grid=(12, 8),                    # specify amount of (columns, rows)
        spacing=(9, 9),               # distances (mm) between each (column, row)
        diameter=5.23,                     # diameter (mm) of each well on the plate
        depth=50,                       # depth (mm) of each well on the plate
        volume=200)

plate = labware.load(plate_name, '7') # Load well plate

tiprack = labware.load(tip_name, '11')# Load tiprack
#tiprack = labware.load('opentrons-tiprack-300ul', '11')


trough = labware.load('opentrons_6_tuberack_falcon_50ml_conical', '10') # Load reservoir

# pipettes
pipette = instruments.P300_Single(mount='left', tip_racks=[tiprack], max_volume = 200) # Load pipette



#========================== Define parameters ==========================================

tipLength = 51# Length of pipette tips in mm
liqDistTop = 20# Length of pipette tip which should remain above liquid
lowestHeight = 2.5# Lowest height on which to set the pipette tip in mm


Nrow = 6 # Number of rows
Ncolumn = 8 # Number of columns

maxPerWell = 1000# Total liquid volume inside each well
minGlucose = 0# Minimal amount of glucose 
maxGlucose = 0.9*maxPerWell# Maximal amount of glucose
minAcetate = 0# Minimal amount of acetate
maxAcetate = 0.9*maxPerWell# Maximal amount of acetate




#====================================== Notes ============================================
#Run-time: roughly 31 minutes with 300 ul tips. Not tested with 200 ul tips yet



#============================= Set source locations ==================================

glucoseSource = trough.wells('A1')# M9 medium with max concentration of glucose
acetateSource = trough.wells('A2')# M9 medium with max concentration of acetate
M9Source = trough.wells('A3')# M9 medium for diulution
glcCultureSource = trough.wells('B1') # pre-culture grown in glucose
aceCultureSource = trough.wells('B2')# pre-culture grown in acetate


#========================== ADAPT THESE PARAMETERS FOR EACH RUN ========================
#========================== Define tip at which to start pipetting =====================

pipette.start_at_tip(tiprack['A1'])

#============================== Define height of liquid in source reservoirs (currently 50 ml falcons)

liqHeightAceCulture = 30 # mm
liqHeightGlcCulture = 20 # mm # NOT USED!!
liqHeightM9 = 100 # mm
liqHeightGlc = 55 # mm
liqHeightAce = 55 # mm




#============================= Define tip distances from bottom of plates ===================

falconHeight = 115 #Height of a falcon in mm
initHeight = falconHeight # initial height of pipette above reservoirs

tipDistGlcCulture = initHeight # Initial distance of pipette tip to BOTTOM of falcon tube for glucose pre-culture
tipDistAceCulture = initHeight # Initial distance of pipette tip to BOTTOM of falcon tube for acetate pre-culture [mm]
tipDistM9 = initHeight # Initial distance of pipette tip to BOTTOM of falcon tube for M9 [mm]
tipDistGlucose = initHeight # Initial distance of pipette tip to BOTTOM of falcon tube for M9 with glucose [mm]
tipDistAcetate = initHeight # Initial distance of pipette tip to BOTTOM of falcon tube for M9 with acetate [mm]


if liqHeightAceCulture - tipLength > 0:
    tipDistAceCulture = liqHeightAceCulture - tipLength + liqDistTop
else:
    tipDistAceCulture = lowestHeight

if liqHeightGlcCulture - tipLength > 0:
    tipDistGlcCulture =liqHeightGlcCulture - tipLength + liqDistTop
else:
    tipDistGlcCulture = lowestHeight

if liqHeightM9 - tipLength > 0:
    tipDistM9 = liqHeightM9 - tipLength + liqDistTop
else:
    tipDistM9 = lowestHeight

if liqHeightGlc - tipLength > 0:
    tipDistGlucose = liqHeightGlc - tipLength + liqDistTop
else:
    tipDistGlucose = lowestHeight
    
if liqHeightAce - tipLength > 0:
    tipDistAcetate = liqHeightAce - tipLength + liqDistTop
else:
    tipDistAcetate = lowestHeight

currentDist = 0 # Current distance of pipette tip from BOTTOM of falcon. Initialized here and set inside each following for-loop

# In a falcon, 5000 ul of volume change correspond to a height change of roughly 9 mm
# --> Reduce the pipette height accordingly. It will be set to 10 mm per 5000 ul for convenience

volChangeFalcon = 5000 # 5000ul
distChangeFalcon = 10 # Distance change per volChangeFalcon (rough value)
rowHeightReduction =  0# Reduction of pipette tip height after each column. Updated depending on total volume pipetted
rowHeightReduction = 0 # Reduction of pipette tip height after each column. Updated depending on total volume pipetted

tipDistFlower = 20 # Distance of pipette tip to BOTTOM of flower plate; used to avoid cross-contamination

#========================== Distribute culture ==========================



#Calculate reduction per column based on amount of liquid pipetted
rowHeightReduction =  distChangeFalcon * (sum(range(Nrow-2))*(maxPerWell - (maxGlucose-minGlucose))/(Nrow-2)) / volChangeFalcon

# Distribute glcCulture
source = aceCultureSource
pipette.pick_up_tip()# Pick up a tip


for i in range(int(Nrow/2)):
    currentDist = tipDistAceCulture - i*rowHeightReduction    
    if currentDist < lowestHeight:
        currentDist = lowestHeight;

    pipette.distribute(maxPerWell - maxGlucose,
                       source.bottom(currentDist),
                       plate.rows(i)[:-1],
                       new_tip='never')
pipette.drop_tip()# Drop current tip

# Distribute aceCulture
source = aceCultureSource
pipette.pick_up_tip()# Pick up a tip

for i in range(int(Nrow/2), Nrow): 
    currentDist = tipDistAceCulture - (i-Nrow/2)*rowHeightReduction
    if currentDist < lowestHeight:
        currentDist = lowestHeight;

    pipette.distribute(maxPerWell - maxGlucose,
                       source.bottom(currentDist),
                       plate.rows(i)[:-1],
                       new_tip='never')
# acetate-grown pre-culture for last column


""" # NO CULTURE IN LAST COLUMN 
currentDist = tipDistAceCulture - (3) * rowHeightReduction
if currentDist < lowestHeight:
    currentDist = lowestHeight
pipette.distribute(maxPerWell - maxGlucose,
                   source.bottom(currentDist),
                   plate.columns(Ncolumn-1),
                   new_tip='never')
"""
pipette.drop_tip()# Drop current tip


# ================== Gradients for first 7 columns =====================

source = M9Source # Distribute M9 medium for dilution

rowHeightReduction = distChangeFalcon * (sum(range(Ncolumn-2))*(maxGlucose-minGlucose)/(Ncolumn-2))/volChangeFalcon

pipette.pick_up_tip()
for i in range(Nrow):
    currentDist = tipDistM9 - i * rowHeightReduction
    if currentDist < lowestHeight:
        currentDist = lowestHeight;

    ericsVolumeGradient(minGlucose,
                        maxGlucose,
                        len(plate.rows(i)[:-1]),
                        pipette,
                        source,
                        plate.rows(i)[:-1],
                        currentDist,
                        tipDistFlower
                        )
pipette.drop_tip()


source = glucoseSource # Dilute glucose

pipette.pick_up_tip()
for i in range(int(Nrow/2)):
    currentDist = tipDistGlucose - i * rowHeightReduction
    if currentDist < lowestHeight:
        currentDist = lowestHeight

    ericsVolumeGradient(maxGlucose,
                        minGlucose,
                        len(plate.rows(i)[:-1]),
                        pipette,
                        source,
                        plate.rows(i)[:-1],
                        currentDist,
                        tipDistFlower
                        )
pipette.drop_tip()

    
source = acetateSource # Dilute acetate
pipette.pick_up_tip()
for i in range(int(Nrow/2), Nrow):
    currentDist = tipDistAcetate - (i-int(Nrow/2)) * rowHeightReduction
    if currentDist < lowestHeight:
        currentDist = lowestHeight

    ericsVolumeGradient(maxAcetate,
                        minAcetate,
                        len(plate.rows(i)[:-1]),
                        pipette,
                        source,
                        plate.rows(i)[:-1],
                        currentDist,
                        tipDistFlower
                        )
pipette.drop_tip()        
#====================== Last column ========================================
        
currentColumn = plate.columns(Ncolumn-1)

source = M9Source # Distribute M9



currentDist = tipDistM9 - 6 * rowHeightReduction
if currentDist < lowestHeight:
    currentDist = lowestHeight
pipette.distribute(
        maxPerWell - maxGlucose/2 ,
        source,
        currentColumn
        )



source = glucoseSource # Dilute glucose
currentDist = tipDistGlucose - 3*rowHeightReduction
if currentDist < lowestHeight:
    currentDist = lowestHeight

pipette.pick_up_tip()
ericsVolumeGradient(maxGlucose/2,
                    minGlucose,
                    len(currentColumn),
                    pipette,
                    source,
                    currentColumn,
                    currentDist,
                    tipDistFlower
                    )
pipette.drop_tip()


source = acetateSource # Dilute acetate
currentDist = tipDistAcetate - 3*rowHeightReduction
if currentDist < lowestHeight:
    currentDist = lowestHeight

pipette.pick_up_tip()
ericsVolumeGradient(minAcetate,
                    maxAcetate/2,
                    len(currentColumn),
                    pipette,
                    source,
                    currentColumn,
                    currentDist,
                    tipDistFlower
                    )
pipette.drop_tip()

"""
Division by zero error when using built-in gradient method:
What works: 
    transfer or distribute with specified volume and either destination.bottom() or no bottom()
    transfer with volume gradient and no bottom() for destination
What does not work:
    transfer with volume gradient tuple and destination.bottom()
Reason: 
    create_volume_gradient uses len(destination); len(destination) = Nwells,
    but len(destination.bottom(x)) = 2. This is because it is a list [destination, coordinates]
"""
