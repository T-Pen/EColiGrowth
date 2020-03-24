from opentrons import robot, labware, instruments, protocol_api
import numpy as np

# metadata
metadata = {
    'protocolName': 'fillhgt5000_test',
    'authors': 'Tobias Pfennig (topfe101@hhu.de)',
    'description': '\
        Protocol aimed at testing the fillght_falcon50000 function after determining the adjustment factor',
}


# ============================= fillhgt function =============================

def fillhgt_falcon50000(volume, diameter=27, cone_height=14.5, cone_tipdiameter=5):
    # measurements: 30 x 150 mm

    cone_slope = (diameter / 2 - cone_tipdiameter / 2) / cone_height
    cone_maxvol = 0.33 * np.pi * cone_height * ((diameter / 2) ** 2 + (diameter / 2) * (cone_tipdiameter / 2) + (cone_tipdiameter / 2) ** 2)
    cone_heighttip = cone_tipdiameter / (2 * cone_slope)
    cone_voltip = 0.33 * np.pi * cone_heighttip * (cone_tipdiameter / 2) ** 2

    ##catch errors##
    if volume > 50000:
        raise ValueError("volume exceeds 50 ml falcon volume")
    if volume < 0:
        raise ValueError("negative volumes not possible")

    ##distribute volume to modeled parts##
    if volume <= cone_maxvol:
        cone_vol = volume
    else:
        cone_vol = cone_maxvol

    cyl_vol = volume - cone_vol

    ##cone tip##
    # volume: 1/3 * pi * h * r^2, r = m * h
    if cone_vol == cone_maxvol:
        fillhgt_cone = cone_height
    else:
        cone_fullvol = volume + cone_voltip
        fillhgt_conefull = ((3 * cone_fullvol) / (np.pi * cone_slope ** 2)) ** 0.33
        fillhgt_cone = fillhgt_conefull - cone_heighttip

    ##top cylinder##
    if cyl_vol == 0:
        fillhgt_cyl = 0
    else:
        fillhgt_cyl = cyl_vol / (np.pi * (diameter / 2) ** 2)

    ##add heights##
    fillhgt = fillhgt_cone + fillhgt_cyl
    
    return fillhgt


# =============================== Load labware ===============================

tip_name = 'Biozym-tiprack-200ul'       # Add Biozym 200 ul tiprack to labware library
if tip_name not in labware.list():
    custom_tip = labware.create(
        tip_name,                       # name of you labware
        grid=(12, 8),                   # specify amount of (columns, rows)
        spacing=(9, 9),                 # distances (mm) between each (column, row)
        diameter=5.23,                  # diameter (mm) of each well on the plate
        depth=50,                       # depth (mm) of each well on the plate
        volume=200)

tiprack = labware.load(tip_name, '9')   # Load tiprack

trough = labware.load('opentrons_6_tuberack_falcon_50ml_conical', '1')  # Load reservoir

pipette = instruments.P300_Single(mount='right', tip_racks=[tiprack], max_volume=200)  # Load pipette


# =========================== Set source locations ============================

EmptyFalcon = trough.wells('B1')        #50 ml falcon left EMPTY for dispension of water

WaterFalcon = trough.wells('B2')        #50 ml falcon with water for testing transferring and immersion depth

pipette.start_at_tip(tiprack['A1'])


# ============================ Define Parameters =============================

liqHeightAdjust = 0.5                   # added to calculated fillheight to adjust for position of the falcon in mm (from fillhgt50000_calibrate)
tipDepth = 3                            # desired tipdepth in liquid AFTER aspiring in mm
lowestHeight = 5                        # lowest accepted tip height above the falcon bottom in mm

WaterFalcon_vol = 50000                 #initial water volume of the water filled falcon in ul


# ============================== Transfer water ===============================

pipette.pick_up_tip()

pipette.transfer(10, EmptyFalcon.bottom(liqHeightAdjust), EmptyFalcon.bottom(fillhgt_falcon50000(50000)+liqHeightAdjust), new_tip='never')
robot.pause()

volumes = np.arange(0,WaterFalcon_vol,200)[::-1]

for vol in volumes:
    
    currentDist = fillhgt_falcon50000(vol) + liqHeightAdjust - tipDepth
    if currentDist < lowestHeight:
        break
    
    pipette.transfer(200,
                     WaterFalcon.bottom(currentDist),
                     EmptyFalcon.top(),
                     new_tip = 'never')

robot.home()