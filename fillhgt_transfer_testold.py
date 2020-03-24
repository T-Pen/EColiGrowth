from opentrons import robot, labware, instruments, protocol_api, util
import numpy as np

# metadata
metadata = {
    'protocolName': 'OT_testpipetting',
    'authors': 'Tobias Pfennig (topfe101@hhu.de)',
    'description': '\
        Protocol aimed at testing he accuracy and precision of the OT2 before using it for pipetting.\
        A 96-well plate will be pipetted with dye on top of hand-pipetted water.\
        The left half will be pipetted by the OT2, the right half by hand as control.\
        The first two columns of each half receive a static dye baseline, he following three columns\
        are pipetted as a gradient. The last column acts as a blank with only water.',
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


def transfer_fillhgt(volume, source, dest, pipette, fillhgtfunc, fillvol, pipParams, **kwargs):
    """
    A funcion for transferring static volumes from one source to one or more destinations.
    During transfer the fillheight in the source is calculated after every pipetting.
    This fillheight is used to ensure minimal tip immersion of around tipDepth.
    """
    
    liqHeightAdjust, tipDepth, lowestHeight = pipParams
    
    source = source.top() # force inclusion of coordinates in transfer plan
    
    
    #cited from "pipette.transfer", substituting 'self' with 'pipette'
    kwargs['mode'] = kwargs.get('mode', 'transfer')

    touch_tip = kwargs.get('touch_tip', False)
    if touch_tip is True:
        touch_tip = -1
    kwargs['touch_tip'] = touch_tip

    tip_options = {
        'once': 1,
        'never': 0,
        'always': float('inf')
    }
    tip_option = kwargs.get('new_tip', 'once')
    tips = tip_options.get(tip_option)
    if tips is None:
        raise ValueError('Unknown "new_tip" option: {}'.format(tip_option))

    plan = pipette._create_transfer_plan(volume, source, dest, **kwargs)
    #end citation
    
    
    #adjust pipetting plan
    for i in range(len(plan)):
        curr_vol = plan[i]["aspirate"]["volume"]
        
        fillvol -= curr_vol  # update fillvolume to after the following pipetting
        
        currdist = fillhgtfunc(fillvol) + liqHeightAdjust - tipDepth
        if currdist < lowestHeight: 
            currdist = lowestHeight 
        
        newLoc = (plan[i]["aspirate"]["location"][0],
                  util.vector.Vector(
                          plan[i]["aspirate"]["location"][1][0],
                          plan[i]["aspirate"]["location"][1][1],
                          currdist))
        
        plan[i]["aspirate"]["location"] = newLoc
    
    #cited from "pipette.transfer"
    pipette._run_transfer_plan(tips, plan, **kwargs)
    
    return(pipette)


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


#================================== Pipette '==================================

pipette.pick_up_tip()

transfer_fillhgt(10000,
                 WaterFalcon,
                 EmptyFalcon.top(),
                 pipette,
                 fillhgt_falcon50000,
                 50000,
                 [2,5,5,])