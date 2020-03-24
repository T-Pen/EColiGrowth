#%%
from opentrons import robot, labware, instruments, protocol_api, util
import numpy as np

def transfer_fillhgt(volume, source, dest, pipette, fillhgtfunc, fillvol, pipParams, **kwargs):
    """
    Changed version of default "transfer" function.
    Transfers static volumes from one source to one or more destinations.
    During transfer the fillheight in the source is calculated after every pipetting.
    This fillheight is used to ensure minimal tip immersion of around tipDepth.
    """
    
    global _new_fillvol
    
    liqHeightAdjust, tipDepth, lowestHeight = pipParams #unpack parameters
    
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
        
        _new_fillvol = fillvol #save the updated fillvol as _new_fillvol to use further
        
    #cited from "pipette.transfer"
    pipette._run_transfer_plan(tips, plan, **kwargs)
    
    return(pipette)
 #%%
test = transfer_fillhgt(50000,
                 DyeSource,
                 plate.wells("A1"),
                 pipette,
                 fillhgt_falcon50000,
                 50000,
                 [2,5,5,])