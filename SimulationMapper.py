#%%
import collections, re
import numpy as np
from prettytable import PrettyTable

#================================= Parameters =================================

PC = False

if PC:
    header = "D:/"
else:
    header = "C:/Users/"

filepath = "tpeng/OneDrive/OD__Universität/BA__Axmann/protocol_simulations/"
filename = "P03_OTV1_GlcGradientONECELL_BigFalconRack_out"

#roundDigits = 2

ContainerDims = {3:(6,8),
                 5:(4,4),
                 6:(2,3)}


#==============================================================================

out_file = open(header+filepath+filename+".txt", "r")

deckList = collections.defaultdict(lambda:collections.defaultdict(lambda:collections.defaultdict(int)))
aspirate = dict()
dispense = dict()

def printSlot(slot=None, outputObject=False):
    
    # get simulated slots
    simSlots = sorted(list(deckTables.keys()))
    simSlots = [str(x) for x in simSlots]
    simSlots = ", ".join(simSlots)
    
    if slot is None:
        print("simulated slots: " + simSlots)
        return None
    
    slotTable = deckTables.get(slot)
    
    if slotTable is None:
        raise ValueError("no simulation in that slot\n\t\t\tsimulated slots: " + simSlots)
    elif outputObject:
        return(slotTable)
    else:
        print(slotTable)
        return None
#%%
for line in out_file:
    
    start = re.search("Aspirating|Dispensing", line)
    
    if start is None:
        continue
    else:
        line = line[start.start():]
    
    line_splt = re.split(" |\n", line)
    
    if line_splt[0] == "Aspirating":
        aspirate[" vol"] = float(line_splt[1])
        aspirate["well"] = line_splt[5]
        aspirate["slot"] = int(line_splt[7][1:-1])
        aspirate["label"] = str(aspirate["slot"]) + "-" + aspirate["well"]
        
        deckList[aspirate["slot"]][aspirate["well"]][" vol"] -= aspirate[" vol"]
        
        continue
    
    elif line_splt[0] == "Dispensing":
        dispense[" vol"] = float(line_splt[1])
        dispense["well"] = line_splt[5]
        dispense["slot"] = int(line_splt[7][1:-1])
        
        deckList[dispense["slot"]][dispense["well"]][" vol"] += dispense[" vol"]
        deckList[dispense["slot"]][dispense["well"]][aspirate["label"]] += dispense[" vol"]
        
        continue
#%%
slotKeys = list(deckList.keys())
deckVis = collections.defaultdict()
deckTables = collections.defaultdict()

for slotKey in slotKeys:
    slot = deckList[slotKey]
    
    wells = list(slot.keys())
    rows = [x[0] for x in wells]
    cols = [x[1:] for x in wells]
    
    colCheck = re.search(r'\D',"".join(cols))
    
    if colCheck is not None:
        wrongCol = colCheck.start()-1
        wrongWell = wells[wrongCol]
        
        raise ValueError("an incorrect well description was called:'" + wrongWell + "'")
    
    slotDim = ContainerDims.get(slotKey)
    
    rowMax = ord(max(rows)) - 64
    colMax = max([int(x) for x in cols])
    
    if slotDim is not None:
        if type(slotDim) in (tuple, list) :
            rowMax_lst, colMax_lst = slotDim
        else:
            raise ValueError("entries in CotainerDims have to be tuples or lists")
        
        if rowMax_lst < rowMax  or colMax_lst < colMax:
            raise ValueError("dimensions for slot {} in CotainerDims ({},{}) are smaller than found in simulation ({},{})".format(slotKey, rowMax_lst, colMax_lst, rowMax , colMax))
        else:
            rowMax, colMax = rowMax_lst, colMax_lst
    
    
    deckVis[slotKey] = np.full((rowMax, colMax)," ", dtype="U128")
    
    for wellKey in wells:
        well = slot[wellKey]
        
        wellItems = sorted(['{} : {:.2f}'.format(k,v) for k,v in well.items()])
        
        wellItems = "\n".join(wellItems)
        
        row = ord(wellKey[0])-65
        col = int(wellKey[1:])-1
        
        deckVis[slotKey][(row,col)] = wellItems
    
    deckTables[slotKey] = PrettyTable()
    deckTables[slotKey].field_names = ["({})".format(slotKey)] + list(np.arange(1,colMax+1,1))
    
    for rowKey in range(len(deckVis[slotKey])):
        tableRow = deckVis[slotKey][rowKey]
        tableRow = [x+"\n" for x in tableRow]
        tableRow = [chr(65+rowKey)] + tableRow
        
        
        deckTables[slotKey].add_row(tableRow)
