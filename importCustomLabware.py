import json
import opentrons.execute

filePath = "C:/Users/tpeng/OneDrive/OD__Universität/BA__Axmann/FlowerPlate 48 Well Plate 3200 µL.json"

protocol = opentrons.execute.get_protocol_api()
with open(filePath) as labware_file:
    labware_def = json.load(labware_file)
    
well_plate = protocol.load_labware_from_definition(labware_def, 1)