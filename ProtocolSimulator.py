from opentrons import simulate, robot
filepath = "tpeng/OneDrive/OD__Universität/BA__Axmann/PlateScripts/"
outpath = "tpeng/OneDrive/OD__Universität/BA__Axmann/protocol_simulations/"

filename = "P03_OTV1_GlcGradientONECELL_BigFalconRack"

PC = False
print_output = False


if PC:
    header = "D:/"
else:
    header = "C:/Users/"

protocol_file = open(header+filepath+filename+".py")

robot.reset()    
runlog = simulate.simulate(protocol_file, file_name = filename+".py")

protocol_file.close()


if print_output:
    print(simulate.format_runlog(runlog[0]))
else:
    file = open(header+outpath+filename+"_out.txt", "w")
    file.write(simulate.format_runlog(runlog[0]))
    file.close()
   
robot.reset()