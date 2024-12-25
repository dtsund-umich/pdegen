#Generate total population over time for output from pdegen systems, for debugging
import os
import sys

if len(sys.argv) < 2:
    print("Give me a directory of output made by a pdegen system of ODEs")


os.chdir(sys.argv[1])

filenames = os.listdir(".")

totalpop = []

timestamps = []

popdict = {}

for f in filenames:
    flines = open(f,'r').readlines()
    if timestamps == []:
        for line in flines:
            timestamps.append(line.split()[0].strip())
        totalpop = [0]*len(timestamps)

    multiplier = 1
    if "_" in f:
        arglist = f.strip()[:-4].split("_")[1:]
        name = f.strip().split("_")[0]
        if not name in popdict:
            popdict[name] = [0]*len(timestamps)
        
        for arg in arglist:
            argname = arg.split("=")[0]
            numvalue = float(arg.split("=")[1])
            allvalues = []
            for g in filenames:
                if g.strip().split("_")[0] == name:
                    g_arglist = g.strip()[:-4].split("_")[1:]
                    for g_arg in g_arglist:
                        if g_arg.split("=")[0] == argname:
                            allvalues.append(float(g_arg.split("=")[1]))
            
            
            allvalues.sort()
            
            if numvalue == min(allvalues) or numvalue == max(allvalues):
                multiplier *= 0.5
            
            multiplier *= (max(allvalues) - min(allvalues))/(len(allvalues) - 1)
            
    print(f + ": multiplier="+str(multiplier))
    for i in range(len(flines)):
        totalpop[i] += float(flines[i].strip().split(" ")[1]) * multiplier
        if "_" in f:
            popdict[name][i] += float(flines[i].strip().split(" ")[1]) * multiplier
    print(totalpop[0])




writer = open("totalpop.txt",'w')
for i in range(len(timestamps)):
    writer.write(timestamps[i] + " " + str(totalpop[i]) + "\n")
writer.close()

for name in popdict:
    writer = open("total_"+name+".txt",'w')
    for i in range(len(timestamps)):
        writer.write(timestamps[i] + " " + str(popdict[name][i]) + "\n")
    writer.close()

