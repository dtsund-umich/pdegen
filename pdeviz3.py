import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib.widgets import Slider, Button, RadioButtons
from mpl_toolkits.mplot3d import Axes3D


#I think how I will do it is this:
#The first time I get a function name, read all the associated functions files.
#Make an array, graphdatas, each entry of which is the necessary graph data
#at the appropriate timestep.
#For functions that are just ODE functions, this will be a number.
#For functions that are PDE functions of one non-time input, a line graph.
#For PDE functions with two non-time inputs, a 3D landscape.
#For PDE functions with higher numbers of input variables...
#I don't know yet.
#Maybe some sort of two-at-a-time deal with sliderbars underneath?
#Will need to figure out an interface for that.
def sortkey(value,index):
    #Trim off the .txt from the end and sort by whichever key value we want
    return float(value[:-4].split("_")[index].split("=")[1])

class pdeoutput:
    def __init__(self, fname, farglist, allfiles, numtimepoints):
        self.name = fname
        self.argnames = []
        for token in farglist:
            self.argnames.append(token.split("=")[0])
        
        #Identify the names of the relevant files.
        relevant_files = []
        for f in allfiles:
            prefix = f[:-4]
            tokens = prefix.split("_")
            #It's possible for an equation name to have an underscore in it.
            #If so, the name will be spread across multiple tokens.
            #But it won't have an = in it, so we just see which tokens don't have them.
            #And the tokens with = will always come last.
            tindex = 0
            
            for token in tokens:
                if "=" in token:
                    break
                tindex += 1
            
            if eqnname == "_".join(tokens[0:tindex]) == fname:
                relevant_files.append(f)


        self.argmins = []
        self.argmaxs = []                
        #Some of the innards of this for loop would break badly on an ODE
        #function, but it should never execute for those.
        for i in range(len(self.argnames)):
            minval = float(farglist[i].split("=")[1])
            maxval = float(farglist[i].split("=")[1])
            for f in relevant_files:
                if float(f[len(fname)+1:-4].split("_")[i].split("=")[1]) < minval:
                    minval = float(f[len(fname)+1:-4].split("_")[i].split("=")[1])
                if float(f[len(fname)+1:-4].split("_")[i].split("=")[1]) > maxval:
                    maxval = float(f[len(fname)+1:-4].split("_")[i].split("=")[1])
            self.argmins.append(minval)
            self.argmaxs.append(maxval)
        
        #Sort the list of relevant files so we can conveniently read the
        #data we want in the order we want it.
        for i in range(fname.count("_")+1,fname.count("_")+len(farglist)+1):
            relevant_files.sort(key=lambda x: sortkey(x,i))
        
        #How many data points "wide" is each input?
        width = 0
        if len(farglist) > 0:
            width = int(round(len(relevant_files)**(1./float(len(farglist)))))

        #Find the increments between function inputs.
        self.argsteps = []
        for i in range(len(self.argnames)):
            self.argsteps.append((self.argmaxs[i] - self.argmins[i]) / (width-1))
        
        
        #Build the vertical coordinate data array.
        self.graphdatas = []
        for i in range(numtimepoints):
            rawlist = []
            #Open each file, read all the lines, take the thing we want from
            #the lines we want.
            #This is extremely inefficient, but compact.
            #If performance is a serious issue, I'll come back and rewrite
            #this part.
            for f in relevant_files:
                rawlist.append(float(open(f,"r").readlines()[i].split()[1].strip()))
            #Data is in a line.  Make it be in a square, cube, etc. as needed.
            
            #FIXME: Only works for 0,1,2 dimensions for now.
            if len(farglist) == 0:
                self.graphdatas.append(rawlist[0])
            if len(farglist) == 1:
                self.graphdatas.append(rawlist)
            if len(farglist) == 2:
                self.graphdatas.append(rawlist)
                #temp = []
                #for j in range(width):
                #    temp.append(rawlist[j*width:(j+1)*width])
                #self.graphdatas.append(temp)
        
        
        
        self.fixedgraphdata = []
        #Build the array describing the X (2D) or XY (3D) axes.
        if len(farglist) == 0:
            templist = []
            for line in open(relevant_files[0],"r").readlines():
                templist.append(float(line.split()[0].strip()))
            self.fixedgraphdata = templist
        if len(farglist) == 1:
            self.fixedgraphdata = list(np.linspace(self.argmins[0],self.argmaxs[0],width))
        if len(farglist) == 2:
            #X data
            self.fixedgraphdata.append(list(np.linspace(self.argmins[0],self.argmaxs[0],width))*width)
            #Y data
            self.fixedgraphdata.append(list(np.repeat(np.linspace(self.argmins[1],self.argmaxs[1],width),width)))
            #ydata = list(np.linspace(self.argmins[1],self.argmaxs[1],width))*width
            #ydata = list(np.transpose(np.array(ydata)))
            #self.fixedgraphdata.append(ydata)
        
        #With the data thus built, a 3D plot should be generatable by doing
        #something like
        #ax.plotwireframe(self.fixedgraphdata[0],self.fixedgraphdata[1],graphdatas[t])
        #where t represents the timestep.
                
        
    def plot_data(self,axes,t):
        axes.clear()
        axes.set_title(self.name)
        if len(self.argnames) == 0:
            axes.plot(self.fixedgraphdata,self.graphdatas)
            axes.set_xlabel("Time")
        if len(self.argnames) == 1:
            axes.plot(self.fixedgraphdata,self.graphdatas[t])
            axes.set_xlabel(self.argnames[0])
        if len(self.argnames) == 2:
            squaresize = int(round(len(self.fixedgraphdata[0])**0.5))
            axes.plot_wireframe(np.reshape(self.fixedgraphdata[0],(-1,squaresize)),
                                np.reshape(self.fixedgraphdata[1],(-1,squaresize)),
                                np.reshape(self.graphdatas[t],(-1,squaresize)))
            axes.set_xlabel(self.argnames[0])
            axes.set_ylabel(self.argnames[1])
        
        
        
        

if len(sys.argv) < 2:
    print("PDEviz usage: python pdeviz.py output_directory")
    print("The argument output_directory should be a directory created by a")
    print("PDE implementation, which in turn should have come from pdegen.py.")
    sys.exit(1)

if not os.path.isdir(sys.argv[1]):
    print("Error: " + sys.argv[1] + " not a directory!")
    sys.exit(1)

os.chdir(sys.argv[1])

filenames = os.listdir(".")

#Figure out the timecourse min, max, step
firstlines = open(filenames[0],"r").readlines()
tmin = 0
tstep = float(firstlines[1].split()[0]) - float(firstlines[0].split()[0])
tmax = float(firstlines[-1].split()[0])

#Read all the files.  Sort the data in a useful manner.
pdes = []
for f in filenames:
    prefix = f[:-4]
    tokens = prefix.split("_")
    #It's possible for an equation name to have an underscore in it.
    #If so, the name will be spread across multiple tokens.
    #But it won't have an = in it, so we just see which tokens don't have them.
    #And the tokens with = will always come last.
    tindex = 0
    for token in tokens:
        if "=" in token:
            break
        tindex += 1
    
    eqnname = "_".join(tokens[0:tindex])
    
    #Is this the first time we've run into this function?
    found = False
    for pde in pdes:
        if pde.name == eqnname:
            found = True
            break
    
    if not found:
        #tokens[tindex:] information is technically redundant with filenames,
        #but makes things a little easier.
        pdes.append(pdeoutput(eqnname, tokens[tindex:], filenames,int(round(tmax/tstep))+1))
        
    

#Figure out the dimensions of the grid of axes.  Number of columns should be
#the square root of the number of graphs we need (rounded up); number of rows
#should be as many rows as needed given that number of columns.
numcols = int(math.ceil(len(pdes)**0.5))
numrows = int(math.ceil(len(pdes)/numcols))
#plt.subplots(x,y) for an x-by-y array of figures
fig, ax = plt.subplots(numrows,numcols)
plt.tight_layout()

#Add a clickable slider bar at the bottom of the figure that we use to choose
#timestamps to look at.
plt.subplots_adjust(bottom=0.25)
axtime = plt.axes([0.25, 0.1, 0.65, 0.03])
time = Slider(axtime, "Time", tstep, tmax, valinit = 0)

#At this point, we know how many "real" sets of axes there are in the graph.
#For 3D plots, what we do is set the corresponding 2D plot to be invisible and
#add a 3D plot in its place.  But we add a new 3D plot each time we update the
#timestamp, so we need to delete the old one each time we do or we leak memory.
#To that end: num_true_axes stores the index of the last set of axes we've added
#so far.  Anything with a higher index is a 3D plot that we should discard on
#each update.
num_true_axes = len(fig.get_axes())

#First have an initial for loop that sets up the axes as necessary.
#The thing defaults to making everything 2D axes, so 3D plots will fail
#unless you specifically set each axes type.
for i in range(numrows*numcols):
    if i >= len(pdes):
        #A set of axes with no corresponding data to graph.  Make it invisible
        #for aesthetic and readability purposes.
        ax[i//numcols][i%numcols].axis('off')
    elif len(pdes[i].argnames) == 2:
        ax[i//numcols][i%numcols].axis('off')
        pdes[i].plot_data(fig.add_subplot(numrows,numcols,i+1,projection='3d'),0)
    else:
        pdes[i].plot_data(ax[i//numcols][i%numcols],0)



def update(val):
    #First, delete the old 3D graphs.
    fakeaxes = fig.get_axes()[num_true_axes:]
    for fake in fakeaxes:
        fig.delaxes(fake)
    
    #Now, update the data in the existing graphs (and in the case of 3D graphs,
    #construct new plots to replace the ones we just nuked).
    t = int(time.val / tstep)
    for i in range(len(pdes)):
        if len(pdes[i].argnames) == 2:
            pdes[i].plot_data(fig.add_subplot(numrows,numcols,i+1,projection='3d'),t)
        else:
            pdes[i].plot_data(ax[i//numcols][i%numcols],t)
    plt.draw()

time.on_changed(update)



plt.show()
