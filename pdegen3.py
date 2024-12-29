import os
import sys
import re


pdelist = []
resolution = 0
class PDE:
    def __init__(self,name,arguments,right_hand):
        self.name = name
        self.rhs = right_hand
        self.arglist = []
        self.argmins = []
        self.argmaxs = []
        if arguments == "":
            pass
        else:
            args = arguments.split(",")
            for arg in args:
                self.arglist.append(arg.split()[0].strip())
                self.argmins.append(float(arg.split()[1].strip()))
                self.argmaxs.append(float(arg.split()[2].strip()))
        #The indices for the ODEs this gets turned into in the output script.
        self.first_index = 0
        self.last_index = 0

def parse_integrals(integral_line):
    #First, check if there are any integrals at all.
    #Especially important because we call this method recursively to deal
    #with possible nested or future integrals.
    if "Int(" not in integral_line:
        return integral_line
    #In order to make sure we track parenthesis balancing correctly, we
    #basically crawl forward from the "Int(", tracking opens and closes, until
    #we get the right ).
    start_index = integral_line.find("Int(")
    num_opens = 1
    end_index = start_index + 4
    while num_opens >= 1 and end_index < len(integral_line):
        if integral_line[end_index] == "(":
            num_opens += 1
        if integral_line[end_index] == ")":
            num_opens -= 1
        end_index += 1
    
    #If this executes, there weren't enough )s to balance the (s.
    if num_opens > 1:
        sys.stderr.write("Error: Unbalanced parentheses in PDE\n")
        sys.stderr.write(integral_line + "\n")
        sys.exit(5)
    
    #Cut the string into three pieces: integral itself and the bit before and
    #after the integral.  
    before_integral = integral_line[:start_index]
    after_integral = integral_line[end_index:]
    integral = integral_line[start_index:end_index]
    #The integral string takes the following form:
    #Int(start,end,expression,dummy,identity)
    #Expression might itself have commas in it.  The other four won't.
    #We split it accordingly.
    start = integral.split(",")[0][4:].strip()
    end = integral.split(",")[1].strip()
    temp = integral.split(",",2)[2]
    dummy = temp.rsplit(",",2)[1].strip()
    #Identity might be just varname instead of varname$index,
    #but this will work anyway.
    #Onle one of the $ or ) splits will actually do anything,
    #but I don't know in advance which it will be.
    identity = temp.rsplit(",",2)[2].strip().split("$")[0].split(")")[0]
    integrand = temp.rsplit(",",2)[0].strip()
    
    #Integrand itself might have integrals in it.  Parse them first.
    integrand = parse_integrals(integrand)
    
    #Find the absolute lower and upper bounds the identity variable could have.
    #This is used to calculate the delta for the Riemann sum.
    identity_lower = 0
    identity_upper = 0
    found = False
    for pde in pdelist:
        for i in range(len(pde.arglist)):
            if identity == pde.arglist[i]:
                found = True
                identity_lower = pde.argmins[i]
                identity_upper = pde.argmaxs[i]
    
    if not found:
        sys.stderr.write("Invalid dummy variable identity: " + identity + "\n")
        sys.exit(5)
    
    delta = (identity_upper - identity_lower)/float(resolution-1)
    
    #Find the indices that most closely match the bounds, if they're numerical.
    #If the bounds are not numerical, they'll be findable just by splitting "$"-wise
    start_index = 0
    end_index = 0
    if "$" in start:
        if start.split("$")[0] == identity:
            start_index = int(start.split("$")[1])
        else:
            sys.stderr.write("Bounds/identity mismatch in PDE:\n")
            sys.stderr.write(integral_line + "\n")
            sys.exit(5)
    else:
        startval = float(start)
        start_index = int(round((startval - identity_lower)/delta))
        if start_index < 0 or start_index > resolution:
            sys.stderr.write("Invalid integral lower bound in PDE:\n")
            sys.stderr.write(integral_line + "\n")
            sys.exit(5)
    if "$" in end:
        if end.split("$")[0] == identity:
            end_index = int(end.split("$")[1])
        else:
            sys.stderr.write("Bounds/identity mismatch in PDE:\n")
            sys.stderr.write(integral_line + "\n")
            sys.exit(5)
    else:
        endval = float(end)
        end_index = int(round((endval - identity_lower)/delta))
        if end_index < 0 or end_index > resolution:
            sys.stderr.write("Invalid integral upper bound in PDE:\n")
            sys.stderr.write(integral_line + "\n")
            sys.exit(5)
    expanded_integral = ""
    #if start_index = end_index, the integral's value is 0.
    if start_index == end_index:
        expanded_integral = "0"
    else:
        expanded_integral = str(delta)+ "*("
        #We add 1 to end_index because we actually do want that to be used
        #since this is being calculated by the trapezoidal rule
        for i in range(start_index,end_index+1):
            #Trapezoidal rule, multiply endpoints by 0.5
            if i == start_index or i == end_index:
                expanded_integral += "0.5*"
            expanded_integral += "(" + re.sub(r"((?<=\W)|^)" + dummy + "((?=\W)|$)",identity+"$"+str(i),integrand) + ")"
            if i != end_index:
                expanded_integral += "+"
        expanded_integral += ")"
        
    
    #Finally: the above will only have dealt with the first integral and any
    #integrals inside it.  Call the method again to deal with any later
    #integrals.
    return parse_integrals(before_integral + expanded_integral + after_integral)
    
    

if len(sys.argv) < 3:
    print("PDEgen usage: python pdegen3.py equation_file.txt outfile_name")
    print("Will produce a Python file outfile_name.py containing an implementation of")
    print("the PDE system described in equation_file.txt.")
    sys.exit(1)

infile = sys.argv[1]
if not os.path.isfile(infile):
    sys.stderr.write(infile + ": No such file found.\n")
    sys.exit(1)

outfile = sys.argv[2]
if os.path.isfile(outfile):
    sys.stderr.write(outfile + " already exists!  Aborting.\n")
    sys.exit(1)



#Initial parsing for the input file.
num_resolution_lines = 0
func_lines = []
partial_lines = []
initial_lines = []
maxtime = ""
num_maxtime_lines = 0
timestep = ""
num_timestep_lines = 0
for line in open(infile,'r').readlines():
    #Remove all comments.
    stripped_line = line.split("#")[0].strip()
    #Ignore empty lines.  Whitespace is allowed, and some lines might be only comments.
    if stripped_line == "":
        continue
    if stripped_line.split(" ")[0] == "resolution":
        num_resolution_lines += 1
        resolution = int(stripped_line.split(" ")[1])
        continue
    if stripped_line.split(" ")[0] == "func":
        func_lines.append(stripped_line[4:].strip())
        continue
    if stripped_line.split(" ")[0] == "partial":
        partial_lines.append(stripped_line[7:].strip())
        continue
    if stripped_line.split(" ")[0] == "initial":
        initial_lines.append(stripped_line[7:].strip())
        continue
    if stripped_line.split(" ")[0] == "maxtime":
        maxtime = stripped_line[7:].strip()
        num_maxtime_lines += 1
        continue
    if stripped_line.split(" ")[0] == "timestep":
        timestep = stripped_line[8:].strip()
        num_timestep_lines += 1
        continue
    print("Invalid line:")
    print(line)
    sys.exit(1)

if num_resolution_lines != 1:
    sys.stderr.write("Input file must have exactly one \"resolution\" keyword line.")
    sys.exit(2)


if num_maxtime_lines > 1:
    sys.stderr.write("Input file must have no more than one \"maxtime\" keyword line.")
    sys.exit(2)

if num_maxtime_lines == 0:
    maxtime = "50.0"

if num_timestep_lines > 1:
    sys.stderr.write("Input file must have no more than one \"timestep\" keyword line.")
    sys.exit(2)

if num_timestep_lines == 0:
    timestep = "0.1"

#Construct helper functions.
helper_functions = ""
for line in func_lines:
    helper_functions += "def " + line.split("=")[0].strip() + ":\n"
    helper_functions += "    return " + line.split("=",1)[1].strip()
    helper_functions += "\n\n\n"


#Parse the PDEs.
for line in partial_lines:
    pde_name = line.split("(")[0]
    if "$" in line:
        sys.stderr.write("Error: " + pde_name + " equation contains disallowed character $.\n")
        sys.exit(3)
        
    #Check for duplicate names.
    for p in pdelist:
        if pde_name == p.name:
            sys.stderr.write("Error: " + pde_name + " used multiple times as PDE function name.\n")
            sys.exit(3)
    pde_args = line.split("(")[1].split(")")[0]
    pde_rhs = line.split("=",1)[1].strip()
    pdelist.append(PDE(pde_name,pde_args,pde_rhs))


#Go down the list of PDEs.  For each PDE, we'll create resolution^dim ODEs,
#where dim is the number of independent variables (besides time) which the
#function uses.
equation_index = 0
filenames = []
eqnlines = []
eqnnames = []
eqnregexs = []
initials = ""
for pde in pdelist:
    pde.first_index = equation_index
    
    #We set up initial conditions at the same time as we parse the PDEs.
    initial_conditions = "0"
    for line in initial_lines:
        if line.find(pde.name+"(") == 0:
            initial_conditions = line.split("=",1)[1].strip()
            break
    
    
    #TODO: Replace varname variables in the ODE's arglist with varname$#, with #
    #indicating discretization slice number.
    #Example: if a PDE has independent variables called "irisk" and "iload",
    #both ranging 0 to 1, and resolution is set to 3, we will get 9 equations.
    #First one will replace irisk and iload with irisk$0 and iload$0.
    #Second will replace them with irisk$1 and iload$0.
    #Third will replace them with irisk$2 and iload$0.
    #Fourth will replace them with irisk$0 and iload$1.
    #And so on, unil we get to:
    #Ninth will replace them with irisk$2 and iload$2.
    #We leave them like this for now because there are two possible fates
    #for these strings.
    #First, if "irisk" was found as an argument to a PDE function, it will have
    #to be excised entirely and turned into the right function call.
    #We might have I(irisk$0,iload$0) ultimately become y[0].
    #Whereas I(irisk$0,iload$1) would become y[3].
    #(Indices higher than 0 and 3 if it's not the first PDE, of course.)
    #If it's just used as a number rather than an argument to a PDE function,
    #though, we instead need to replace it with its numerical value.
    #So irisk$0 becomes 0, and irisk$1 becomes 0.5.
    #XXX: We must do this BEFORE passing to parse_integrals.  parse_integrals
    #relies on this preprocessing work.
    
    #Dummy is never used, but I need to loop a specific number of times.
    for dummy in range(resolution**len(pde.arglist)):
        #Go has a feature where foreach loops get both the token and index.
        #I wish Python did too.
        raw_index = equation_index - pde.first_index
        index_list = []
        parsed_eqn = pde.rhs
        parsed_initials = initial_conditions
        for i in range(len(pde.arglist)):
            index = int(raw_index / (resolution**i)) % resolution
            indexed_var = pde.arglist[i] + "$" + str(index)
            index_list.append(index)
            
            #Replace all instances of arglist[i] with indexed_var
            #Can't just do this with a basic string replacement, though, might hit
            #substrings.
            parsed_eqn = re.sub(r"((?<=\W)|^)" + pde.arglist[i] + "((?=\W)|$)",indexed_var,parsed_eqn)
            parsed_initials = re.sub(r"((?<=\W)|^)" + pde.arglist[i] + "((?=\W)|$)",indexed_var,parsed_initials)
        
            
        initials += parsed_initials + ","
        
        eqnlines.append(parse_integrals(parsed_eqn))
        
        #Figure out what the file name will be for this ODE file output
        eqnname = pde.name
        for i in range(len(pde.arglist)):
            stepsize = (pde.argmaxs[i] - pde.argmins[i]) / (resolution-1)
            eqnname += "_" + pde.arglist[i] + "=" + str(index_list[i] * stepsize + pde.argmins[i])
        eqnnames.append(eqnname)
    
        #Determine a regular expression that can be used as a replacement target for
        #the next parsing pass.
        eqnregex = r"((?<=\W)|^)" + pde.name + "\("
        for i in range(len(pde.arglist)):
            eqnregex += " *" + pde.arglist[i] + "\$" + str(index_list[i]) + " *"
            if i < len(pde.arglist) - 1:
                eqnregex += "," #All but the last argument will have following ,
        eqnregex += "\)"
        eqnregexs.append(eqnregex)
        equation_index += 1
        



    pde.last_index = equation_index - 1

    #As it currently stands, the things added to eqnlines still have unparsed
    #calls to PDE functions.  Can't be helped, need to do that in another pass
    #afterward.  And after that pass we do the number replacements.

        


#Parse things like funcname(foo$index,bar$index) into things like y[num]
for i in range(len(eqnlines)):
    for j in range(len(eqnregexs)):
        eqnlines[i] = re.sub(eqnregexs[j],"y["+str(j)+"]",eqnlines[i])

#Now, parse any remaining foo$index things into numbers.
#This is a stupid and inefficient way to do it, but it will do for now.
for i in range(len(eqnlines)):
    for pde in pdelist:
        for j in range(len(pde.arglist)):
            for k in range(resolution):
                value = str(((pde.argmaxs[j] - pde.argmins[j]) / (resolution-1)) * k)
                eqnlines[i] = re.sub(r"((?<=\W)|^)" + pde.arglist[j]+"\$"+str(k)+"((?=\W)|$)",value,eqnlines[i])


#It's possible for a PDE function to be called in the source file with another
#PDE's arguments.  For example, we could have a situation like:
#first(foo 0 1) = stuff
#second(bar 0 1) = first(bar)
#In that case, at this point the ODE expansion of second's PDE will still have
#unparsed "first" calls in it, because the first parsing pass wouldn't catch
#it with eqnregexs because of the argument mismatch, and the second only
#replaces variables with numbers.  Here we fix this: any function calls
#remaining will have their input numbers parsed to the nearest applicable
#variable indices and then we run the regexes again.
#FIXME: Currently will fail for arguments that aren't just variables.
for i in range(len(eqnlines)):
    for pde in pdelist:
        searchresult = re.search(r"((?<=\W)|^)" + pde.name + "\(",eqnlines[i])
        while searchresult:
            #FIXME: We assume there are no parentheses in function calls.
            #This will need to be fixed if this is to be a production tool,
            #although of course probably this entire implementation will need
            #to be scrapped in that case.
            startindex = searchresult.end()
            endindex = eqnlines[i].find(")",startindex)
            
            before = eqnlines[i][:startindex]
            after = eqnlines[i][endindex:]
            nums = eqnlines[i][startindex:endindex]
            
            numlist = nums.split(",")
            indexed_list = []
            for j in range(len(numlist)):
                value = float(numlist[j])
                index = str(int(round(value / ((pde.argmaxs[j]-pde.argmins[j])/(resolution-1)))))
                indexed_list.append(pde.arglist[j] + "$" + index)
            
            eqnlines[i] = before + ",".join(indexed_list) + after
            #Have to do it this because re.search doesn't take a start index, whereas
            #regexname.search does.  For some reason.
            pattern = re.compile("((?<=\W)|^)" + pde.name + "\(")
            searchresult = pattern.search(eqnlines[i],endindex)
            #searchresult = re.search(r"((?<=\W)|^)" + pde.name + "\(",eqnlines[i])


#Parse things like funcname(foo$index,bar$index) into things like y[num].
#Again.
for i in range(len(eqnlines)):
    for j in range(len(eqnregexs)):
        eqnlines[i] = re.sub(eqnregexs[j],"y["+str(j)+"]",eqnlines[i])
                
            


#Shave off the trailing comma, because everything appended to initials added
#a comma
initials = initials[:-1]
for pde in pdelist:
    for j in range(len(pde.arglist)):
        for k in range(resolution):
            value = str(((pde.argmaxs[j] - pde.argmins[j]) / (resolution-1)) * k)
            initials = re.sub(r"((?<=\W)|^)" + pde.arglist[j]+"\$"+str(k)+"((?=\W)|$)",value,initials)



#Parsing is FINALLY done.  Now time to make the output file.
#FIXME: dirname should either be an argument or something in the source file
output = """#This file was automatically generated by pdegen.py.
#It is not meant to be edited, unless you're sure what you're doing.
#The equations are probably not very readable.

from scipy.integrate import odeint
from numpy import arange
import time
import datetime
import os
import sys

dirname = \"\"

"""

output += helper_functions

output += """

#Initial conditions
y0 = ["""

output += initials

output += """]


names = []
"""

for name in eqnnames:
    output += "names.append(\"" + name + "\")\n"

output += """


#The derivative function for the differential equation system.
def func(y,t):
    return [
"""

for line in eqnlines:
    output += "             " + line + ",\n"

output += """           ]


t = arange(0, """ + maxtime + ", " + timestep + """)

y = odeint(func, y0, t, ixpr=True)

if dirname == \"\":
    dirname = datetime.datetime.fromtimestamp(time.time()).strftime(\'%Y-%m-%d-%H:%M:%S\')
os.makedirs(dirname)
os.chdir(dirname)

for i in range(len(y0)):
    writer = open(names[i]+\".txt\", \'w\')
    for j in range(len(t)):
        writer.write(str(t[j]) + \" \" + str(y[j][i]) + \"\\n\")
    writer.close()
"""

#FIXME: Both max time and timestep should be adjustable in the source file.

writer = open(outfile,'w')
writer.write(output)
