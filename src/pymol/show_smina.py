#pymol plugin that reads the atominfo file of smina and colors
#the passed selection appropriately 

from pymol import cmd
import re

#a class for (somewhat) efficiently matching atoms based on their position
#takes the parsed data of an atom terms file and can be used to lookup the
#energy of an atom based on its position
#this stores data for a single conformer
class MolPosLookup:
    def __makeinttuple(self, t):
        return tuple(map(lambda x: int((10*x)),  t))
    
    def __init__(self,data):
        self.data = dict()
        self.datalist = data
        for val in data:
            pos = val['pos']
            t = self.__makeinttuple(pos)
            self.data[t] = val['energy']
    
    def __distsq(selfself, a, b):
        return sum([(x[0]-x[1])*(x[0]-x[1]) for x in zip(a,b)])
    #find the data item with the closest coordinate to xyz
    def __findclosest(self, xyz):
        return min(self.datalist, key = lambda a:  self.__distsq(a['pos'],xyz))['energy']
    #return energy of atom with xyz position
    #return none if not close to known point
    def energyForPos(self, xyz):
        #currently very inefficient, should at least sort
        #along one dimension
        t = self.__makeinttuple(xyz)
        ret=self.data.get(t,None)
        if ret is None:
            #round off made us miss in the dictionary
            ret = self.__findclosest(xyz)
        return ret

#return a color in a gradient from green to white to red in pymol form
#for a little extra variation, throw yello into the green specturm
def gwr(val,min,max,mid):
    #cap value           
    if val < min:
        val = min
    elif val > max:
        val = max
    
    if val < mid: #green
        midmid = (min+mid)/2.0
        if(val < midmid): #green to yello
            scale = (midmid-val)/(midmid-min)
            return "[%f,1,0]" % (scale)
        else:
            scale = (val-midmid)/(mid-midmid)
            return "[1,1,%f]" % (scale)
    else:
        scale = 1-(val-mid)/(max-mid)
        return "[1,%f,%f]" % (scale,scale)
           
def show(sel, file,min=-1,max=1):
    f = open(file)
    header = f.readline()
    #read in atom term data
    statesdata = []
    data = []
    pregex = re.compile('<(.*),(.*),(.*)>')
    for line in f:        
        d = line.split()
        if len(d) > 3 and d[0] != "atomid":
            pos = d[2]
            m = pregex.match(pos)
            energy = sum(float(e) for e in d[3:])
            data.append({'pos':map(float,m.groups()), 'energy': energy})
        elif line.count("END"):
            #reset for next molecule
            statesdata.append(data);
            data = [];
    
    #now get atom data from pymol
    n_states = cmd.count_states(sel)
    if(n_states != len(statesdata)):
        print "Inconsistent number of states/molecules";
        return;
    
    for i in xrange(1,n_states+1):
        energies = MolPosLookup(statesdata[i-1])
        model = cmd.get_model(sel,i)
        for a in model.atom:
            e = energies.energyForPos(a.coord)
            if e is None:
                print "Missing",i,a.coord
            else:
                color = gwr(e,min,max,0)
                cmd.set_color("smina_color%d" % a.index,color)
                cmd.color("smina_color%d" %a.index,"%s and index %d" % (sel,a.index))
    

#register with pymol
cmd.extend('show_smina', show)
