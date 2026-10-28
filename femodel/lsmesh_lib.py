# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 16:07:48 2015

@author: Konstantinov
v2.0 - Calculation of mass properties: mass, volume,
       center of mass and inertia tensor.
v2.1 - Calc mass props fix.
v2.2 - Fixing inertia matrix to central.
v2.3 - Fixing xy, xz, yz inertia values.
v3.0 - Universal sweeping algorithm for
       any 2D mesh.
v3.1 - Sweep element numberring fix.
v3.2 - Fixed usage of include files.
       Added list of part numbers as
       a property of model object.
v3.3 - Sweep segments.
v3.4 - Sweep node set.
v4.0 - New parsing method. readSegmentSets method for model.
v4.1 - Added model reflect method.
v4.5 - Added scale, translate and rotate for model.
v5.0 - Storing and using etype property of element.
       Added 2D extrusion procedure.
v5.1 - Implemented volume calculation (and inertia properties)
       calculation for axisymmetric shells.
v5.2 - Fixed model add operation.
v5.3 - fixingSegmentsOrientationMethod.
v6.0 - New universal sweeping algorithm.
v6.1 - Swepped segments orientation fixed.
v6.2 - readSegmentSet method fixed.
       Added saveSegmentSet method.
v6.3 - Added BoundBox method for model.
v6.4 - Added readMSH function.
v6.5 - Added saveCCXmesh method. Added some functions.
v6.6 - Added element containsNds and getFaceNum methods.
v6.7 - Added getFacesNumbers method for model class.
v6.8 - Added offsets to saveCCXmesh procedure. Fixed el.nodes property.
v6.9 - Added nShiftRequest argument to sweep2Dmesh procedure.
v7.0 - msh second order elements.
v7.1 - Fixed second ordet tet reading from gmsh format.
v7.2 - nds2Dcount property added to sweep procedure.
v8.0 - new getElementCenterByNodes added. Improved 2D mesh modification on sweeping.
v9.0 - added procedures for gmsh meshing.
v9.2 - adding cleanDegenerated property to model. New gmsh algorithm for gmsh
       quad dominant meshing.
v9.3 - Can choose meshing algorithm in mesh2DquadDominant procedure.
v9.4 - Added projection to z=0 in procedure fix2DMeshForSweeping.
v9.5 - Added getShellNormal method.
v9.6 - Added getNodesOnLines method to model class.
v9.7 - Changed result format for sweep2Dmesh procedure.
v9.8 - fix2Dorientation for high order triangles. gmshMesh general procedure.
v9.9 - Console progress for getFacesNumbers method.
v10.0 - model.uniteParts method.
v10.1 - model.uniteParts by element set method.
v11   - Element quality for 3D bricks (edge lengths and edge angles) -
        getBadElements model method.
        formatSetKeyword function.
v11.1 - Element quality for all solid element types.
v11.2 - New algorithm in fixSegmentsOrientations model method.
v11.3 - New functions for fast searching in list isInSet and isSubSet.
v11.4 - element.getFaceNum method modification.
v12.0 - Parallel element find faces method.
v12.1 - Parallel element find faces method commented.
        3D elements getFaces fixed to external normal.
        Optimization of fixSegmentsOrientations method.
v12.2 - setInd bug fix.
        del_node, nodes, element, elements speed improvement.
v12.3 - Fixed BoundBox model method.
v12.4 - Fixed _nodes, _getNodes and uniqNodes properties of element object.
v12.5 - Fixed getFaces for hight order elements.
v12.6 - getNearestNode model method.
        distTo node method.
v12.6 - getNearestNode model method. Optimized get_element_center model method.
v12.7 - model getExternalNodes method.
v12.8 - fixed model __add__ method.
v12.9 - model __add__ method (+) use numbering shift, __mul__ (*) method just
        unite two models without changing numeration.
v13   - plane box mesh procedure
v14   - added from __future__ import division
v14.1 - Saving model with empty nodes fix
v14.2 - class boundBox, node.isInsideBox, model.nodesInsideBox
        added lsPrePostTools
v14.3 - model.findUnconnectedRegions. model.findNodePairsForTieing
v14.4 - writeSegments fix
v14.5 - writeSegments fileStream argument support
v15.0 - lsdyna_model class
v15.1 - gmsh Tet mesh element size based on curvature.
v15.2 - gmsh Tet mesh element size based on curvature - optional. Default - False
v15.3 - externalSegmentsForAleFill procedure
v15.6 - fixed getFaces for pyramid elements
v15.7 - separeted segments and nodes lists in mesh2Dsweep (use set of lines
        instead list)
v15.8 - fixed read_set function
06.05.19 - extrude2dmesh ndset and segments list and dict
06.05.19 - model.getExternalNodes method new implementation (fast)
          (numpy >1.13.0 requared)
06.05.19 - extrude2dmesh ndset and segments list and dict
08.05.19 - fixed model.getExternalNodes method
13.05.2019 - model.getMeshStatistics mehod, point.angle fix
31.05.2019 - model.nodes_to_nprecord method
05.06.2019 - list of element angles can be passed to sweep2dmesh:
             for example: fdi = [5,5,10,5,5]
           - list of element heights can be passed to extrude2dmesh,
             i.e. depth = [0.1, 0.1, 0.2, 0.1]
10.06.2019 - model.findDuplicateNodes
           - model.sewDuplicateNodes
           - default sew nodes when meshing 2D quad dominant
           - fixing elements with 3 nodes on axis when prepearing mesh for sweeping.
13.06.2019 - fixed model.sewDuplicateNodes (model.allPoints renumbering)
28.06.2019 - sewDuplicateNodes fix
15.09.2019 - sortNodesForPlaneElement function, create_cohesive_layer
28.09.2019 - sort_8_nodes_for_hex
09.10.2024 - find model.findUnconnectedRegions x10 speedup
"""

from __future__ import print_function, generators, division, with_statement
import numpy as np
import math
import os
from math import sqrt, sin, cos, pi, acos
from copy import deepcopy
from bisect import bisect_left
#import pp
from os import system, remove, curdir, listdir
from collections import defaultdict
from itertools import chain, combinations, permutations
from functools import partial

lsPrePostPath='prepost'

class lsPrePost(object):
    def __init__(self, prepostPath=lsPrePostPath):
        self.prePost=prepostPath
    def createExternalNodesSetCommand(self, nodeSet=1):
        return """setnode
genselect target node
genselect clear
genselect 3dsurf on
genselect allvis
setnode createset {0:d} 1 0 0 0 0
genselect clear
""".format(nodeSet)

    def createExternalNodesSetByPartCommand(self, nodeSets=1, parts=1):
        if type(nodeSets)==int: nodeSets=[nodeSets]
        if type(parts)==int: parts=[parts]
        s="""setnode
genselect target node
genselect 3dsurf on
"""
        for ns, p in zip(nodeSets, parts):
            s+="""genselect clear
genselect node add part {0:d}/0
setnode createset {1:d} 1 0 0 0 0
""".format(p, ns)
        return s

    def openKeywordCommand(self, filePath):
        return "open keyword \"{0:s}\"\n".format(filePath)

    def saveKeywordCommand(self, filePath):
        return """save keywordbylongfmt 0
save keyword \"{0:s}\"
""".format(filePath)

    def delElementsInsideBoxCommand(self, p1, p2):
        return """delelement unrefn 2
delelement target element
delelement clean 1
genselect element add box in {0:10.3e} {1:10.3e} {2:10.3e} {3:10.3e} {4:10.3e} {5:10.3e}
delelement delete
delelement accept
genselect clear
""".format(*(p1+p2))

    def createNodeSetInsideBoxCommand(self, p1, p2, setNum=1):
        return """setnode
genselect target node
genselect clear
genselect node add box in {1:10.3e} {2:10.3e} {3:10.3e} {4:10.3e} {5:10.3e} {6:10.3e}
setnode createset {0:d} 1 0 0 0 0
genselect clear
""".format(setNum, *(p1+p2))

    def runCommand(self, command):
        f=open('tmp.cfile', 'w')
        f.write(command)
        f.close()
        cl=self.prePost+' -nographics c=tmp.cfile'
        print(cl)
        system(cl)
        remove('tmp.cfile')

    def extractDataFromNodout(self, nodoutPath, rezPath, nodes='all', components=[1]):
        command="""ascii nodout open "{0:s}"
ascii nodout plot {1:s} {2:s}
xyplot 1 savefile ms_csv "{3:s}" 1 all
"""
        cmps='/'.join(map(str, components))
        if nodes=='all':
            nds='all'
        else:
            nds='/'.join(map(str, nodes))
        runCommand(command.format(nodoutPath, cmps, nds, rezPath))

def shift(l, idx):
    return l[idx:]+l[:idx]

def splitByN(l, n):
    if n>=len(l):
        return [l]
    nn=len(l)//n
    rez=[]
    for i in range(nn):
        rez.append(l[i*n:(i+1)*n])
    if len(l)%n:
        rez.append(l[nn*n:])
    return rez

def splitByNseq(s, n):
    while sum(n)+n[-1]<=len(s):
        n.append(n[-1])
    pos=[sum(n[:i+1]) for i in range(len(n))]
    pos.insert(0,0)
    rez=[]
    for i in range(len(pos)-1):
        rez.append(s[pos[i]:pos[i+1]].strip())
    return rez

def setInd(n, s, isSorted=False):
    ss=s
    if not isSorted:
        ss=deepcopy(s)
        ss.sort()
    idx=bisect_left(ss, n)
    if idx==len(ss):
        idx=-1
    if idx==0 and ss[0]!=n:
        idx=-1
    if s[idx]!=n:
        idx=-1
    return idx

def isInSet(n, s, isSorted=False):
    return setInd(n, s, isSorted)!=-1

def isSubSet(s1, s2, isSorted=False):
    ss=s2
    if not isSorted:
        ss=deepcopy(s2)
        ss.sort()
    rez=True
    for n in s1:
        rez=rez and (isInSet(n, ss, isSorted=True)!=-1)
    return rez


def formatSetKeyword(data, num=1, kw='*SET_NODE_LIST'):
    n=len(data)
    if not n: return ""
    rez=kw+'\n{:10d}\n'.format(num)
    fs=("{:>10}"*8+"\n")*(n//8)
    if n%8:
        fs+="{:>10}"*(n%8)+"\n"
    rez+=fs.format(*data)
    return rez

class ConsoleProgress(object):
    def __init__(self, total, scaleLength=50, text='Progress=',
                 symbol='-#'):
        self.current = 0
        self.total = total
        self.chekRedraw = lambda i: not (i % (total//scaleLength))
        self.scaleLength = scaleLength
        self.text = text
        self.symbol = symbol

    def redraw(self):
        if self.chekRedraw(self.current):
            print(progress(self.current, self.total, self.scaleLength, self.text, self.symbol), end='')#, flush=True)
        self.current+=1

    def finalize(self):
        print(progress(self.total, self.total, self.scaleLength, self.text, self.symbol))


def progress(cur, total, scaleLength=20, text='Progress=', symbol='-#'):
    if total!=0:
        prog=int(100.*cur/total)
    else:
        prog=100
    rez='\r{0}{1:4}% '.format(text, prog)
    numD=int(prog*scaleLength/100.)
    rez+=symbol[1]*numD+symbol[0]*(scaleLength-numD)
    rez+=' Total={0:10}, current={1:10}'.format(total,cur)
    return rez

class Point(object):
    def __init__(self, x, y, z):
        """
        Point(x,y,z) - create point(vector) in 3D
        """
        self.x = x
        self.y = y
        self.z = z

    def norm(self):
        """
        Length of the vector
        """
        return sqrt(self.x*self.x+self.y*self.y+self.z*self.z)

    def normalized(self):
        """
        Return vector of unity length
        """
        l = self.norm()
        if l == 0:
            return self
        p = Point(self.x/l, self.y/l, self.z/l)
        return p

    def __repr__(self):
        return '(%f,%f,%f)' % (self.x, self.y, self.z)

    def dot(self, p2):
        rez = self.x*p2.x+self.y*p2.y+self.z*p2.z
        return rez

    def cross(self, p2):
        x = self.y*p2.z-self.z*p2.y
        y = self.z*p2.x-self.x*p2.z
        z = self.x*p2.y-self.y*p2.x
        return Point(x, y, z)

    def scaled(self, a):
        return Point(self.x*a, self.y*a, self.z*a)

    def __add__(self, p2):
        return Point(self.x+p2.x, self.y+p2.y, self.z+p2.z)

    def __sub__(self, p2):
        return Point(self.x-p2.x, self.y-p2.y, self.z-p2.z)

    def __neg__(self):
        return Point(-self.x, -self.y, -self.z)

    def shuff(self, p1, p2):
        return (self.cross(p1)).dot(p2)

    def _getData(self):
        return np.array([self.x, self.y, self.z])

    def _setData(self, coords):
        coords = list(coords) + [0]*(3-len(coords))
        self.x, self.y, self.z = coords

    data = property(_getData, _setData)

    def angle(self, p2, direction_vector=None):
        p1_u = self.normalized()
        p2_u = p2.normalized()
        rez = np.arccos(
            np.clip(np.dot(p1_u.data, p2_u.data), -1.0, 1.0))*180./np.pi
        if direction_vector:
            if self.cross(p2).dot(direction_vector) < 0:
                rez = 360-rez
        return rez
class Quaternion(object):
    data=[0,0,0,0]
    R=[[0,0,0],[0,0,0],[0,0,0]]
    def __init__(self, angle, v):
        """
        Quaternion(angle, vector)
        Creates quanternion with angle - angle, abount axis - vector.
        """
        self.a=angle
        self.v=v.normalized()
        a2=angle/2.
        self.data[0]=cos(a2)
        sina2=sin(a2)
        d=self.v.data
        for i in range(3):
            self.data[i+1]=d[i]*sina2
        self.calculateRotationMatrix()
    def calculateRotationMatrix(self):
        w=self.data[0]
        x=self.data[1]
        y=self.data[2]
        z=self.data[3]
        self.R[0][0]=1-2*y*y-2*z*z
        self.R[0][1]=2*x*y-2*z*w
        self.R[0][2]=2*x*z+2*y*w
        self.R[1][0]=2*x*y+2*z*w
        self.R[1][1]=1-2*x*x-2*z*z
        self.R[1][2]=2*y*z-2*x*w
        self.R[2][0]=2*x*z-2*y*w
        self.R[2][1]=2*y*z+2*x*w
        self.R[2][2]=1-2*x*x-2*y*y
    def rotatePoint(self, p):
        rez=[0,0,0]
        d=p.data
        for i in range(3):
            rez[i]=0
            for j in range(3):
                rez[i]+=self.R[i][j]*d[j]
        return Point(rez[0],rez[1],rez[2])

    def _getData(self):
        return [self.data[i] for i in range(4)]
    quaterionData=property(_getData)

    def _getRotationMatrix(self):
        rez=[]
        for i in range(3):
            rez.append([])
            for j in range(3):
                rez[i].append(self.R[i][j])
        return rez
    rotationMatrix=property(_getRotationMatrix)

def searchKeywords(data):
    comments=[]
    for i, l in enumerate(data):
        if l.startswith('$'):
            comments.append(i)
    comments.reverse()
    for i in comments:
        data.pop(i)
    keywords=[]
    indexes=[]
    for i, l in enumerate(data):
        if l.startswith('*'):
            keywords.append(l.upper()[:-1])
            indexes.append(i)
    return keywords, indexes

def readData(fname, procedIncludes=False):
    data=open(fname,'r').readlines()
    dirName=os.path.dirname(os.path.realpath(fname))+os.path.sep
    keywords, indexes = searchKeywords(data)
    if procedIncludes and keywords.count('*INCLUDE'):
        incFiles=[]
        for i, kw in enumerate(keywords):
            if kw=='*INCLUDE':
                incFiles.extend(data[indexes[i]+1:indexes[i+1]])
        for f in incFiles:
            print('Proced includes ', f[:-1])
            data.extend(readData(dirName+f[:-1], True))
    return data

def calcV4points(points):
    p=[np.array(pp) for pp in points]
    m=np.matrix([p[1]-p[0], p[2]-p[0], p[3]-p[0]])
    return abs(np.linalg.det(m)/6.)

def read_set(fname, stype='NODE', num=None):
    if isinstance(num, int):
        num=[num]
    rez={}
    f=open(fname, 'r')
    l=next(f)
    while 1:
        if l.upper().startswith('*SET_'+stype.upper()):
            while 1:
                try:
                    nset=int(l[:10])
                    break
                except:
                    pass
                l=next(f)
            if num:
                if not nset in num:
                    l=next(f)
                    continue
            rez[nset]=[]
            l=next(f)
            while not l.startswith('*'):
                try:
                    rez[nset]+=[int(ll) for ll in l.split()]
                except:
                    pass
                l=next(f)
            while rez[nset][-1]==0:
                rez[nset].pop()
            continue
        try:
            l=next(f)
        except:
            break
    f.close()
    return rez

class boundBox(object):
    def __init__(self, xmin=0, xmax=1, ymin=0, ymax=1, zmin=0, zmax=1):
        self.minX=xmin
        self.maxX=xmax
        self.minY=ymin
        self.maxY=ymax
        self.minZ=zmin
        self.maxZ=zmax
    def __init__(self, d):
        self.minX=d['minX']
        self.maxX=d['maxX']
        self.minY=d['minY']
        self.maxY=d['maxY']
        self.minZ=d['minZ']
        self.maxZ=d['maxZ']
    def isPointInside(self, x, y, z):
        rez=x>=self.minX and x<=self.maxX
        rez=rez and y>=self.minY and y<=self.maxY
        rez=rez and z>=self.minZ and z<=self.maxZ
        return rez
    def offsetBounds(self, dist):
        self.minX-=dist
        self.maxX+=dist
        self.minY-=dist
        self.maxY+=dist
        self.minZ-=dist
        self.maxZ+=dist
    def __repr__(self):
        s="[{}, {}\n".format(self.minX, self.maxX)
        s+=" {}, {}\n".format(self.minY, self.maxY)
        s+=" {}, {}]\n".format(self.minZ, self.maxZ)
        return s

class node(object):
    __slots__ = ['n', 'x', 'y', 'z']
    def __init__(self, n=1, crds=[0.0, 0.0, 0.0]):
        self.n=n
        self.x=crds[0]
        self.y=crds[1]
        self.z=crds[2]

    def get_crds(self):
        return [self.x, self.y, self.z]

    def set_crds(self, data):
        self.x=data[0]
        self.y=data[1]
        self.z=data[2]
    crds=property(get_crds, set_crds)

    def distTo(self, x, y, z):
        return sqrt((x-self.x)**2+(y-self.y)**2+(z-self.z)**2)

    def isInsideBox(self, box):
        return box.isPointInside(self.x, self.y, self.z)

class element(object):
    __slots__ = ['n', '_nodes', 'part', 'etype']
    def __init__(self, n=1, nodes=[], part=1, etype='*ELEMENT'):
        self.n=n
        while 0 in nodes:
            nodes.remove(0)
        self._nodes=nodes
        self.part=part
        self.etype=etype

    def getNodesCount(self):
        return len(set(self.uniqNodes))
    nodesCount = property(getNodesCount)

    def getUniqNodes(self):
        nodes=self._nodes
        for nn in nodes:
            while nodes.count(nn)-1:
                nodes.remove(nn)
        return nodes
    uniqNodes=property(getUniqNodes)

    def getFaces(self):
        nodes=self.uniqNodes
        if self.etype=="*ELEMENT_SOLID":
            if self.nodesCount==8:
                finds=[[1,4,3,2],
                       [5,6,7,8],
                       [2,3,7,6],
                       [1,5,8,4],
                       [3,4,8,7],
                       [1,2,6,5]]
            if self.nodesCount==6:
                finds=[[1,4,3,2],
                       [1,5,6,4],
                       [2,3,6,5],
                       [1,2,5],
                       [4,6,3]]
            if self.nodesCount==4:
                finds=[[3,2,1],
                       [2,4,1],
                       [3,4,2],
                       [1,4,3]]
            if self.nodesCount==10:
                finds=[[3,2,1,6,5,7],
                       [2,4,1,8,10,5],
                       [3,4,2,9,8,6],
                       [1,4,3,10,9,7]]

        if self.etype=="*ELEMENT_SHELL":
            if self.nodesCount==4:
                finds=[[1,2],
                       [2,3],
                       [3,4],
                       [4,1]]
            if self.nodesCount==3:
                finds=[[1,2],
                       [2,3],
                       [3,1]]
            if self.nodesCount==6:
                finds=[[1,2,4],
                       [2,3,5],
                       [3,1,6]]
        rez=[]
        for idx in finds:
            rez.append([nodes[i-1] for i in idx])
        return rez

    def getEdges(self):
        nodes=self.uniqNodes
        if self.etype=="*ELEMENT_SOLID":
            if self.nodesCount==8:
                finds=[[1,2],
                       [2,3],
                       [3,4],
                       [4,1],
                       [5,6],
                       [6,7],
                       [7,8],
                       [8,5],
                       [1,5],
                       [2,6],
                       [3,7],
                       [4,8]]
            if self.nodesCount==6:
                finds=[[1,2],
                       [2,3],
                       [3,4],
                       [4,1],
                       [5,6],
                       [1,5],
                       [2,5],
                       [3,6],
                       [4,6]]
            if self.nodesCount==4 or self.nodesCount==10:
                finds=[[1,2],
                       [2,3],
                       [3,1],
                       [1,4],
                       [2,4],
                       [3,4]]
        if self.etype=="*ELEMENT_SHELL":
            if self.nodesCount==4:
                finds=[[1,2],
                       [2,3],
                       [3,4],
                       [4,1]]
            if self.nodesCount==3 or self.nodesCount==6:
                finds=[[1,2],
                       [2,3],
                       [3,1]]
        rez=[]
        for idx in finds:
            rez.append([nodes[i-1] for i in idx])
        return rez

    def getEdgesAngles(self):
        nodes=self.uniqNodes
        if self.etype=="*ELEMENT_SOLID":
            if self.nodesCount==8:
                finds=[ [1,9], [1,-4],  [9,-4],
                        [-1,2], [2,10], [-1,10],
                       [-2,3], [-2,11], [3,11],
                       [-3,12], [4,12], [-3,4],
                       [-9,5], [5,-8], [-8,-9],
                       [-5,6], [6,-10], [-5,-10],
                       [-6,7], [-6,11], [7,-11],
                       [-7,8], [8,-12], [-7,-12]]
            if self.nodesCount==6:
                finds=[ [1,6], [1,-4],  [6,-4],
                        [-1,2], [2,7], [-1,7],
                       [-2,3], [-2,8], [3,8],
                       [-3,9], [4,9], [-3,4],
                       [-6,-7], [5,-7], [5,-6],
                       [-9,-8], [-5,-8], [-5,-9]]
            if self.nodesCount==4 or self.nodesCount==10:
                finds=[ [1,-3], [1,4],  [4,-3],
                        [-1,2], [5,2], [-1,5],
                       [-2,3], [-2,6], [3,6],
                       [-4,-5], [-5,-6], [-6,-4]]
        if self.etype=="*ELEMENT_SHELL":
            if self.nodesCount==4:
                finds=[[1,-4],
                       [-1,2],
                       [-2,3],
                       [4,-3]]
            if self.nodesCount==3 or self.nodesCount==6:
                finds=[[1,-3],
                       [2,-1],
                       [3,-2]]
        return finds

    def _getNodes(self):
        rez=self._nodes
        if self.etype=='*ELEMENT_SOLID':
            if len(rez)==4:
                rez+=[rez[-1]]*4
            if len(rez)==6:
                rez.insert(-1, rez[-2])
                rez.append(rez[-1])
        if self.etype=='*ELEMENT_SHELL':
            if len(rez)==3:
                rez+=[rez[-1]]
        return rez
    def _setNodes(self, nodes):
        self._nodes=nodes

    nodes=property(_getNodes, _setNodes)

    def containsNds(self, nds=1):
        nodes=self.uniqNodes
        if type(nds)==int:
            return nds in nodes
        rez=True
        for n in nds:
            rez=rez and (n in nodes)
        return rez

    def getFaceNum(self, face):
        if not self.containsNds(face):
            return 0
        allFaces=self.getFaces()
        face=set(face)
        rez=0
        for i, f in enumerate(allFaces):
            if face==set(f):
                rez=i+1
                break
        return rez

class lsdyna_model:
    __slots__ = ['kwds', 'nodes', 'solids', 'shells', 'progress']

    def __init__(self, fname = None, procIncludes = False, progress = False):
        """
        :param fname: имя k-файла с задачей LS-DYNA
        :param procIncludes: Если True, то обрабатываются файлы из карты *INCLUDE
        :param progress: Показывать прогресс выполнения длительных задач
        """
        self.kwds = defaultdict(list)
        self.nodes = None
        self.solids = None
        self.shells = None
        self.progress = progress
        if fname:
            self.read_data(fname, procIncludes)
            self.proc_mesh()
            print(self)

    def read_data(self, fname, includes = False):
        """
        Метод считывает информацию из k-файла с учетом вложений.

        :param fname: Имя k-файла
        :param includes: Если True, то обрабатываются файлы из карты *INCLUDE
        :return:
        """
        self.read_one_file(fname)
        if includes:
            while self.kwds['*INCLUDE']:
                ff = self.kwds['*INCLUDE'].pop()
                for _ in ff.split('\n')[:-1]:
                    if _[0] == '$':
                        continue
                    self.read_one_file(_)
            self.kwds.pop('*INCLUDE')

    def read_one_file(self, fname):
        """
        Метод чтения информации из отдельного файла.

        :param fname: имя файла
        :return:
        """
        print('Reading data from {}...'.format(fname), end='\t', flush=True)
        with open(fname, 'r') as k:
            block = []
            for l in k:
                if l.startswith('*'):
                    self.proc_block(block)
                    block = [l]
                else:
                    block.append(l)
            self.proc_block(block)
        print('Done...')

    def proc_block(self, block):
        """
        Метод разбирает блок ключевого слова и сохраняет
        содержание карточки в словарь, ключом является
        keyword.

        :param block: Массив строк.
        :return:
        """
        if block:
            if block[0][0] == '$':
                pass
            elif block[0].upper().strip() in ['*KEYWORD', '*END']:
                pass
            else:
                self.kwds[block[0].upper().strip()].append(''.join(block[1:]))

    def proc_nodes(self):
        """
        Разбирает карточки *NODE и создает слооварь узлов
        self.nodes={номер_узла: узел, ...}

        :return:
        """
        print('Processing nodes...')
        lines = ''.join(self.kwds['*NODE']).split('\n')
        self.nodes = {}
        cprogress = ConsoleProgress(len(lines)) if self.progress else None
        for i, l in enumerate(lines):
            if not l:
                continue
            if l[0] == '$':
                continue
            if "," in l:
                l += "," * (3 - l.count(","))
                ll = l.split(",")
            else:
                l += " " * (8 + 3 * 16 - len(l))
                ll = [l[0:8], l[8:24], l[24:40], l[40:56]]
            nd = node(int(ll[0]), [float(ll[ii]) for ii in range(1, 4)])
            self.nodes[nd.n] = nd
            if cprogress: cprogress.redraw()
        if cprogress: cprogress.finalize()
        del self.kwds['*NODE']

    def proc_elements(self):
        """
        Разбирает карточки *ELEMENT_SOLID и *ELEMENT_SHELL и создает
        словари вида:
        self.solids (self.shells) = {номер_части: {номер_элемента: элемент, ...}, ...}

        :return:
        """
        self.solids = defaultdict(dict)
        self.shells = defaultdict(dict)
        elements = {'solids': ['*ELEMENT_SOLID', self.solids],
                    'shells': ['*ELEMENT_SHELL', self.shells]
                    }
        for etype in elements:
            lines = ''.join(self.kwds[elements[etype][0]]).split('\n')
            print('Processing {}...'.format(etype))
            cprogress = ConsoleProgress(len(lines)) if self.progress else None
            for i, l in enumerate(lines):
                if not l:
                    continue
                if l[0] == '$':
                    continue
                if "," in l:
                    ll = l.split(",")
                else:
                    ll = l.split()
                tmp = list(map(int, ll))
                el = element(tmp[0], tmp[2:], tmp[1], etype)
                elements[etype][1][tmp[1]][el.n] = el
                if cprogress: cprogress.redraw()
            if cprogress: cprogress.finalize()
            del self.kwds[elements[etype][0]]

    def proc_mesh(self):
        """
        Разбор карточек с узлами и элементами.

        :return:
        """
        self.proc_nodes()
        self.proc_elements()

    def __repr__(self):
        rez='Number of keywords: {}'.format(len(self.kwds))
        if self.nodes:
            rez+='\nNumber of nodes: {}'.format(len(self.nodes))
        if self.solids:
            rez += '\nNumber of solid parts: {}'.format(len(self.solids))
            for p, s in self.solids.items():
                rez+='\n\tPart {} consists of {} elements'.format(p, len(s))
        if self.shells:
            rez += '\nNumber of shell parts: {}'.format(len(self.shells))
            for p, s in self.shells.items():
                rez+='\n\tPart {} consists of {} elements'.format(p, len(s))
        return rez

    def save_nodes(self, fname, nShift, mode='a', end=False):
        print('Saving nodes...')
        if not self.nodes:
            return
        if type(fname)==str:
            f=open(fname, mode)
        else:
            f=fname
        if mode=='w':
            f.write("*KEYWORD\n")
        f.write('*NODE\n')
        for n, c in self.nodes.items():
            f.write("%8d%16.9e%16.9e%16.9e\n" % tuple([n+nShift]+c.get_crds()))
        if end:
            f.write("*END")
        if type(fname) == str: f.close()

    def save_elements(self, fname, nShift=0, eShift=0, pShift=0, mode='a', end=False):
        print('Saving elements...')
        if type(fname)==str:
            f=open(fname, mode)
        else:
            f=fname
        if mode=='w':
            f.write("*KEYWORD\n")
        if self.solids:
            f.write('*ELEMENT_SOLID\n')
            for p, ps in self.solids.items():
                for e in ps.values():
                    nds=e.nodes
                    fstr='%8d'*((len(nds)+2))+'\n'
                    f.write(fstr % tuple([e.n+eShift, e.part+pShift]+[nd+nShift for nd in nds]))
        if self.shells:
            f.write('*ELEMENT_SHELL\n')
            for p, ps in self.shells.items():
                for e in ps.values():
                    nds = e.nodes
                    fstr='%8d'*((len(nds)+2))+'\n'
                    f.write(fstr % tuple([e.n+eShift, e.part+pShift]+[nd+nShift for nd in nds]))
        if end:
            f.write("*END")
        if type(fname) == str: f.close()

    def save_model(self, fname, nShift=0, eShift=0, pShift=0):
        with open(fname, 'w') as f:
            f.write('*KEYWORD\n')
            for kw in self.kwds:
                for data in self.kwds[kw]:
                    f.write('{}\n'.format(kw))
                    f.write(data)
            self.save_nodes(f, nShift=nShift)
            self.save_elements(f, nShift=nShift, eShift=eShift, pShift=pShift)
            f.write('*END\n')

class model:
    def __init__(self, fname=None, includes=False):
        self.cards=[]
        self.prts=set([])
        self.nds={}
        self.els={}
        self.dens={}
        self.data=[]
        self.keywords=[]
        self.indexes=[]
        self.tol=1e-10
        self.segmentSets={}
        if fname!=None:
            self.getDataFromFile(fname, includes)
            self.read_elements()
            self.read_nodes()

    def findDuplicateNodes(self, tol=1e-6):
        nds = self.getExternalNodes()
        pairs = defaultdict(list)
        checked = []
        for n1,n2 in combinations(nds, 2):
            if (n1 in checked) and (n2 in checked): continue
            if self.getNdDistance(n1,n2)<=tol:
                pairs[n1].append(n2)
                checked.append(n1)
                checked.append(n2)
        print(len(pairs), 'duplicated nodes have been found...')
        return pairs

    def sewDuplicateNodes(self, tol=1e-6):
        dup_nodes = self.findDuplicateNodes(tol)
        for n1, ns in dup_nodes.items():
            for n2 in ns:
                if not n2 in self.nds: continue
                self.nds.pop(n2)
                for e in self.els.values():
                    while n2 in e.nodes:
                        e.nodes[e.nodes.index(n2)]=n1
                if hasattr(self, 'allLines'):
                    for segments in self.allLines.values():
                        for segment in segments:
                            while n2 in segment:
                                segment[segment.index(n2)]=n1
                if hasattr(self, 'allPoints'):
                    for p in self.allPoints:
                        if self.allPoints[p]==n2:
                            self.allPoints[p]=n1

    def getDataFromFile(self, fname, includes):
        self.fname=fname
        self.data=readData(fname, includes)
        self.keywords, self.indexes=searchKeywords(self.data)

    def read_nodes(self):
        print('Reading nodes...')
        if self.keywords.count('*NODE'):
            for i, kw in enumerate(self.keywords):
                if kw=='*NODE':
                    for l in self.data[self.indexes[i]+1: self.indexes[i+1]]:
                        if l.count(","): l+=","*(3-l.count(",")); ll=l.split(",")
                        else:
                            l+=" "*(8+3*16-len(l))
                            ll=[l[0:8],l[8:24],l[24:40],l[40:56]]
                        nd=node(int(ll[0]), [float(ll[ii]) for ii in range(1,4)])
                        self.nds[nd.n]=nd

    def read_elements(self, pts=None):
        print('Reading elements...')
        for i, kw in enumerate(self.keywords):
            if kw[:8]=="*ELEMENT":
                for l in self.data[self.indexes[i]+1: self.indexes[i+1]]:
                    if l.count(","): ll=l.split(",")
                    else: ll=l.split()
                    tmp=[]
                    for i in ll: tmp.append(int(i))
                    if pts!=None and not (tmp[1] in pts): continue
                    el=element(tmp[0], tmp[2:], tmp[1], kw)
                    self.prts.add(tmp[1])
                    self.els[el.n]=el

    def save_nodes(self, fname, nShift, mode='a', end=False):
        print('Saving nodes...')
        if not self.nds:
            return
        f=open(fname, mode)
        if mode=='w':
            f.write("*KEYWORD\n")
        f.write('*NODE\n')
        for n, c in self.nds.items():
            f.write("%8d%16.9e%16.9e%16.9e\n" % tuple([n+nShift]+c.get_crds()))
        if end:
            f.write("*END")
        f.close()

    def save_elements(self, fname, nShift=0, eShift=0, pShift=0, mode='a', end=False):
        print('Saving elements...')
        f=open(fname, mode)
        if mode=='w':
            f.write("*KEYWORD\n")
        solids, shells = self.getByType()
        if solids:
            f.write('*ELEMENT_SOLID\n')
            for n in solids:
                e=self.els[n]
                nds=e.nodes
                fstr='%8d'*((len(nds)+2))+'\n'
                f.write(fstr % tuple([n+eShift, e.part+pShift]+[nd+nShift for nd in e.nodes]))
        if shells:
            f.write('*ELEMENT_SHELL\n')
            for n in shells:
                e=self.els[n]
                nds=e.nodes
                fstr='%8d'*((len(nds)+2))+'\n'
                f.write(fstr % tuple([n+eShift, e.part+pShift]+[nd+nShift for nd in e.nodes]))
        if end:
            f.write("*END")
        f.close()

    def save_model(self, fname, nShift=0, eShift=0, pShift=0, end=True, mode='w'):
        f=open(fname, mode)
        if mode=='w': f.write('*KEYWORD\n')
        for i, kw in enumerate(self.keywords):
            if kw in ['*INCLUDE', '*END', '*EN', '*NODE', '*ELEMENT_SHELL',
                      '*ELEMENT_SOLID', '*KEYWORD']: continue
            f.writelines(self.data[self.indexes[i]:self.indexes[i+1]])
        f.close()
        self.save_nodes(fname, nShift)
        self.save_elements(fname, nShift, eShift, pShift)
        if end:
            f=open(fname, 'a')
            f.write('*END')
        f.close()

    def saveSegmentSet(self, n, fname, nShift=0, mode='w', end='True'):
        if n not in self.segmentSets: return
        if not self.segmentSets[n]: return
        fout=open(fname, mode)
        if mode=='w':
            fout.write("*KEYWORD\n")
        fout.write("*SET_SEGMENT_LIST\n")
        fout.write("{0:10d}\n".format(n))
        for s in self.segmentSets[n]:
            if len(s)==3:
                s.append(s[-1])
            fs='%10d'*len(s)+'\n'
            fout.write(fs % tuple([ss+nShift for ss in s]))
        if end:
            fout.write('*END')
        fout.close()

    def get_element_center(self, en):
#        if not en in self.els.keys():
#            return
        try:
            el=self.els[en] if isinstance(en, int) else en
        except:
            return 0
        #print(el)
        nds=el._nodes
        mean_crds=[0,0,0]
        for n in nds:
            crds=self.nds[n].crds
            for i in range(3):
                mean_crds[i]+=crds[i]
        nn=len(nds)
        return [mm/nn for mm in mean_crds]

    def getShellNormal(self, en):
        e=self.els[en]
        n1=Point(*self.nds[e.nodes[0]].crds)
        n2=Point(*self.nds[e.nodes[1]].crds)
        n3=Point(*self.nds[e.nodes[2]].crds)
        return ((n2-n1).cross(n3-n1)).normalized()

    def getElementCenterByNodes(self, en, snds=[]):
        if not en in list(self.els.keys()):
            return None
        nds=set(self.els[en].nodes).intersection(set(snds))
        if not nds:
            return None
        mean_crds=[0,0,0]
        for n in nds:
            crds=self.nds[n].crds
            for i in range(3):
                mean_crds[i]+=crds[i]
        nn=len(nds)
        return [mm/nn for mm in mean_crds]

    def del_element(self, en):
    #       if not en in self.els.keys():
    #           return
        self.els.pop(en)

    def del_elements(self, ens):
        for en in ens:
#            if not en in self.els.keys():
#                continue
            self.els.pop(en)

    def del_node(self, nn):
#        if not nn in self.nds.keys():
#            return
        self.nds.pop(nn)

    def del_nodes(self, nns):
        for nn in nns:
    #           if not nn in self.nds.keys():
    #               continue
            self.nds.pop(nn)

    def get_elements_nodes(self, ens=None):
        if not ens:
            ens=list(self.els.keys())
        rez=[]
        for en in ens:
            rez+=self.els[en].nodes
        return set(rez)

    def get_unref_nodes(self):
        allnds=self.get_elements_nodes()
        rez=[]
        for nn in list(self.nds.keys()):
            if not nn in allnds:
                rez.append(nn)
        return rez

    def getElementsByParts(self, parts=None, values=False):
        rezKeys=[]
        rezValues=[]
        if type(parts)==int:
            parts=[parts]

        for e in list(self.els.values()):
            if e.part in parts:
                rezKeys.append(e.n)
                rezValues.append(e)
        return rezValues if values else rezKeys

    def get_parts_mats(self):
        parts={}
        mats={}
        for i, kw in enumerate(self.keywords):
            if kw=='*PART':
                l=self.data[self.indexes[i]+2]
                parts[int(l[0:10])]=int(l[20:30])
            if kw.startswith('*MAT'):
                di=1 if not 'TITLE' in kw else 2
                l=self.data[self.indexes[i]+di]
                mats[int(l[0:10])]=float(l[10:20])
        return parts, mats

    def get_densities(self):
        print('Reading material densities...')
        parts, mats = self.get_parts_mats()
        for p, m in parts.items():
            self.dens[p]=mats[m]

    def getByType(self):
        solids=[]
        shells=[]
        for e in list(self.els.values()):
            if e.etype=='*ELEMENT_SOLID':
                solids.append(e.n)
            if e.etype=='*ELEMENT_SHELL':
                shells.append(e.n)
        return solids, shells
    def getElementArea(self, en):
        nc=self.els[en].nodesCount
        if self.els[en].etype!='*ELEMENT_SHELL':
            return 0
        p=[]
        for i in range(nc):
            pp=Point(*self.nds[self.els[en].nodes[i]].crds)
            p.append(pp)
        S=abs(0.5*(p[0]-p[1]).cross(p[2]-p[1]).norm())
        if nc==4:
            S+=abs(0.5*(p[0]-p[3]).cross(p[2]-p[3]).norm())
        return S

    def getElementVolume(self, en):
        V=0
        nc=self.els[en].nodesCount
        if self.els[en].etype=='*ELEMENT_SOLID':
            if nc==4:
                idxs=[0,1,2,3]
                points=[self.nds[self.els[en].nodes[i]].crds for i in idxs]
                V+=calcV4points(points)

            if nc==6:
                idxs=[0,1,3,4]
                points=[self.nds[self.els[en].nodes[i]].crds for i in idxs]
                V+=calcV4points(points)
                idxs=[4,1,3,2]
                points=[self.nds[self.els[en].nodes[i]].crds for i in idxs]
                V+=calcV4points(points)
                idxs=[4,3,2,7]
                points=[self.nds[self.els[en].nodes[i]].crds for i in idxs]
                V+=calcV4points(points)

            if nc==8:
                idxs=[0,4,5,7]
                points=[self.nds[self.els[en].nodes[i]].crds for i in idxs]
                V+=calcV4points(points)
                idxs=[2,6,5,7]
                points=[self.nds[self.els[en].nodes[i]].crds for i in idxs]
                V+=calcV4points(points)
                idxs=[0,2,3,7]
                points=[self.nds[self.els[en].nodes[i]].crds for i in idxs]
                V+=calcV4points(points)
                idxs=[0,2,5,7]
                points=[self.nds[self.els[en].nodes[i]].crds for i in idxs]
                V+=calcV4points(points)
                idxs=[0,1,2,5]
                points=[self.nds[self.els[en].nodes[i]].crds for i in idxs]
                V+=calcV4points(points)
        if self.els[en].etype=='*ELEMENT_SHELL':
            p=[]
            for i in range(nc):
                pp=Point(*self.nds[self.els[en].nodes[i]].crds)
                p.append(pp)
            S=abs(0.5*(p[0]-p[1]).cross(p[2]-p[1]).norm())
            if nc==4:
                S+=abs(0.5*(p[0]-p[3]).cross(p[2]-p[3]).norm())
            c=self.get_element_center(en)
            V=2*pi*c[0]*S
        return V

    def get_items_on_axis(self):
        nox=[]
        for num, nd in self.nds.items():
            if abs(nd.crds[0])<self.tol:
                nox.append(num)
                self.nds[num].x=0
        eox=[[],[]]
        noxs=set(nox)
        for ne, el in self.els.items():
            nns=set(el.uniqNodes)
            nint=noxs.intersection(nns)
            if len(nint)==1:
                eox[0].append(ne)
            if len(nint)==2:
                eox[1].append(ne)
        return nox, eox

    def twist(self, y0, dl, y1=1000., direction=1.):
        print('Twisting mesh...')
        for nn, nd in self.nds.items():
            if nd.y<=y0:
                continue
            if nd.y<=y1:
                ang=direction*2.0*math.pi/dl*(nd.y-y0)
            else:
                ang=direction*2.0*math.pi/dl*(y1-y0)
            ca=math.cos(ang)
            sa=math.sin(ang)
            x=nd.x
            z=nd.z
            nd.x=x*ca+z*sa
            nd.z=-x*sa+z*ca

    def __add__(self, mesh2):
        newMesh=model()
        newMesh.nds=self.nds.copy()
        newMesh.els=self.els.copy()
        nshift=0
        eshift=0
        newMesh.prts=self.prts.union(mesh2.prts)
        try:
            nshift=max(self.nds.keys())
        except:
            pass
        try:
            eshift=max(self.els.keys())
        except:
            pass
        for nn, nd in mesh2.nds.items():
            newnode=node(nn+nshift, nd.crds)
            newMesh.nds[nn+nshift]=newnode
        for en, el in mesh2.els.items():
            newnodes=[nn+nshift for nn in el.nodes]
            newelement=element(en+eshift, nodes=newnodes, part=el.part, etype=el.etype)
            newMesh.els[en+eshift]=newelement
        return newMesh

    def __mul__(self, mesh2):
        newMesh=model()
        newMesh.nds=self.nds.copy()
        newMesh.els=self.els.copy()
        nshift=0
        eshift=0
        newMesh.prts=self.prts.union(mesh2.prts)
        for nn, nd in mesh2.nds.items():
            newnode=node(nn, nd.crds)
            newMesh.nds[nn]=newnode
        for en, el in mesh2.els.items():
            newnodes=[nn for nn in el.nodes]
            newelement=element(en, nodes=newnodes, part=el.part, etype=el.etype)
            newMesh.els[en]=newelement
        return newMesh

    def fix2Dorientation(self):
        n=0
        for en, el in self.els.items():
            ns=el.uniqNodes
            if ((self.nds[ns[1]].x-self.nds[ns[0]].x)*
                (self.nds[ns[2]].y-self.nds[ns[0]].y)-
                (self.nds[ns[2]].x-self.nds[ns[0]].x)*
                (self.nds[ns[1]].y-self.nds[ns[0]].y))<0:
                if len(ns)<=4:
                    ns.reverse()
                if len(ns)==6:
                    ns1=ns[0:3]
                    ns2=ns[3:5]
                    ns1.reverse()
                    ns2.reverse()
                    ns=ns1+ns2+[ns[5]]
                self.els[en].nodes=ns
                n+=1
        print('Repaired %d element directions' % (n,))

    def swap2Dorientation(self):
        for el in list(self.els.values()):
            ns=el.nodes
            ns.reverse()
            el.nodes=ns

    def calcMassProps(self, parts=[]):
        '''
        Return mass, volume, center of mass and inertia tensor.
        '''
        if not self.dens:
            self.get_densities()
        if not parts:
            els = self.els.keys()
        else:
            els = self.getElementsByParts(parts)
        print('Calculating mass properties...')
        m=0
        I=np.array([[0,0,0],[0,0,0],[0,0,0]], dtype=np.float64)
        v=0
        cm=np.array([0,0,0], dtype=np.float64)
        count=len(els)
        proceded=0
        progress0=0
        idx=[[1,2],[0,2],[0,1]]
        print('Calculating mass, volume and center of mass...')
        emss={}
        ecntrs={}
        pr = ConsoleProgress(count)
        for en in els:
            e = self.els[en]
            proceded+=1.
#            progress=int(round(proceded/count*100.))
#            changed=progress0!=progress
#            progress0=progress
#            if not progress%10 and changed:
            pr.current = proceded
            pr.redraw()
                #print(progress, '% done...')
#                pr.chekRedraw()
            # if not e.part in parts and len(parts)>0:
            #     continue
            ev=self.getElementVolume(en)
            ec=np.array(self.get_element_center(en))
            em=ev*self.dens[e.part]
            cm=cm+em*ec
            m+=em
            v+=ev
            emss[en]=em
            ecntrs[en]=ec
        pr.finalize()
        cm=cm/m
        print('Calculating inertia tensor...')
        for en in list(emss.keys()):
            ec=ecntrs[en]-cm
            em=emss[en]
            for i in range(3):
                I[i][i]=I[i][i]+(ec[idx[i][0]]**2+ec[idx[i][1]]**2)*em
                for j in range(3):
                    if j>=i:
                        continue
                    I[i][j]=I[i][j]-ec[i]*ec[j]*em
                    I[j][i]=I[i][j]
        return m, v, cm, I

    def readSegmentSets(self):
        self.segmentSets={}
        if self.keywords.count('*SET_SEGMENT_LIST') or self.keywords.count('*SET_SEGMENT_LIST_TITLE') or self.keywords.count('*SET_SEGMENT'):
            for i, kw in enumerate(self.keywords):
                if '*SET_SEGMENT' in kw:
                    nn=2 if 'TITLE' in kw else 1
                    lines=self.data[self.indexes[i]+nn:self.indexes[i+1]]
                    setNum=int(lines[0].split()[0])
                    segments=[]
                    for l in lines[1:]:
                        n=min(len(l)//10,4)
                        segment=[]
                        for i in range(n):
                            segment.append(int(l[i*10:(i+1)*10]))
                        segments.append(segment)
                    self.segmentSets[setNum]=segments

    def reflect(self, norm=0):
        nshift=max(self.nds.keys())
        eshift=max(self.els.keys())
        for n in list(self.nds.values()):
            if abs(n.crds[norm])<=self.tol:
                continue
            crds=n.crds
            crds[norm]*=-1
            newn=node(n.n+nshift, crds)
            self.nds[newn.n]=newn
        for e in list(self.els.values()):
            nns=list(e.uniqNodes)
            for i, n in enumerate(nns):
                if not abs(self.nds[n].crds[norm])<=self.tol:
                    nns[i]+=nshift
            if e.etype=='*ELEMENT_SHELL':
                nns.reverse()
            if e.etype=='*ELEMENT_SOLID':
                if len(nns)==4:
                    tmp=nns[0:3]
                    tmp.reverse()
                    nns=tmp+[nns[-1]]
                if len(nns)==8:
                    tmp1=nns[:4]
                    tmp2=nns[4:]
                    tmp1.reverse()
                    tmp2.reverse()
                    nns=tmp1+tmp2
            newe=element(e.n+eshift, nodes=nns, part=e.part, etype=e.etype)
            self.els[newe.n]=newe

    def scale(self, sx=1., sy=1., sz=1.):
        for n in list(self.nds.keys()):
            self.nds[n].x*=sx
            self.nds[n].y*=sy
            self.nds[n].z*=sz

    def translate(self, dx=0., dy=0., dz=0.):
        for n in list(self.nds.keys()):
            self.nds[n].x+=dx
            self.nds[n].y+=dy
            self.nds[n].z+=dz

    def rotate(self, angle=0., axis=[0,0,1]):
        q=Quaternion(angle*pi/180., Point(*axis))
        for n, nd in self.nds.items():
            self.nds[n].crds=q.rotatePoint(Point(*nd.crds)).data

    def fixSegmentsOrientations(self, segments, prts=[]):
        print('Fixing segments orientations...')
        n=0
        if not prts:
            prts=self.prts
        segmentNodes=set(chain(*segments))
        elements=self.getElementsByParts(prts)
        allSegments=[]
        for en in elements:
            if not segmentNodes.intersection(set(self.els[en]._nodes)):
                continue
            allSegments.extend(self.els[en].getFaces())
        segmentsSorted=deepcopy(segments)
        allSegmentsSorted=deepcopy(allSegments)
        for i, s in enumerate(allSegmentsSorted):
            allSegmentsSorted[i].sort()
        for i, s  in enumerate(segmentsSorted):
            segmentsSorted[i].sort()
        rez=[]
        for i, s in enumerate(segmentsSorted):
            if s in allSegmentsSorted:
                rez.append(allSegments[allSegmentsSorted.index(s)])
                if shift(segments[i], segments[i].index(rez[-1][0]))!=rez[-1]:
                    n+=1
            else:
                print(s)
        print("Fixed ", n, " segments orientations...")
        return rez

    def nsetSortByY(self, ns):
        N=len(ns)
        for i in range(N-1):
            for j in range(i+1,N):
                if self.nds[ns[j]].y>self.nds[ns[i]].y:
                    ns[i], ns[j] = ns[j], ns[i]

    def esetSortByCentrY(self, eset):
        cy=[]
        for en in eset:
            cy.append(self.get_element_center(en)[1])
        N=len(eset)
        for i in range(N-1):
            for j in range(i+1,N):
                if cy[j]>cy[i]:
                    cy[i], cy[j] = cy[j], cy[i]
                    eset[i], eset[j] = eset[j], eset[i]

    def esetSortByY(self, eset, nset):
        cy=[]
        for en in eset:
            cy.append(self.getElementCenterByNodes(en, nset)[1])
        N=len(eset)
        for i in range(N-1):
            for j in range(i+1,N):
                if cy[j]>cy[i]:
                    cy[i], cy[j] = cy[j], cy[i]
                    eset[i], eset[j] = eset[j], eset[i]

    def _getBoundBox(self, prts=[]):
        if not prts:
            nodes=list(self.nds.keys())
        else:
            els=self.getElementsByParts(prts)
            nodes=self.get_elements_nodes(els)
        maxX=-1e6
        maxY=-1e6
        maxZ=-1e6
        minX=1e6
        minY=1e6
        minZ=1e6
        for nn in nodes:
            n=self.nds[nn]
            maxX=max(maxX, n.x)
            maxY=max(maxY, n.y)
            maxZ=max(maxZ, n.z)
            minX=min(minX, n.x)
            minY=min(minY, n.y)
            minZ=min(minZ, n.z)
        return {'minX': minX, 'maxX': maxX,
                'minY': minY, 'maxY': maxY,
                'minZ': minZ, 'maxZ': maxZ
                }

    boundBox=property(_getBoundBox)

    def getNodesInsideBox(self, box, prts=[]):
        if not prts:
            nodes=list(self.nds.keys())
        else:
            els=self.getElementsByParts(prts)
            nodes=self.get_elements_nodes(els)
        rez=[]
        for nn in nodes:
            if self.nds[nn].isInsideBox(box):
                rez.append(nn)
        return rez

    def saveCCXmesh(self, fname, solidType='C3D8', shellType='CAX4', plane=False,
                    poffset=0, noffset=0, eoffset=0, mode='w', nset='allnodes'):
        fout=open(fname, mode)
        fout.write('*Node, nset={0}\n'.format(nset))
        d=2 if plane else 3
        fs='%8d,'+ ','.join(d*[' %16.8e'])+'\n'
        for nd in list(self.nds.values()):
            fout.write(fs % tuple([nd.n+noffset]+nd.crds[:d]))
        for p in self.prts:
            els=self.getElementsByParts([p])
            tp=solidType if self.els[els[0]].etype=='*ELEMENT_SOLID' else shellType
            fout.write('*Element, type=%s, elset=Part%d\n' % (tp, p+poffset))
            for e in els:
                el=self.els[e]
                nc=el.getNodesCount()
                fs=', '.join(['%8d']*(nc+1))+'\n'
                if nc>8:
                    fs=fs[:30]+'\n'+fs[30:]
                fout.write(fs % tuple([el.n+eoffset]+[nn+noffset for nn in el.nodes[:nc]]))
        fout.close()

    def getFacesNumbers(self, faces, elements='all'):
        if elements=='all':
            elements=list(self.els.keys())
        if type(faces[0])!=list:
            faces=[faces]
        rez=[]
        n1=len(faces)
        n2=len(self.els)
        total=n1*n2
        for i, f in enumerate(faces):
            for j, en in  enumerate(elements):
                e=self.els[en]
                n=e.getFaceNum(f)
                if n:
                    rez.append([e.n,n])
            print(progress(n2*(i+1), total), end=' ')
        print()
        return rez

#    def getFacesNumbersParallel(self, faces, elements='all'):
#        if elements=='all':
#            elements=self.els.keys()
#        if type(faces[0])!=list:
#            faces=[faces]
#        serv=pp.Server()
#        ncpus=max(serv.get_ncpus()-1,1)
#        print('Serching parallel.', ncpus, 'workers...')
#        NN=len(elements)//ncpus+1
#        ee=splitByN(elements, NN)
#        jobs=[serv.submit(self.getFacesNumbers, (faces, eee), (), ('lsmesh_lib',)) for eee in ee]
#        serv.wait()
#        serv.print_stats()
#        rez=[]
#        for j in jobs:
#            rez+=j()
#        return rez

    def getFacesDict(self, faces):
        rez={}
        for i in range(6):
            rez[i+1]=[]
        if type(faces[0])!=list:
            faces=[faces]
        n1=len(faces)
        n2=len(self.els)
        total=n1*n2
        for i, f in enumerate(faces):
            for j, e in  enumerate(self.els.values()):
                n=e.getFaceNum(f)
                if n:
                    rez[n].append(e.n)
            print(progress(n2*(i+1), total), end=' ')
        print()
        return rez

    def cleanDegenerated(self):
        n=0
        deletedNodes=[]
        for e in list(self.els.values()):
            v=self.getElementVolume(e.n)
            if v<=self.tol**2:
                n+=1
                self.del_element(e.n)
        if n:
            print("Found and delete %d elements with negative or zero zolume." % (n, ))
            deletedNodes=self.get_unref_nodes()
            self.del_nodes(deletedNodes)
        return deletedNodes

    def getNodesOnLines(self, lines=[]):
        if not hasattr(self, 'allLines'):
            return None
        if not lines:
            lines=list(self.allLines.keys())
        rez=[]
        for l in lines:
            for ll in self.allLines[l]:
                rez.extend(ll)
        return set(rez)

    def uniteParts(self, newPartNum=1, elements=[]):
        if not elements:
            elements=list(self.els.keys())
        for e in elements:
            self.els[e].part=newPartNum
        self.prts=set(self.partsNumbers())

    def partsNumbers(self):
        rez=[]
        for e in list(self.els.values()):
            if not e.part in rez:
                rez.append(e.part)
        return rez

    def getEdgeLengths(self, en):
        e=self.els[en]
        rez=[]
        for ns in e.getEdges():
            rez.append((Point(*self.nds[ns[0]].crds)-Point(*self.nds[ns[1]].crds)).norm())
        return rez

    def getEdgeAngles(self, en):
        e=self.els[en]
        rez=[]
        edges=e.getEdges()
        for ee in e.getEdgesAngles():
            v1=(Point(*self.nds[edges[abs(ee[0])-1][1]].crds)-Point(*self.nds[edges[abs(ee[0])-1][0]].crds)).scaled(np.sign(ee[0]))
            v2=(Point(*self.nds[edges[abs(ee[1])-1][1]].crds)-Point(*self.nds[edges[abs(ee[1])-1][0]].crds)).scaled(np.sign(ee[1]))
            rez.append(v1.angle(v2))
        return rez

    def getBadElements(self, aspect=None, angles=None):
        bad=[]
        for en in list(self.els.keys()):
            if aspect:
                d=self.getEdgeLengths(en)
                asp=max(d)/min(d)
                if asp>aspect:
                    print('Element number: ', en, ' aspect ratio: ', asp)
                    bad.append(en)
                    continue
            if angles:
                a= self.getEdgeAngles(en)
                if min(a)<angles[0] or max(a)>angles[1]:
                    print('Element number: ', en, ' Min, max angles: ', min(a), max(a))
                    bad.append(en)
#        print('Total number of bad elements: ', len(bad))
        return bad

    def getNdDistance(self, n1, n2):
        nd1=self.nds[n1]
        nd2=self.nds[n2]
        return sqrt((nd1.x-nd2.x)**2+(nd1.y-nd2.y)**2+(nd1.z-nd2.z)**2)

    def getNearestNode(self, x,y,z, nds=[]):
        if not nds:
            nds=self.nds.keys()
        d=1e20
        rez=0
        for n in list(nds):
            dd=self.nds[n].distTo(x,y,z)
            if dd<d:
                d=dd
                rez=n
        return rez

    def getNearestElement(self, x,y,z):
        d=1e20
        rez=0
        #N=len(self.els)
        for i, e in enumerate(self.els.values()):
            xx,yy,zz=self.get_element_center(e)
            dd=sqrt((xx-x)**2+(yy-y)**2+(zz-z)**2)
            if dd<d:
                d=dd
                rez=e.n
            #print(progress(i, N-1), end='')
        return rez

    def getExternalNodes(self, parts=[]):
        if np.__version__<'1.13.0':
            print('Error! Version of numpy >=1.13.0 is required!')
        if not parts:
            elements=self.els.values()
        else:
            elements=self.getElementsByParts(parts, values=True)
        allFaces=[]
        for e in elements:
            tmp=e.getFaces()
            for tt in tmp:
                tt.sort()
                tt+=[0]*(4-len(tt))
            allFaces.extend(tmp)
        f, c = np.unique(allFaces, return_counts=True, axis=0)
        rez = set(chain(*f[c==1].tolist()))
        rez-={0}
        return list(rez)


    def externalNodesbyPrepost(self, parts=[], sShift=100):
#        if not parts:
#            parts=self.prts
        self.save_model('tmp.k')
        prepost=lsPrePost()
        command=prepost.openKeywordCommand('tmp.k')
        num=[]
        if not parts:
            command+=prepost.createExternalNodesSetCommand(sShift)
            num=sShift
        else:
            num=list(map(lambda x: x+sShift, parts))
            command+=prepost.createExternalNodesSetByPartCommand(num, parts)
        command+=prepost.saveKeywordCommand('tmp.k')
        prepost.runCommand(command)
        s=read_set('tmp.k', num=num)
        remove('tmp.k')
        return s

    def findUnconnectedRegions2(self, parts=[]):
        if not parts:
            parts=self.prts
        e=self.getElementsByParts(parts)
        regions=[]
        ee=e.pop()
        region=[ee]
        nel=1
        nnds=set(self.els[ee].nodes)
        N=len(e)
        while(len(e)):
            jj=[]
            for i in range(len(e)):
                nndss=set(self.els[e[i]].nodes)
                if nnds.intersection(nndss):
                    nnds.update(nndss)
                    region.append(e[i])
                    jj.append(i)
            jj.reverse()
            for i in jj:
                e.pop(i)
            if len(e)==nel:
                regions.append(region)
                ee=e.pop()
                region=[ee]
                nel=1
                nnds=set(self.els[ee].nodes)
            nel=len(e)
            print(progress(N-nel, N), end='')
        print(progress(N, N))
        regions.append(region)
        return regions

    def findUnconnectedRegions(self, parts: list[int]=[])->list[list[int]]:
        if not parts:
            parts=self.prts
        e=self.getElementsByParts(parts)
        nodes = []
        elements = []
        N = len(e)
        pr = max(N//20, 1)
        for ii, ee in enumerate(e):
            e = self.els[ee]
            found = []
            ns = set(e.nodes)
            for i, n in enumerate(nodes):
                if ns.intersection(n):
                    found.append(i)
            if not found:
                nodes.append(ns)
                elements.append(set([e.n]))
            else:
                nodes[found[0]] = ns.union(nodes[found[0]])
                elements[found[0]].add(e.n)
                if len(found)>1:
                    for i in found[:0:-1]:
                        nodes[found[0]] = nodes.pop(i).union(nodes[found[0]])
                        elements[found[0]] = elements.pop(i).union(elements[found[0]])
            if ii % pr==0:
                print(progress(ii, N), end='')
        print(progress(N, N))
        return elements

    def findNodePairsForTieing(self, parts1=[], parts2=[], tol=1e-6, startSetN=800):
        if not parts1:
            parts1=list(self.prts)
        if not parts2:
            parts2=list(self.prts)
        if type(parts1)==int: parts1=[parts1]
        if type(parts2)==int: parts2=[parts2]
        partsPairs=[]
        for p2 in parts2:
            if p2 in parts1:
                parts1.remove(p2)
            for p1 in parts1:
                partsPairs.append([p1,p2])
        eNodes=self.externalNodesbyPrepost(list(parts1)+list(parts2), startSetN)
        cnstr=[]
        for pair in partsPairs:
            print('Connecting part {} to part {}\n'.format(*pair))
            b1=boundBox(self._getBoundBox(pair[0]))
            b2=boundBox(self._getBoundBox(pair[1]))
            b1.offsetBounds(1.1*tol)
            b2.offsetBounds(1.1*tol)
            ns1=[]
            ns2=[]
            for n in eNodes[pair[0]+startSetN]:
                if self.nds[n].isInsideBox(b2):
                    ns1.append(n)
            for n in eNodes[pair[1]+startSetN]:
                if self.nds[n].isInsideBox(b1):
                    ns2.append(n)
            N=len(ns1)
            for i, n in enumerate(ns1):
                nn=self.getNearestNode(*self.nds[n].crds, nds=ns2)
                if self.nds[n].distTo(*self.nds[nn].crds)<=tol and n!=nn:
                    cnstr.append([n,nn])
                    eNodes[pair[0]+startSetN].remove(n)
                    eNodes[pair[1]+startSetN].remove(nn)
                    ns2.remove(nn)
                print(progress(i+1,N), end='')
            print(progress(N,N))
        coupledCard='*CONSTRAINED_NODE_SET\n{:10d}         7\n'
        i=0
        kw=""
        for c in cnstr:
            i+=1
            kw+='*SET_NODE_LIST\n{:10d}\n'.format(startSetN+i)
            kw+='{:10d}{:10d}\n'.format(*c)
            kw+=coupledCard.format(startSetN+i)
        return cnstr, kw

    def getMeshStatistics(self):
        rez = {}
        rez['number_of_nodes'] = len(self.nds)
        solids, shells = self.getByType()
        rez['number_of_shells'] = len(shells)
        rez['number_of_shells'] = len(shells)
        rez['number_of_solids'] = len(solids)
        rez['number_or_parts'] = len(self.prts)
        min_a = 1e6
        max_a = -1e6
        min_l = 1e6
        max_l = -1e6
        max_ratio = -1e-6
        for en in self.els.keys():
            aa = self.getEdgeAngles(en)
            ll = self.getEdgeLengths(en)
            min_a = min(min_a, min(aa))
            max_a = max(max_a, max(aa))
            min_l = min(min_l, min(ll))
            max_l = max(max_l, max(ll))
            max_ratio = max(max_ratio, max(ll)/min(ll))
        rez['min_edges_angle'] = min_a
        rez['max_edges_angle'] = max_a
        rez['min_edge_length'] = min_l
        rez['max_edge_length'] = max_l
        rez['max_edges_ratio'] = max_ratio
        return rez
    def nodes_to_nprecord(self, nodes=[]):
        if type(nodes)==int:
            nodes=[nodes]
        if not nodes:
            nodes = self.nds.keys()
        dt = np.dtype({'names': ['n', 'x', 'y', 'z'], 'formats': [np.uint8, np.float, np.float, np.float]})
        rez=[]
        for n in nodes:
            nd = self.nds[n]
            rez.append(tuple([n]+nd.crds))
        return np.rec.array(rez, dtype=dt)


def fix2DmeshForSweeping(oldmsh):
    msh=deepcopy(oldmsh)
    for nn in list(msh.nds.keys()):
        msh.nds[nn].crds[2]=0
#    msh.del_nodes(msh.get_unref_nodes())
    allNodes=msh.get_elements_nodes()
    Ncount=max(allNodes)
    nodes_on_axis, elements_on_axis=msh.get_items_on_axis()
    mshecount=max(msh.els.keys())
    for en in elements_on_axis[0]:
        if msh.els[en].nodesCount==4:
            for i in range(4):
                if msh.els[en].nodes[i] in nodes_on_axis: break
            idx=[kk-4 if kk>3 else kk for kk in range(i,i+4)]
            e=msh.els[en]
            tmpnodes=e.uniqNodes
            nns=[tmpnodes[idx[0]], tmpnodes[idx[1]], tmpnodes[idx[3]], tmpnodes[idx[3]]]
            newe=element(en, part=e.part, nodes=nns, etype=e.etype)
            msh.els[en]=newe
            nns=[tmpnodes[idx[1]], tmpnodes[idx[2]], tmpnodes[idx[3]], tmpnodes[idx[3]]]
            newe=element(mshecount+1, part=e.part, nodes=nns, etype='*ELEMENT_SHELL')
            msh.els[mshecount+1]=newe
            mshecount+=1
    msh.nsetSortByY(nodes_on_axis)
    allElsOnAxis=[]
    allElsOnAxis.extend(elements_on_axis[0])
    allElsOnAxis.extend(elements_on_axis[1])
    msh.esetSortByY(allElsOnAxis, nodes_on_axis)
    NN=len(nodes_on_axis)
    NE=len(allElsOnAxis)
    tops=[0]
    bottoms=[NN-1]
    for i in range(1, NN-2):
        found=False
        for j in range(0, NE):
            e=msh.els[allElsOnAxis[j]]
            if e.containsNds([nodes_on_axis[i], nodes_on_axis[i+1]]):
                found=True
                break
        if not found:
            tops.append(i+1)
            bottoms.append(i)
#    print tops, bottoms
    for en in elements_on_axis[0]:
        msh.els[en].nodes=msh.els[en].uniqNodes
        for i in range(3):
            if msh.els[en].nodes[i] in nodes_on_axis: break
        ndind=i
        ndonax=msh.els[en].nodes[i]
        ndinset=nodes_on_axis.index(ndonax)
        elinset=allElsOnAxis.index(en)
        NE=len(allElsOnAxis)
        NN=len(nodes_on_axis)
        if ndinset in tops:
#            print 't', ndinset, ndonax
            searchDiap=list(range(elinset+1,NE))
            nnd=node(Ncount+1, [0, 0.5*(msh.nds[nodes_on_axis[ndinset]].y+msh.nds[nodes_on_axis[ndinset+1]].y),0])
            msh.nds[nnd.n]=nnd
            msh.els[en].nodes.insert(ndind+1, nnd.n)
            elements_on_axis[1].append(en)
            nodes_on_axis.insert(ndinset+1, nnd.n)
            msh.els[en].nodes.pop()
            Ncount+=1
            NN+=1
        elif ndinset in bottoms:
#            print 'b', ndinset
            searchDiap=list(range(0,elinset))
            nnd=node(Ncount+1, [0, 0.5*(msh.nds[nodes_on_axis[ndinset]].y+msh.nds[nodes_on_axis[ndinset-1]].y),0])
            msh.nds[nnd.n]=nnd
            msh.els[en].nodes.insert(ndind, nnd.n)
            elements_on_axis[1].append(en)
            nodes_on_axis.insert(ndinset, nnd.n)
            Ncount+=1
            NN+=1
        else:
#            print 'm', ndinset
            searchDiap=list(range(elinset+1,NE))
            nnd=node(Ncount+1, [0, 0.5*(msh.nds[nodes_on_axis[ndinset]].y+msh.nds[nodes_on_axis[ndinset+1]].y),0])
            ny=0.5*(msh.nds[nodes_on_axis[ndinset]].y+msh.nds[nodes_on_axis[ndinset-1]].y)
            msh.nds[ndonax].crds=[0,ny,0]
            msh.nds[nnd.n]=nnd
            msh.els[en].nodes.insert(ndind+1, nnd.n)
            msh.els[en].nodes.pop()
            elements_on_axis[1].append(en)
            nodes_on_axis.insert(ndinset+1, nnd.n)
            Ncount+=1
            NN+=1
        for iii, tp in enumerate(tops):
            if tp>ndinset:
                tops[iii]+=1
        for iii, bt in enumerate(bottoms):
            if bt>ndinset:
                bottoms[iii]+=1
        for een in searchDiap:
            e=msh.els[allElsOnAxis[een]]
            while ndonax in e.nodes:
                ii=e.nodes.index(ndonax)
                e.nodes[ii]=Ncount
        een = max(msh.els.keys())+1
        new_els = {}
        for e in msh.els:
            try:
                aa = msh.getEdgeAngles(e)
            except:
                print(e)
                continue
            if max(aa)>=180:
                el = msh.els[e]
                if el.nodesCount!=4:
                    continue
                idx=aa.index(max(aa))
                _=shift(el.nodes, idx)
                el.nodes=_[:3]
                new_els[een] = element(een, _[2:]+[_[0]], part=el.part, etype=el.etype)
                een+=1
        msh.els.update(new_els)
    return msh

def sweep2Dmesh(msh, divs=10, fi0=0, dfi=360., segments=[], ndset=[]):
    N=divs
    newN={}
    if type(dfi)==list:
        angles = [fi0]
        for ddfi in dfi:
            angles.append(angles[-1]+ddfi)
        N = len(dfi)
        dfi=sum(dfi)
        if abs(dfi-360.)<1e-3:
            angles.pop()
    else:
        angles=[]
        fi=dfi/N
        diap = list(range(0, N)) if dfi==360 else list(range(0, N+1))
        for k in diap:
            angles.append(fi0+fi*k)
    angles = [a*math.pi/180. for a in angles]
    msh=fix2DmeshForSweeping(msh)
    allNodes=msh.get_elements_nodes()
    Ncount=max(allNodes)
    newModel=model()
    nodes_on_axis, elements_on_axis=msh.get_items_on_axis()
    diap = list(range(0, N)) if dfi==360 else list(range(0, N+1))
    for n in list(msh.nds.keys()):
        if n in nodes_on_axis:
            newModel.nds[n]=node(n, msh.nds[n].crds)
            continue
        for k, a in enumerate(angles):
            x=msh.nds[n].x*math.cos(a)
            y=msh.nds[n].y
            z=msh.nds[n].x*math.sin(a)
            nn=n+k*Ncount
            newModel.nds[nn]=node(nn, [x, y, z])
    nen=1
    alleonax=elements_on_axis[0]+elements_on_axis[1]
    for en in list(msh.els.keys()):
        if not en in alleonax:
            newN[en]=nen
            for k in range(N):
                kn=k+1
                if kn==N and dfi==360: kn=0
                nodes=msh.els[en].nodes
                nc=msh.els[en].nodesCount
                if nc not in [3,4]:
                    continue
                if nc==4:
                    nns=[nodes[0]+k*Ncount, nodes[1]+k*Ncount, nodes[2]+k*Ncount, nodes[3]+k*Ncount,
                         nodes[0]+kn*Ncount, nodes[1]+kn*Ncount, nodes[2]+kn*Ncount, nodes[3]+kn*Ncount]
                if nc==3:
                    nns=[nodes[0]+k*Ncount, nodes[2]+k*Ncount, nodes[2]+kn*Ncount, nodes[0]+kn*Ncount,
                         nodes[1]+k*Ncount, nodes[1]+k*Ncount, nodes[1]+kn*Ncount, nodes[1]+kn*Ncount]
                newModel.els[nen]=element(nen, nodes=nns, part=msh.els[en].part,
                                          etype='*ELEMENT_SOLID')
                nen+=1
    for en in elements_on_axis[1]:
        newN[en]=nen
        cnt=msh.els[en].nodesCount
        for j in range(cnt):
            if msh.els[en].nodes[j] in nodes_on_axis: break
        if not (msh.els[en].nodes[j+1] in nodes_on_axis): j=cnt-1
        idx=[kk-cnt if kk>cnt-1 else kk for kk in range(j,j+cnt)]
        for k in range(N):
            kn=k+1
            if kn==N and dfi==360: kn=0
            nodes=msh.els[en].nodes
            nc=msh.els[en].nodesCount
            if nc not in [3,4]:
                continue
            if nc==4:
                nns=[nodes[idx[1]], nodes[idx[2]]+k*Ncount, nodes[idx[3]]+k*Ncount, nodes[idx[0]],
                     nodes[idx[2]]+kn*Ncount, nodes[idx[2]]+kn*Ncount, nodes[idx[3]]+kn*Ncount, nodes[idx[3]]+kn*Ncount]
            if nc==3:
                nns=[nodes[idx[0]], nodes[idx[1]], nodes[idx[2]]+k*Ncount]+[nodes[idx[2]]+kn*Ncount]*5
            newModel.els[nen]=element(nen, nodes=nns, part=msh.els[en].part,
                                      etype='*ELEMENT_SOLID')
            nen+=1
    newModel.prts=msh.prts
    newModel.nds2Dcount=Ncount
    if segments:
        if type(segments)==list:
            newSegments=[]
            for s in segments:
                if (s[0] in nodes_on_axis) and (s[1] in nodes_on_axis):
                    continue
                if s[0] in nodes_on_axis:
                    for k in range(N):
                        kn=k+1
                        if kn==N and dfi==360:
                            kn=0
                        newSegments.append([s[0], s[1]+k*Ncount, s[1]+kn*Ncount, s[1]+kn*Ncount])
                if s[1] in nodes_on_axis:
                    for k in range(N):
                        kn=k+1
                        if kn==N and dfi==360:
                            kn=0
                        newSegments.append([s[1], s[0]+k*Ncount, s[0]+kn*Ncount, s[0]+kn*Ncount])
                if not s[0] in nodes_on_axis and not s[1] in nodes_on_axis:
                    for k in range(N):
                        kn=k+1
                        if kn==N and dfi==360:
                            kn=0
                        newSegments.append([s[0]+k*Ncount, s[1]+k*Ncount, s[1]+kn*Ncount,
                                            s[0]+kn*Ncount])
        if type(segments)==dict:
            newSegments={}
            for ss in segments:
                newSegments[ss]=[]
                for s in segments[ss]:
                    if (s[0] in nodes_on_axis) and (s[1] in nodes_on_axis):
                        continue
                    if s[0] in nodes_on_axis:
                        for k in range(N):
                            kn=k+1
                            if kn==N and dfi==360:
                                kn=0
                            newSegments[ss].append([s[0], s[1]+k*Ncount, s[1]+kn*Ncount, s[1]+kn*Ncount])
                    if s[1] in nodes_on_axis:
                        for k in range(N):
                            kn=k+1
                            if kn==N and dfi==360:
                                kn=0
                            newSegments[ss].append([s[1], s[0]+k*Ncount, s[0]+kn*Ncount, s[0]+kn*Ncount])
                    if not s[0] in nodes_on_axis and not s[1] in nodes_on_axis:
                        for k in range(N):
                            kn=k+1
                            if kn==N and dfi==360:
                                kn=0
                            newSegments[ss].append([s[0]+k*Ncount, s[1]+k*Ncount, s[1]+kn*Ncount,
                                                s[0]+kn*Ncount])
        newModel.segmentsSet=newSegments

    if ndset:
        if type(ndset)==list:
            newNodes=[]
            for n in ndset:
                if n in nodes_on_axis:
                    newNodes.append(n)
                else:
                    for k in diap:
                        newNodes.append(n+k*Ncount)
        if type(ndset)==dict:
            newNodes={}
            for ns in ndset:
                newNodes[ns]=[]
                for n in ndset[ns]:
                    if n in nodes_on_axis:
                        newNodes[ns].append(n)
                    else:
                        for k in diap:
                            newNodes[ns].append(n+k*Ncount)
        newModel.nodesSet=newNodes

    newModel.newElNumbers=newN

    return newModel

def extrude2Dmesh(msh, N=10., depth=1., normal=2, segments=None, ndset=None):
    N = int(N) 
    if type(depth)==list:
        h = [0]
        for _ in depth:
            h.append(h[-1]+_)
        N = len(depth)
        depth = sum(depth)
    else:
       dz=depth/N
       h = [dz*k for k in range(int(N)+1)]
    Ncount=max(msh.nds.keys())
    newModel=model()
    if depth<0:
        #depth*=-1
        msh.swap2Dorientation()
    for nd in list(msh.nds.values()):
        for k, hh in enumerate(h):
            crds=nd.crds
            crds[normal]+=hh
            newn=node(k*Ncount+nd.n, crds)
            newModel.nds[newn.n]=newn
    nen=1
    for en in list(msh.els.keys()):
        for k in range(int(N)):
            kn=k+1
            nodes=msh.els[en].nodes
            nc=msh.els[en].nodesCount
            if nc==4:
                nns=[nodes[0]+k*Ncount, nodes[1]+k*Ncount, nodes[2]+k*Ncount, nodes[3]+k*Ncount,
                     nodes[0]+kn*Ncount, nodes[1]+kn*Ncount, nodes[2]+kn*Ncount, nodes[3]+kn*Ncount]
            if nc==3:
                nns=[nodes[0]+k*Ncount, nodes[2]+k*Ncount, nodes[2]+kn*Ncount, nodes[0]+kn*Ncount,
                     nodes[1]+k*Ncount, nodes[1]+k*Ncount, nodes[1]+kn*Ncount, nodes[1]+kn*Ncount]
            newModel.els[nen]=element(nen, nodes=nns, part=msh.els[en].part, etype='*ELEMENT_SOLID')
            nen+=1
    newModel.prts=list(msh.prts)
    if segments:
        if type(segments)==list:
            newSegments=[]
            for s in segments:
                for k in range(N):
                    kn=k+1
                    newSegments.append([s[0]+k*Ncount, s[1]+k*Ncount, s[1]+kn*Ncount,
                                        s[0]+kn*Ncount])
        if type(segments)==dict:
            newSegments={}
            for ns in segments:
                newSegments[ns]=[]
                for n in segments[ns]:
                    for k in range(N):
                        kn=k+1
                        newSegments[ns].append([s[0]+k*Ncount, s[1]+k*Ncount, s[1]+kn*Ncount,
                                            s[0]+kn*Ncount])
                        #if depth<0: newSegments[-1].reverse()
        newModel.segmentsSet=newSegments
    if ndset:
        if type(ndset)==list:
            newNodes=[]
            for n in ndset:
                for k in range(N+1):
                    newNodes.append(n+k*Ncount)
        if type(ndset)==dict:
            newNodes = {}
            for ns in ndset:
                newNodes[ns]=[]
                for n in ndset[ns]:
                    for k in range(N+1):
                        newNodes[ns].append(n+k*Ncount)
        newModel.nodesSet=newNodes

    return newModel

def readMSH(mshFileName, etype='solids'):
#    etypes={'triangle': 2, 'tetra': 4, 'prism': 6, 'piramid': 7, 'quad': 3, 'hexa': 5}
    allowedElements={'solids': [4,5,11],
                     'shells': [2,3,9],
                     'lines' : [1, 8],
                     'all'   : [2,3,4,5,9,11]}
    msh=open(mshFileName, 'r')
    lsMesh=model()
    for l in msh:
        if l.startswith('$Nodes'):
            break
    n=int(next(msh))
    for i in range(n):
        l=msh.next().split()
        nn=int(l[0])
        crds=[float(ll) for ll in l[1:]]
        newNode=node(nn, crds)
        lsMesh.nds[nn]=newNode
    for l in msh:
        if l.startswith('$Elements'):
            break
    n=int(next(msh))
    for i in range(n):
        l=msh.next().split()
        en=int(l[0])
        et=int(l[1])
        part=int(l[4])
        nds=[int(ll) for ll in l[5:]]
        if et in allowedElements[etype]:
            if et in allowedElements['solids']:
                elt='*ELEMENT_SOLID'
            if et in allowedElements['shells']:
                elt='*ELEMENT_SHELL'
            if et in allowedElements['lines']:
                elt='*ELEMENT_LINE'
            if et==11:
                nds[9], nds[8] = nds[8], nds[9]
            nel=element(en, nodes=nds, part=part, etype=elt)
            lsMesh.els[en]=nel
            lsMesh.prts.add(part)
    msh.close()
    return lsMesh

def readSegmentsFromMsh(msh_file='mesh.msh', tasktype='2D'):
    segmentTypes={'1D': [15], '2D': [1,8], '3D': [2,3,9,16]}
    segments={}
    msh=open(msh_file, 'r')
    for l in msh:
        if l.startswith('$Elements'):
            break
    n=int(next(msh))
    for i in range(n):
        l=next(msh)
        ll=l.split()
        if int(ll[1]) in segmentTypes[tasktype]:
            key=int(ll[4])
            if key not in segments:
                segments[key]=[]
            segments[key].append(list(map(int, ll[5:])))
    msh.close()
    return segments

def readNodesFromMsh(msh_file='mesh.msh'):
    nodes={}
    msh=open(msh_file, 'r')
    for l in msh:
        if l.startswith('$Elements'):
            break
    n=int(next(msh))
    for i in range(n):
        l=next(msh)
        ll=l.split()
        if ll[1]=='15':
            nodes[int(ll[4])]=int(ll[5])
    msh.close()
    return nodes

def writeSegments(fname, segments, num=1, mode='a', offset=0, end=True):
    from io import TextIOWrapper
    if type(fname) in [TextIOWrapper, file]:
        mode='fileStreem'
        fout=fname
    else:
        fout=open(fname, mode)
    if mode=='w':
        fout.write('*KEYWORD\n')
    fout.write('*SET_SEGMENT_LIST\n%10d\n' % (num,))
    for s in segments:
        s+=[s[-1]]*(4-len(s))
        fs='%10d'*len(s)+'\n'
        if offset:
            s=[nn+offset for nn in s]
        fout.write(fs % tuple(s))
    if end:
        fout.write('*END')
    if mode!='fileStreem':
        fout.close()

def writeNodeSet(fname, nodes, num=1, mode='a', offset=0, end=True):
    fout=open(fname, mode)
    if mode=='w':
        fout.write('*KEYWORD\n')
    fout.write('*SET_NODE_LIST\n%10d\n' % (num,))
    i=0
    for nd in nodes:
        fout.write('%10d' % (nd+offset,))
        i+=1
        if i==8:
            fout.write('\n')
            i=0
    if i!=0:
        fout.write('\n')
    if end:
        fout.write('*END')
    fout.close()

def getSegmentsOnSurface(allSegments, surface):
    if type(surface)==int:
        surface=[surface]
    rez=[]
    for s in surface:
        rez.extend(allSegments[s])
    return rez

def getNodesOnSurface(allSegments, surface):
    rez=[]
    for cs in getSegmentsOnSurface(allSegments, surface):
        rez.extend(cs)
    return list(set(rez))

getSegmentsOnLines=getSegmentsOnSurface
getNodesOnLines=getNodesOnSurface

#---------------------------------------
gmsh_mesher_options='''Geometry.Tolerance = geometryTolerance;
Geometry.OCCSewFaces = 0;
Mesh.Algorithm = meshingAlgorithm;
Mesh.Algorithm3D = 1;
Mesh.AngleSmoothNormals = 30;
Mesh.AnisoMax = 5;
Mesh.AllowSwapAngle = 10;
Mesh.BdfFieldFormat = 1;
Mesh.Binary = 0;
Mesh.Bunin = 0;
Mesh.Lloyd = 0;
Mesh.SmoothCrossField = 0;
Mesh.CgnsImportOrder = 1;
Mesh.ChacoArchitecture = 1;
Mesh.ChacoEigensolver = 1;
Mesh.ChacoEigTol = 0.001;
Mesh.ChacoGlobalMethod = 1;
Mesh.ChacoHypercubeDim = 0;
Mesh.ChacoLocalMethod = 1;
Mesh.ChacoMeshDim1 = 1;
Mesh.ChacoMeshDim2 = 1;
Mesh.ChacoMeshDim3 = 1;
Mesh.ChacoPartitionSection = 1;
Mesh.ChacoSeed = 7654321;
Mesh.ChacoVMax = 250;
Mesh.ChacoParamINTERNAL_VERTICES = 0;
Mesh.ChacoParamREFINE_MAP = 1;
Mesh.ChacoParamREFINE_PARTITION = 0;
Mesh.ChacoParamTERMINAL_PROPOGATION = 0;
Mesh.CharacteristicLengthExtendFromBoundary = 1;
Mesh.CharacteristicLengthFactor = 1;
Mesh.CharacteristicLengthMin = 0;
Mesh.CharacteristicLengthMax = maxElementSize;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.CharacteristicLengthFromPoints = 1;
Mesh.Clip = 0;
Mesh.ColorCarousel = 1;
Mesh.CpuTime = 18.53560638427734;
Mesh.DrawSkinOnly = 0;
Mesh.Dual = 0;
Mesh.ElementOrder = 1;
Mesh.Explode = 1;
Mesh.FlexibleTransfinite = 0;
Mesh.NewtonConvergenceTestXYZ = 0;
Mesh.Format = 10;
Mesh.Hexahedra = 1;
Mesh.HighOrderNumLayers = 6;
Mesh.HighOrderOptimize = 0;
Mesh.HighOrderPoissonRatio = 0.33;
Mesh.HighOrderThresholdMin = 0.1;
Mesh.HighOrderThresholdMax = 2;
Mesh.HighOrderOptPrimSurfMesh = 0;
Mesh.LabelSampling = 1;
Mesh.LabelType = 0;
Mesh.LcIntegrationPrecision = 1e-09;
Mesh.Light = 1;
Mesh.LightLines = 1;
Mesh.LightTwoSide = 1;
Mesh.Lines = 0;
Mesh.LineNumbers = 0;
Mesh.LineWidth = 1;
Mesh.MeshOnlyVisible = 0;
Mesh.MetisAlgorithm = 1;
Mesh.MetisEdgeMatching = 3;
Mesh.MetisRefinementAlgorithm = 3;
Mesh.MinimumCirclePoints = 7;
Mesh.MinimumCurvePoints = 3;
Mesh.MshFileVersion = 2.2;
Mesh.MshFilePartitioned = 0;
Mesh.PartitionHexWeight = 1;
Mesh.PartitionPrismWeight = 1;
Mesh.PartitionPyramidWeight = 1;
Mesh.PartitionQuadWeight = 1;
Mesh.PartitionTetWeight = 1;
Mesh.PartitionTriWeight = 1;
Mesh.NbHexahedra = 0;
Mesh.NbNodes = 473;
Mesh.NbPartitions = 1;
Mesh.NbPrisms = 0;
Mesh.NbPyramids = 0;
Mesh.NbQuadrangles = 176;
Mesh.NbTetrahedra = 0;
Mesh.NbTriangles = 14;
Mesh.Normals = 0;
Mesh.NumSubEdges = 2;
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 0;
Mesh.Partitioner = 2;
Mesh.Points = 0;
Mesh.PointNumbers = 0;
Mesh.PointSize = 4;
Mesh.PointType = 0;
Mesh.Prisms = 1;
Mesh.Pyramids = 1;
Mesh.Quadrangles = 1;
Mesh.QualityInf = 0;
Mesh.QualitySup = 0;
Mesh.QualityType = 2;
Mesh.RadiusInf = 0;
Mesh.RadiusSup = 0;
Mesh.RandomFactor = 1e-09;
Mesh.IgnorePartitionBoundary = 0;
Mesh.RecombinationAlgorithm = 0;
Mesh.RecombineAll = 1;
Mesh.Recombine3DAll = 0;
Mesh.DoRecombinationTest = 1;
Mesh.RecombinationTestHorizStart = 1;
Mesh.RecombinationTestNoGreedyStrat = 0;
Mesh.RecombinationTestNewStrat = 0;
Mesh.RemeshAlgorithm = 0;
Mesh.RemeshParametrization = 0;
Mesh.RefineSteps = 2;
Mesh.Remove4Triangles = 1;
Mesh.ReverseAllNormals = 0;
Mesh.SaveAll = 0;
Mesh.SaveElementTagType = 1;
Mesh.SaveParametric = 0;
Mesh.SaveGroupsOfNodes = 0;
Mesh.ScalingFactor = 1;
Mesh.SecondOrderExperimental = 0;
Mesh.SecondOrderIncomplete = 0;
Mesh.SecondOrderLinear = 0;
Mesh.Smoothing = 1;
Mesh.SmoothNormals = 1;
Mesh.SmoothRatio = 1.8;
Mesh.SubdivisionAlgorithm = 0;
Mesh.SurfaceEdges = 1;
Mesh.SurfaceFaces = 0;
Mesh.SurfaceNumbers = 0;
Mesh.SwitchElementTags = 0;
Mesh.Tangents = 0;
Mesh.Tetrahedra = 1;
Mesh.ToleranceEdgeLength = 0;
Mesh.ToleranceInitialDelaunay = 1e-08;
Mesh.Triangles = 1;
Mesh.VolumeEdges = 1;
Mesh.VolumeFaces = 0;
Mesh.VolumeNumbers = 0;
Mesh.Voronoi = 0;
Mesh.ZoneDefinition = 0;
'''
def mesh2DQuadDominant(cadFile, maxsize=1e10, tol=1e-6, gmshPath="",
                       meshingAlgorithm=9, plane=True, sew=True):
    mesher_options=gmsh_mesher_options
    mesher_options=mesher_options.replace('meshingAlgorithm', str(meshingAlgorithm))
    mesher_options=mesher_options.replace('maxElementSize', str(maxsize))
    mesher_options=mesher_options.replace('geometryTolerance', str(tol))
    f=open('mesh.geo', 'w')
    f.write(mesher_options)
    f.write('Merge \'%s\';\n' % (cadFile,))
    f.close()
    os.system(gmshPath+'gmsh.exe mesh.geo -2')
    mesh2D=readMSH('mesh.msh', etype='shells')
    if plane:
        mesh2D.fix2Dorientation()
    allPoints=readNodesFromMsh()
    allLines=readSegmentsFromMsh()
    os.remove('mesh.geo')
    os.remove('mesh.msh')
    mesh2D.allPoints=allPoints
    mesh2D.allLines=allLines
    print("Removing temporary files...")
    for f in os.listdir(os.curdir):
        if f[-4:].upper()=='.POS':
            os.remove(f)
    if sew:
        mesh2D.sewDuplicateNodes(tol)
    return mesh2D

def mesh2Dsweep(cadFile, segmentsLines=[], nodesLines=[], maxsize=1e10,
                angleDivs=40, tol=1e-6, gmshPath="", splLines=[], sew=True):
    mesh2D = mesh2DQuadDominant(cadFile, maxsize=maxsize, tol=tol,
                                gmshPath=gmshPath, sew=sew)
    if splLines:
        splitElements(mesh2D, lines=splLines)
    pointNodes=mesh2D.allPoints
    allSegments=mesh2D.allLines
    mesh2D.tol=tol
    if segmentsLines==set([]):
        segmentsLines=set(allSegments.keys())
    if nodesLines==set([]):
        nodesLines=set(allSegments.keys())
    if type(segmentsLines)==list:
        partSegments=getSegmentsOnLines(allSegments, segmentsLines)
        partSegments=mesh2D.fixSegmentsOrientations(partSegments)
    if type(segmentsLines)==set:
        partSegments={}
        for ps in segmentsLines:
            partSegments[ps]=mesh2D.allLines[ps]
            partSegments[ps]=mesh2D.fixSegmentsOrientations(partSegments[ps])
    if type(nodesLines)==list:
        partNodes=getNodesOnLines(allSegments, nodesLines)
    if type(nodesLines)==set:
        partNodes={}
        for pn in nodesLines:
            partNodes[pn]=list(set(chain(*mesh2D.allLines[pn])))
    rez=sweep2Dmesh(mesh2D, divs=angleDivs,
                    segments=partSegments,
                    ndset=partNodes)
    rez.pointNodes=pointNodes
    rez.mesh2D=mesh2D
    return rez

def mesh3Dtet(cadFile, segmentsLines=[], nodesLines=[], maxsize=0.5e-3,
              tol=1e-6, gmshPath=""):
    gmshTetCommand="{0} - {1} -o mesh.msh -3 -tol {2} -clmax {3}"
    os.system(gmshTetCommand.format(gmshPath+'gmsh.exe',
                                    cadFile, tol, maxsize))
    mesh3D=readMSH("mesh.msh")
    allSegments=readSegmentsFromMsh(msh_file='mesh.msh', tasktype='3D')
    allLines=readSegmentsFromMsh(msh_file='mesh.msh', tasktype='2D')
    allPoints=readNodesFromMsh()
    os.remove('mesh.msh')
    mesh3D.allPoints=allPoints
    mesh3D.allLines=allLines
    mesh3D.allSegments=allSegments
    return mesh3D

def splitElements(mesh, lines=[]):
    lNodes=mesh.getNodesOnLines(lines)
    onePlyElements=[]
    for e in list(mesh.els.values()):
        if set(e.nodes).issubset(lNodes):
            onePlyElements.append(e.n)
    mesh.esetSortByCentrY(onePlyElements)
    nn=max(mesh.nds.keys())+1
    ee=max(mesh.els.keys())+1
    e1n=onePlyElements[0]
    e1=mesh.els[e1n]
    fcs=e1.getFaces()
    for i, f in enumerate(fcs):
        if set(f).issubset(set(mesh.els[onePlyElements[1]].nodes)):
            break
    fn=[ii%4 for ii in [i, i+1, i+2, i+3]]
    mesh.nds[nn]=node(nn, mesh.get_element_center(e1n))
    mesh.nds[nn+1]=node(nn+1, mesh.getElementCenterByNodes(e1n, fcs[i]))
    mesh.els[e1n].nodes=fcs[fn[3]]+[nn+1, nn]
    mesh.els[ee]=element(ee, fcs[fn[1]]+[nn, nn+1], part=e1.part, etype=e1.etype)
    mesh.els[ee+1]=element(ee+1, fcs[fn[2]]+[nn], part=e1.part, etype=e1.etype)
    nn+=1
    ee+=2
    N=len(onePlyElements)
    for j in range(1, N-1):
        e1n=onePlyElements[j]
        e1=mesh.els[e1n]
        nns=e1.nodes
        fcs=e1.getFaces()
        for i, f in enumerate(fcs):
            if set(f).issubset(set(mesh.els[onePlyElements[j+1]].nodes)):
                break
        fn=[ii%4 for ii in [i, i+1, i+2, i+3]]
        mesh.nds[nn+1]=node(nn+1, mesh.getElementCenterByNodes(e1n, fcs[i]))
        mesh.els[e1n].nodes=fcs[fn[3]]+[nn+1, nn]
        mesh.els[ee]=element(ee, fcs[fn[1]]+[nn, nn+1], part=e1.part, etype=e1.etype)
        ee+=1
        nn+=1
    e1n=onePlyElements[-1]
    e1=mesh.els[e1n]
    fcs=e1.getFaces()
    for i, f in enumerate(fcs):
        if set(f).issubset(set(nns)):
            break
    fn=[ii%4 for ii in [i, i+1, i+2, i+3]]
    mesh.nds[nn+1]=node(nn+1, mesh.get_element_center(e1n))
    mesh.els[e1n].nodes=fcs[fn[1]]+[nn+1, nn]
    mesh.els[ee]=element(ee, fcs[fn[3]]+[nn, nn+1], part=e1.part, etype=e1.etype)
    mesh.els[ee+1]=element(ee+1, fcs[fn[2]]+[nn+1], part=e1.part, etype=e1.etype)

def gmshMesh(cadFile, maxsize=0.5e-3, minsize=1e-12, tol=1e-6, gmshPath="", order=1, dim=2,
             smoothingSteps=0, recombine=False, curvSize=False):
    gmshCommand="{0} - {1} -o mesh.msh -{2} -tol {3} -clmax {4} -smooth {5} -order {6} -clmin {7}".format(
        os.path.join(gmshPath,'gmsh.exe'), cadFile, dim, tol, maxsize, smoothingSteps, order, minsize)
    if curvSize:
        gmshCommand+=" -clcurv"
    if dim==2 and recombine:
        gmshCommand+=' -string \"Mesh.RecombineAll = 1;\" -algo delquad'
    os.system(gmshCommand)
    et='solids' if dim==3 else 'shells'
    mesh=readMSH("mesh.msh", et)
    allLines=readSegmentsFromMsh(msh_file='mesh.msh', tasktype='2D')
    allPoints=readNodesFromMsh()
    if dim==3:
        allSegments=readSegmentsFromMsh(msh_file='mesh.msh', tasktype='3D')
    else:
        allSegments=allLines
    os.remove('mesh.msh')
    mesh.allPoints=allPoints
    mesh.allLines=allLines
    mesh.allSegments=allSegments
    return mesh

def planeBoxMesh(box=((0,0), (1,1)), esize=1.0, part=1):
    def calcNumDeltas(dx, esize):
        num=int(max(round(dx/esize),1))
        dn=dx/num
        return num, dn

    ix, dx=calcNumDeltas(box[1][0]-box[0][0], esize)
    iy, dy=calcNumDeltas(box[1][1]-box[0][1], esize)
    m=model()
    for j in range(iy+1):
        for i in range(ix+1):
            n=node(i+1+j*(ix+1), [box[0][0]+dx*i, box[0][1]+dy*j, 0])
            m.nds[n.n]=n
    for j in range(iy):
        for i in range(ix):
            e=element(i+1+j*ix, part=part, nodes=[i+1+j*(ix+1), i+2+j*(ix+1), i+2+(j+1)*(ix+1), i+1+(j+1)*(ix+1)], etype='*ELEMENT_SHELL')
            m.els[e.n]=e
    return m

def externalSegmentsForAleFill(cad, elsize, tol=1e-6, offset=0,
                                noffset=0):
    """Функция создает наборы внешних сегментов на основе 2D или 3D cad модели
        В случае 2D геометрии она прокручивается относительно оси OY.
        Возвращает экземпляр класса model, содержащий узлы и наборы сегментов в
        model.segmentSets

    Arguments:
        cad {str} -- путь к cad файлу (step, iges)
        elsize {float} -- максимальный размер конечного элемента

    Keyword Arguments:
        tol {float} -- геометрический допуск (default: {1e-6})
        offset {int} -- смещение нумерации наборов (default: {0}). Для части с номером N
        будет создан набор сегментов с номером N+offset - внешняя поверхность части
        noffset {int} -- смещение нумерации узлов (default: {0}). Узлы в сегментах будут
        иметь номера: исходный номер узла+noffset
    """
    from math import pi, floor
    m = gmshMesh(cad, dim=3, maxsize=elsize, tol=tol, smoothingSteps=2)
    if not m.els:
        lines = []
        for l in m.allLines:
            n1 = m.allLines[l][0][0]
            n2 = m.allLines[l][0][1]
            if abs(m.nds[n1].x) >= tol or abs(m.nds[n2].x) >= tol:
                lines.append(l)
        bb = m.boundBox
        N = int(floor(2*pi*bb['maxX']/elsize))
        m = mesh2Dsweep(cad, segmentsLines=lines,
                           maxsize=elsize, angleDivs=N, tol=tol)
        for i, s in enumerate(m.segmentsSet):
            if len(s)!=len(set(s)):
                m.segmentsSet[i] = s[:-1]
                s=s[:-1]
        m.segmentSets[500] = m.segmentsSet
    else:
        m.segmentSets[500] = list(chain(*m.allSegments.values()))
#    m.save_model('test.k')
    for p in m.prts:
        ss = deepcopy(m.segmentSets[500])
        seg = []
        els = m.getElementsByParts(p)
        nodes = m.get_elements_nodes(els)
        for sss in ss:
            seg_on_part = True
            for n in sss:
                seg_on_part = seg_on_part and (n in nodes)
            if seg_on_part:
                seg.append(sss)
        seg = m.fixSegmentsOrientations(seg, p)
        m.segmentSets[offset+p] = seg
    m.segmentSets.pop(500)
    rez = model()
    nn = noffset+1
    for p in m.segmentSets:
        newnn={}
        rez.segmentSets[p]=[]
        for sss in m.segmentSets[p]:
            for nd in sss:
                if not nd in newnn:
                    newnn[nd]= nn
                    rez.nds[nn] = node(nn, m.nds[nd].crds)
                    nn+=1
            rez.segmentSets[p].append([newnn[nd] for nd in sss])
    return rez

def sortNodesForPlaneElement(nodes):
    dtype = np.dtype(
        {
            'names': ['n', 'x', 'y', 'a', 'r'],
            'formats': [np.int32, np.float64, np.float64, np.float64, np.float64]
        }
    )
    data = [(n.n, n.x, n.y, 0, 0) for n in nodes]
    data = np.array(data, dtype=dtype)
    data.sort(order='y')
    vx = data['x']-data['x'][0]
    vy = data['y']-data['y'][0]
    data['r'] = vx*vx+vy*vy
    data['a'] = np.arctan2(vy, vx)
    data[1:].sort(order=['a', 'r'])
    return data['n']

def create_cohesive_layer(in_model, out_model, cohesive_part=5, shrink=0, estart=0, nstart=0, part=None,
                          fillCorners=False, tol=1e-6,
                          symm=None, symm_tol=0):
    edges = defaultdict(lambda: [])
    m1 = model(in_model)
    m2 = model()
    print('inserting layer...')
    if not estart:
        estart = 1+max(m1.els.keys())
    if not nstart:
        nstart = 1+max(m1.nds.keys())
    splitted_nodes = defaultdict(lambda: [])
    splitted_edges = defaultdict(lambda: [])
    if part == None:
        all_els = m1.els.values()
    else:
        all_els = m1.getElementsByParts(part, values=True)
    nn = nstart
    en = estart
    for e in all_els:
        ns = []
        ec = m1.get_element_center(e.n)
        nnn = {}
        for n in e.nodes:
            nd = m1.nds[n]
            ns.append(nn)
            if shrink:
                m2.nds[nn] = node(nn, [ec[0] + shrink * (nd.x - ec[0]),
                                       ec[1] + shrink * (nd.y - ec[1]),
                                       ec[2] + shrink * (nd.z - ec[2])])
            else:
                m2.nds[nn] = node(nn, nd.crds)
            splitted_nodes[n].append(nn)
            nnn[n] = nn
            nn += 1
        el = element(en, ns, e.part, e.etype)
        for f in e.getFaces():
            ff = f.copy()
            ff.sort()
            fs = '-'.join(['{}'] * len(f))
            edges[fs.format(*ff)].append([nnn[fff] for fff in f])
        m2.els[en] = el
        en += 1
        if e.etype.upper()=='*ELEMENT_SOLID':
            for f in e.getEdges():
                ff = f.copy()
                ff.sort()
                fs = '-'.join(['{}'] * len(f))
                splitted_edges[fs.format(*ff)].append([nnn[fff] for fff in f])
    for es in edges.values():
        if len(es) < 2:
            continue
        if e.etype == '*ELEMENT_SHELL':
            ns = es[0] + es[1]
            ns.reverse()
        if e.etype == '*ELEMENT_SOLID':
            ff = es[1].copy()
            ff.reverse()
            nnn = m2.getNearestNode(*m2.nds[es[0][0]].crds, es[1])
            ff = shift(ff, ff.index(nnn))
            ns = es[0] + ff
        m2.els[en] = element(en, ns, cohesive_part, e.etype)
        en += 1
        
    if fillCorners:
        els_to_sort = []
        sorter = {'*ELEMENT_SHELL': sortNodesForPlaneElement, '*ELEMENT_SOLID': partial(sort_8_nodes_for_hex, tol=tol)}[e.etype.upper()]
        min_nodes_count = {'*ELEMENT_SHELL': [4], '*ELEMENT_SOLID': [8]}
        print('filling conners...')
        for n_old, ns_new in splitted_nodes.items():
            if len(ns_new) in min_nodes_count[e.etype.upper()]:
                m2.els[en] = element(en, ns_new, cohesive_part, e.etype)
                els_to_sort.append(en)
                en += 1
        print('filling edges...')
        for se in splitted_edges.values():
            if len(se)==4:
                m2.els[en] = element(en, list(chain(*se)), cohesive_part, e.etype)
                els_to_sort.append(en)
                en+=1
        print('sorting nodes...')
        p = ConsoleProgress(len(els_to_sort))
        for i, e in enumerate(els_to_sort):
            nds = m2.els[e].nodes
            nds = sorter([m2.nds[n] for n in nds])
            m2.els[e].nodes = nds
            p.current=i
            p.redraw()
    p.finalize()
    if symm!=None and symm_tol:
        print('processing symmetry plane...')
        for n in m2.nds:
            if abs(m2.nds[n].crds[symm])<=symm_tol:
                crds = m2.nds[n].crds
                crds[symm]=0
                m2.nds[n].crds = crds
    prts = m1.partsNumbers()
    prts.remove(part)
    els = m1.getElementsByParts(prts)
    nds = m1.get_elements_nodes(els)
    for n in nds:
        m2.nds[n] = m1.nds[n]
    for e in els:
        m2.els[e] = m1.els[e]
    m2.cards = m1.cards
    m2.save_model(out_model)

def sort_4_nodes(nodes, direction_vector):
    dtype = np.dtype(
        {
            'names': ['n', 'a'],
            'formats': [np.int32, np.float64]
        }
    )
    centr = np.zeros(3)
    for n in nodes:
        centr += np.array(n.crds)
    centr /= 4.
    c = Point(*centr)
    v = []
    v.append(Point(*nodes[0].crds)-c)
    v.append(Point(*nodes[1].crds)-c)
    v.append(Point(*nodes[2].crds)-c)
    v.append(Point(*nodes[3].crds)-c)
    a = []
    n = []
    for i in range(1, 4):
        n.append(nodes[i].n)
        a.append(v[0].angle(v[i], direction_vector))
    data = np.array(list(zip(n, a)), dtype=dtype)
    data.sort(order='a')
    return [nodes[0].n]+list(data['n'])

def sort_8_nodes_for_hex(nodes, tol=1e-6):
    def plane(point, normal):
        return lambda x, y, z: (x-point.x)*normal.x+(y-point.y)*normal.y+(z-point.z)*normal.z
    nodes = {n.n: n for n in nodes}
    success = False
    numbers = nodes.keys()
    for pnts in permutations(numbers, 4):
        pp = [Point(*nodes[p].crds) for p in pnts]
        nrm = (pp[1]-pp[0]).cross(pp[2]-pp[0]).normalized()
        pl1 = plane(pp[0], nrm)
        if abs(pl1(*pp[3].data)) <= tol:
            pnts2 = list(set(numbers) - set(pnts))
            check = [pl1(*nodes[p].crds) for p in pnts2]
            if np.all(np.array(check) > 0):
                normal = 1.
            elif np.all(np.array(check) < 0):
                normal = -1.
            else:
                normal = 0
            if normal != 0:
                success = True
                break
    if not success:
        print('Can\'t find 4 nodes in one plane')
        return None
    pnts = sort_4_nodes([nodes[p] for p in pnts], nrm.scaled(normal))
    pp = [Point(*nodes[p].crds) for p in pnts]
    pp2 = [Point(*nodes[p].crds) for p in pnts2]
    # if len(pnts2)==2:
    #     pp_ = Point(*nodes[pnts2[0]].crds)
    #     nrm = (pp[0]-pp_).cross(pp[3]-pp_).normalized()
    #     pl_ = plane(pp_, nrm)
    #     check = np.array( [pl_(*p_.data) for p_ in [pp[1], pp[2], pp2[1]] ])
    #     if not(np.all(check>0) or np.all(check<0)):
    #         pnts2.reverse()
    #     return pnts+pnts2
    success = False
    for pnts3 in permutations(pnts2, 2):
        pp3 = [Point(*nodes[p].crds) for p in pnts3]
        nrm = (pp[1]-pp[0]).cross(pp3[0]-pp[0]).normalized()
        pl3 = plane(pp[0], nrm)
        if abs(pl3(*pp3[1].data)) <= tol:
            pnts4 = [pnts[3], pnts[2]]+list(set(pnts2) - set(pnts3))
            check = [pl3(*nodes[p].crds) for p in pnts4]
            if np.all(np.array(check) > 0):
                normal = 1.
            elif np.all(np.array(check) < 0):
                normal = -1.
            else:
                normal = 0
            if normal != 0:
                success = True
                break
    if not success:
        print('increase tolerance...')
        return None      
    pnts3 = sort_4_nodes([nodes[pnts[0]], nodes[pnts[1]], nodes[pnts3[0]], nodes[pnts3[1]]],
        (pp[1]-pp[0]).cross(pp3[0]-pp[0]))[2:]
    pnts3.reverse()
    idx = pnts2.index(pnts3[0])
    pnts2[0], pnts2[idx] = pnts2[idx], pnts2[0]
    nrm2 = (pp2[1]-pp2[0]).cross(pp2[2]-pp2[0])
    normal = 1 if plane(pp2[0], nrm2)(*pp[0].data)<0 else -1
    pnts2 = sort_4_nodes([nodes[p] for p in pnts2], nrm2.scaled(normal))
    # pp4 = [Point(*nodes[p].crds) for p in pnts4]
    # pnts4 = sort_4_nodes([nodes[pnts[3]], nodes[pnts[2]], nodes[pnts4[2]], nodes[pnts4[3]]],
    #     (pp[2]-pp[3]).cross(pp4[2]-pp[3]))[2:]

    return pnts+pnts2

if __name__=='__main__':
#    m = model('5.k')
#    print(list(m.els[300968].nodes))
#    nns = sort_8_nodes_for_hex(list(m.nds.values()))
#    print(nns)
    # nodes = [
    #     node(1, [0,0,0]),
    #     node(2, [1,0,0]),
    #     node(3,[1,1,0]),
    #     node(4,[0,1,0]),
    #     node(5, [0,0,1]),
    #     node(6, [1,0,1]),
    #     node(7,[1,1,1]),
    #     node(8,[0,1,1]),
    # ]
    # import random
    # random.shuffle(nodes)
    # nds = sort_8_nodes_for_hex(nodes)
    # print(nds)
    create_cohesive_layer('solid-model.k', 'korp+coh.k',
                shrink=0.9, part=1, cohesive_part=11,
                nstart=200000, estart=200000, fillCorners=True, tol=0.03,
                symm=2, symm_tol=1)
    # from random import shuffle
    # m = model('test_mesh.k')
    # for e in m.els:
    #     nds = m.els[e].nodes
    #     shuffle(nds)
    #     m.els[e].nodes = nds
    # m.save_model('sh_mesh.k')
    # for e in m.els:
    #     nds = m.els[e].nodes
    #     nds = sort_8_nodes_for_hex([m.nds[n] for n in nds])
    #     m.els[e].nodes = nds

    # m.save_model('tst.k')
#    en = 256153
#    m = model('1.k')
#    e = m.els[en]
#    print(e.nodes)
#    print(sort_8_nodes_for_hex([m.nds[n] for n in e.nodes], tol=1e-5))