
########## cartisian PDB manipulations #############
### usage "manpdb inputpdb outpdb rotate/trans axis theta/length"##
### For Example: #####
###python manpdb.py tmp.pdb tmprot.pdb rotate 1,0,0 180 #

import math
from math import pi
import itertools
import numpy as np

def square(vec):
    return vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]

def root(vec,sqrt=math.sqrt):
    sqr = vec[0]**2 + vec[1]**2 + vec[2]**2
    return sqrt(sqr)

def unit_vec(vec):
    sqr = root(vec);
    vec[0] = vec[0]/sqr;
    vec[1] = vec[1]/sqr;
    vec[2] = vec[2]/sqr;
    return vec
    
def unit_norm(vec1,vec2):
    norm = [0,0,0];
    norm[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
    norm[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
    norm[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
    norm = unit_vec(norm);
    return norm

def rot_vec(vec,axis,theta):
    "rotate vec along an axis"
    rotmat = [[0,0,0],
              [0,0,0],
              [0,0,0]]
   
    rotmat[0][0] = math.cos(theta) + axis[0]*axis[0]*(1-math.cos(theta));
    rotmat[0][1] = axis[0]*axis[1]*(1-math.cos(theta)) - axis[2]*math.sin(theta)
    rotmat[0][2] = axis[0]*axis[2]*(1-math.cos(theta)) + axis[1]*math.sin(theta)
    rotmat[1][0] = axis[0]*axis[1]*(1-math.cos(theta)) + axis[2]*math.sin(theta)
    rotmat[1][1] = math.cos(theta) + axis[1]*axis[1]*(1-math.cos(theta));
    rotmat[1][2] = axis[1]*axis[2]*(1-math.cos(theta)) - axis[0]*math.sin(theta)
    rotmat[2][0] = axis[0]*axis[2]*(1-math.cos(theta)) - axis[1]*math.sin(theta)
    rotmat[2][1] = axis[1]*axis[2]*(1-math.cos(theta)) + axis[0]*math.sin(theta)
    rotmat[2][2] = math.cos(theta) + axis[2]*axis[2]*(1-math.cos(theta));

    v_new = [0,0,0];
    v_new[0] = rotmat[0][0]*vec[0] + rotmat[0][1]*vec[1] + rotmat[0][2]*vec[2];
    v_new[1] = rotmat[1][0]*vec[0] + rotmat[1][1]*vec[1] + rotmat[1][2]*vec[2];
    v_new[2] = rotmat[2][0]*vec[0] + rotmat[2][1]*vec[1] + rotmat[2][2]*vec[2];
    
    return v_new

def trans_vec(vec,unit_vec,theta):
    vec_new = [vec[s]+unit_vec[s]*theta for s in xrange(3)]
    return vec_new

import sys

if __name__ == '__main__':
    usage = "manpdb inputpdb outpdb rotate/trans axis theta/length"
    if len(sys.argv) < 6:
        print (usage)
        sys.exit(0)

    p = open(sys.argv[1],'r')
    lines = p.readlines()
    p.close()
    p = open(sys.argv[2],'w')


    #compute com for rotation first#
    if sys.argv[3] == 'rotate':
        x=0;y=0;z=0
        nat = 0
        for i in range(len(lines)):
            if not lines[i].startswith("ATOM  "):
                pass
            else:
                vec = [float(s) for s in lines[i][30:55].split()]
                x += vec[0]; y += vec[1]; z += vec[2]
                nat += 1
        com = [x/nat,y/nat,z/nat]

    for i in range(len(lines)):
        if lines[i][0:4]!="ATOM":
            p.write(lines[i])
        else:
            axis = [float(s) for s in sys.argv[4].split(',')]
            axis = unit_vec(axis)

            if sys.argv[3] == 'rotate':
                vec = [float(s) for s in lines[i][30:55].split()]
                vec = [vec[s]-com[s] for s in range(3)]
                newvec = rot_vec(vec,axis,float(sys.argv[5])/180*pi)
                newvec = [newvec[s]+com[s] for s in range(3)]
            elif sys.argv[3] == 'trans':
                vec = [float(s) for s in lines[i][30:55].split()]
                newvec = trans_vec(vec,axis,float(sys.argv[5]))
            else:
                print ("only roate and trans mode are supported")
            x,y,z = newvec
            new = lines[i][0:30] + "%8.3f%8.3f%8.3f"%(x,y,z) + lines[i][54:]
            p.write(new)
    p.close()


