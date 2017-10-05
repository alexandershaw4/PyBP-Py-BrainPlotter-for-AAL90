#!/usr/bin/env python2
'''
Example user script for PyBP

AS
'''

import sys
sys.path.append("/Users/Alex/code/PyBP/")

from PyBP import *
import numpy as np
import os


# get template / example files
thisdir = os.path.dirname(os.path.realpath(__file__))


# both overlay & network
#-----------------------
g   = thisdir + "/inflated_8k.gii"  # the brain surface (gii)
v,f = GetMesh(g)                    # get vertices & faces

# load a t-vec of length 90
funcfile = thisdir + "/TVec.txt"
o = np.loadtxt(funcfile)
y = alignoverlay(v,f,o) 

AALv,L   = GetAAL()  # fetch AAL vertices & labels
file     = thisdir + "/exnet.edge"
A        = np.loadtxt(file)  # load .edge file (90x90 matrix)
xx,yy,vv = CtoN(A,AALv)

ax, p3dcollec, fig = PlotMeshwOverlay(v,f,y,.05)
PlotNet(xx,yy,vv,ax,fig)
fig.show()


#from matplotlib import animation
#def rotate(angle):
#    ax.view_init(azim=angle)
#
#rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0,362,2),interval=50)
#rot_animation.save('rotation1.gif', dpi=150, writer='imagemagick')

#
## Read a gifti surface and plot
##-------------------------------
#
#g = thisdir + "/spm.surf.gii" # the brain surface (gii) // try 'inflated_8k.gii'
#v,f = GetMesh(g)   # get vertices & faces
#ax, p3dcollec, fig = PlotMesh(v,f)
#
#
#
## read .edge file (90x90 .csv)
##-------------------------------
#
#g   = thisdir + "/inflated_8k.gii" # the brain surface (gii) 
#v,f = GetMesh(g)   # get vertices & faces
#AALv,L = GetAAL()  # fetch AAL vertices & labels
#
## Read .edge file
#file = "/Users/Alex/Dropbox/KET-PMP-GABZOL/GMM_Paired/thresh_50/KET_Alpha_signal.edge"
#A = np.loadtxt(file)
#
#
## (OR) Generate toy 90x90 connectivity matrix
#A = GenTestNet()
#
#
## convert matrix to sets of AAL vertices
#xx,yy,vv = CtoN(A,AALv)
#
#ax, p3dcollec, fig = PlotMesh(v,f)
#PlotNet(xx,yy,vv,ax,fig)
#fig.show()
#
## save output
##fig.savefig('testfig.png', format='png', dpi=800)
#
#
#
## plot an interpolated 90 element array as functional overlay
##--------------------------------------------------------------
#
#g = thisdir + "/inflated_8k.gii" # the brain surface (gii)
#v,f = GetMesh(g)                 # get vertices & faces
#
## load a t-vec of length 90
#funcfile = "/Users/Alex/Desktop/TVec.txt"
#o = np.loadtxt(funcfile)
#
## (or generate a test overlay)
#o = GenTestOverlay()       # Use a pretend 90-element overlay
#
#
#y = alignoverlay(v,f,o)    # map vals to brain mesh using ICP
#PlotMeshwOverlay(v,f,y,.2) 
#






