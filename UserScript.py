#!/usr/bin/env python3
'''
Example user script for PyBP

AS
'''

import sys
sys.path.append("/Users/Alex/code/PyBP/")

import PyBP
import numpy as np
import os


# define path where all the data are
#------------------------------------------------------------------------------
#thisdir = os.path.dirname(os.path.realpath(__file__))
thisdir = "/Users/Alex/Dropbox/SourceMesh_TestDatas"

# Get the mesh
#------------------------------------------------------------------------------
g   = thisdir + "/BrainMesh_Ch2.gii"  # the brain surface (gii)
v,f = PyBP.GetMesh(g)                      # get vertices & faces

# define the source positions from an atlas
#------------------------------------------------------------------------------
AALv,AALf,AALi = PyBP.definesources('/Users/Alex/Dropbox/AAL_Template/AAL_Labels.nii')

# load ACTIVITY: a vec of length 90 (matching regions in ROI nifti) 
#------------------------------------------------------------------------------
funcfile = thisdir + "/OverlayData.txt"
o = np.loadtxt(funcfile)

# now compute the functional values in the atlas
# i.e. make the vertices of each parcel have the correct value
#------------------------------------------------------------------------------
overlay = PyBP.ComputeROIParcels(AALv,AALi,o)

# Now project the source activity onto the mesh faces
#------------------------------------------------------------------------------
#y = PyBP.alignoverlay_icp(v,f,overlay,AALv) 
nv,nf,y = PyBP.alignoverlay_raycast(v,f,overlay,AALv) 

# Finally plot the overlay values on the glass brain
#------------------------------------------------------------------------------
pts,fig = PyBP.PlotMeshwOverlay(v,f,y,.05)


# NETWORK
#------------------------------------------------------------------------------
file     = thisdir + "/exnet.edge"
A        = np.loadtxt(file)  # load .edge file (90x90 matrix)
xx,yy,vv = PyBP.CtoN(A,AALv)


PyBP.PlotNet(xx,yy,vv)
fig.show()


from matplotlib import animation
def rotate(angle):
    ax.view_init(azim=angle)

rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0,362,2),interval=50)
rot_animation.save('NewkRandStatsRotate.gif', dpi=150, writer='imagemagick')

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






