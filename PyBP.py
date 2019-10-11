#! /usr/bin/env python

'''
PyBP: PyBrainPlot for cortical mesh

Includes plotting functions for surface, network and overlays.
Functions for mesh inflation, overlay registration, ...


AS2017
'''
# imports
from __future__ import absolute_import, division, print_function

import os
#os.environ['ETS_TOOLKIT'] = 'qt4' #'qt4'
#os.environ['QT_API'] = 'pyqt'



#import matplotlib
import numpy as np
#import nibabel as nb
from nibabel import gifti
import nibabel as nib
import numpy.matlib
import operator
from matplotlib import cm
import os
from scipy.interpolate import interp1d
from scipy.spatial import Delaunay
from skimage import measure
import sys

#from matplotlib import pyplot as plt
#import matplotlib as mpl
from mayavi import mlab
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d

__all__ = ["GetMesh","PlotMesh","GetAAL","CtoN","PlotNet","cdist","spherefit",
           "minpoints","maxpoints","alignoverlay_icp","PlotMeshwOverlay",
           "GenTestNet","GenTestOverlay","meshadj","normalize_v3","meshnormals",
           "inflate","ReadNiftiOverlay","template","orthogonal_proj",
           "curvature","parseVertices","definesources","ComputeROIParcels",
           "alignoverlay_raycast"]

thisdir = os.path.dirname(os.path.realpath(__file__))

def GetMesh(g):
    # read a gifti (.gii) surface file
    
    #g = '/mesh/cortex_20484.surf.gii'
    v, f = gifti.read(g).getArraysFromIntent(1008)[0].data, \
           gifti.read(g).getArraysFromIntent(1009)[0].data

    Centre = spherefit(v)
    v[0] = v[0] - Centre[0]
    v[1] = v[1] - Centre[1]
    v[2] = v[2] - Centre[2]
    #v = fixboundary(v)

    return v, f

def template():
    g = resource_path('NewSmoothed.gii')
    v,f = GetMesh(g)
    return v,f
    
def meshadj(v,f):
    # Compute adjacency of mesh, given vertices and faces
    N = len(v)
    f = np.int32(f)
    x = np.concatenate((f[:,0],f[:,0],f[:,1],f[:,1],f[:,2],f[:,2]),axis=0)
    y = np.concatenate((f[:,1],f[:,2],f[:,0],f[:,2],f[:,0],f[:,1]),axis=0)
    A = np.empty([N,N])
    
    A[x,y] = 1
    return A
    
def normalize_v3(arr):
    # Normalize a numpy array of 3 component vectors shape=(n,3) 
    lens = np.sqrt( arr[:,0]**2 + arr[:,1]**2 + arr[:,2]**2 )
    arr[:,0] /= lens
    arr[:,1] /= lens
    arr[:,2] /= lens                
    return arr
    
def meshnormals(v,f):
    norm = np.zeros( v.shape, dtype=v.dtype )
    tris = v[f]          
    n = np.cross( tris[::,1 ] - tris[::,0]  , tris[::,2 ] - tris[::,0] )
    n = normalize_v3(n)

    norm[ f[:,0] ] += n
    norm[ f[:,1] ] += n
    norm[ f[:,2] ] += n
    meshnorms = normalize_v3(norm)
    return meshnorms

def parseVertices(vertices):
    # Remove dupicate vertices
    # Perform lex sort and get sorted data
    sorted_idx = np.lexsort(vertices.T)
    sorted_data =  vertices[sorted_idx,:]
    # Get unique row mask
    mask = np.any(np.diff(sorted_data,axis=0),1)
    row_mask = np.append([True],mask)
    # Get unique rows
    vertices = sorted_data[row_mask]
    
    return vertices

def curvature(v,f):
    
    A  = meshadj(v,f)
    C  = np.inner(A,v.T)
    N  = meshnormals(v,f)
    s  = np.sign(sum(N.T*C.T))
    sq = np.sqrt(sum((C*C).T))
    Curv = s * sq
    return Curv
    

def inflate(v,f):
    minxyz = np.array([v[:,0].min(),v[:,1].min(),v[:,2].min()])
    maxxyz = np.array([v[:,0].max(),v[:,1].max(),v[:,2].max()])
    
    T = np.floor(len(v) * 0.003 -2)
    b = 0.5;
    w = 0.05;
    A = meshadj(v,f)
    
    for i in range(int(T)):
        N = meshnormals(v,f)
        mv = np.inner(A,v.T) - v
        
        d = sum((mv*N).T)
        rd = np.matlib.repmat(d,3,1)
        v = v + b * (w*rd.T*N + (1-w)*mv)
        v = np.mean((maxxyz-minxyz)/(v.max() - v.min())) * v
    
    Centre = spherefit(v)
    v[0] = v[0] - Centre[0]
    v[1] = v[1] - Centre[1]
    v[2] = v[2] - Centre[2]
    v = fixboundary(v)
    return v

def orthogonal_proj(zfront, zback):
    a = (zfront+zback)/(zfront-zback)
    b = -2*(zfront*zback)/(zfront-zback)
    # -0.0001 added for numerical stability as suggested in:
    # http://stackoverflow.com/questions/23840756
    return numpy.array([[1,0,0,0],
                        [0,1,0,0],
                        [0,0,a,b],
                        [0,0,-0.0001,zback]])
    
def PlotMesh(v,f):
    # plot a surface (vert & faces) as a 3D patch (trisurf)
    fig = mlab.figure(1, bgcolor=(0, 0, 0))
    pts = mlab.triangular_mesh(v[:,0], v[:,1], v[:,2], f,color=(1,1,1),opacity=0.3)
    mlab.get_engine().scenes[0].scene.x_plus_view()
    mlab.view(0., 0.)
    return pts, fig #ax, p3dcollec, fig

def GetAAL():
    # read the AAL90 source vertex and label lists
    str = "Loading AAL vertices from package directory: %s\n" % thisdir
    print(str)
    pth = resource_path('AALv.txt')
    AALv = np.genfromtxt(pth, delimiter=',')
    #AALv = np.genfromtxt(thisdir + '/AALv.txt', delimiter=',') # source vertex data

    lab = resource_path('aal_labels.txt')
    f = open(lab,'r')
    #f = open(thisdir + '/aal_labels.txt', 'r')
    x = f.readlines()    
    return AALv, x

def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    try:
        # PyInstaller creates a temp folder and stores path in _MEIPASS
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")

    return os.path.join(base_path, relative_path)

def CtoN(A,AALv):
    # convert 90x90 matrix to to-from vertex lists
    xx = np.array([[0,0,0]])
    yy = np.array([[0,0,0]])
    vv = np.array(0)
    for i in range(90):
        a  = A[i]
        to = np.argwhere(a != 0)
        
        for j in range(len(to)):
            xx = np.concatenate( (xx,AALv[to[j]]) , axis=0)
            yy = np.concatenate( (yy,[AALv[i]]) , axis=0)
            vv = np.append(vv,A[i][to[j]])
    xx = xx[1:] # remove empty that was just to initiate array
    yy = yy[1:]
    vv = vv[1:]
    return xx, yy, vv

def PlotNet(xx,yy,vv):#,ax,fig):
    
    
    jet = cm.get_cmap('jet') 
    cNorm  = cm.colors.Normalize(vmin=vv.min(), vmax=vv.max())
    scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)

    for i in range(len(xx)):
        colorVal = scalarMap.to_rgba(vv[i])
        colorVal = colorVal[0:3]
        mlab.plot3d([xx[i][0], yy[i][0]],[xx[i][1], yy[i][1]],[xx[i][2], yy[i][2]],color=colorVal,line_width=10,tube_radius=2)
        mlab.points3d(xx[i][0],xx[i][1],xx[i][2],color=(1,0,0),scale_factor=5)
        mlab.points3d(yy[i][0],yy[i][1],yy[i][2],color=(1,0,0),scale_factor=5)
        
    #mlab.colorbar(title='Edge')
    
    #jet = cm.get_cmap('jet') 
    #cNorm  = cm.colors.Normalize(vmin=vv.min(), vmax=vv.max())
    #scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)
    
    
#    for i in range(len(xx)):
#        colorVal = scalarMap.to_rgba(vv[i])
#        ax.plot([xx[i][0], yy[i][0]],[xx[i][1], yy[i][1]],[xx[i][2], yy[i][2]],c=colorVal)
#        ax.scatter(xx[i][0],xx[i][1],xx[i][2])
#        ax.scatter(yy[i][0],yy[i][1],yy[i][2])
#    
#    fig.subplots_adjust(bottom=0.25)
#    ax1 = fig.add_axes([0.05, 0.10, 0.9, 0.1])
#
#    norm = cm.colors.Normalize(vmin=vv.min(), vmax=vv.max())
#    cb1  = cm.colorbar.ColorbarBase(ax1, cmap=cm.jet,
#                                norm=norm,
#                                orientation='horizontal')
#    cb1.set_label('Edge strength')

    return mlab

#def PlotLabels(xx,yy,ax):
#    for i in range(len(xx)):
#        ax.text(xx[i], y, z, label, zdir)
        
def cdist(mv,v):
    #
    # D = sum( (mv-v).^2 )
    # 
    
    dv = np.matlib.repmat(v,mv.shape[0],1)
    sq = np.square(mv-dv)
    D  = sq.sum(axis=1)
    return D

def spherefit(X):
    # fit sphere to vertex list and find middle
    
    A = np.array([[np.mean(X[:,0]*(X[:,0]-np.mean(X[:,0]))),
                2*np.mean(X[:,0]*(X[:,1]-np.mean(X[:,1]))),
                2*np.mean(X[:,0]*(X[:,2]-np.mean(X[:,2])))],
                [0,
                np.mean(X[:,1]*(X[:,1]-np.mean(X[:,1]))),
                np.mean(X[:,1]*(X[:,2]-np.mean(X[:,2])))],
                [0,0,
                np.mean(X[:,2]*(X[:,2]-np.mean(X[:,2])))]])
    A = A+A.T
    B = np.array([ [np.mean((np.square(X[:,0])+np.square(X[:,1])+np.square(X[:,2]))*(X[:,0]-np.mean(X[:,0])))],
                   [np.mean((np.square(X[:,0])+np.square(X[:,1])+np.square(X[:,2]))*(X[:,1]-np.mean(X[:,1])))],
                   [np.mean((np.square(X[:,0])+np.square(X[:,1])+np.square(X[:,2]))*(X[:,2]-np.mean(X[:,2])))]])
    Centre = np.linalg.lstsq(A,B,rcond=None) # avoid FutureWarning
    Centre = Centre[0]      # final solution is approx matlab unique solution
    return Centre

def fixboundary(v):
    
    Av,l = GetAAL()
    
    m = np.amin(Av,axis=0)
    M = np.amax(Av,axis=0)
    
    m = m * 1.2
    M = M * 1.2
    
    v[:,0] = m[0] + [M[0]-m[0]] * (v[:,0] - v[:,0].min()) / (v[:,0].max() - v[:,0].min() )
    v[:,1] = m[1] + [M[1]-m[1]] * (v[:,1] - v[:,1].min()) / (v[:,1].max() - v[:,1].min() )   
    v[:,2] = m[2] + [M[2]-m[2]] * (v[:,2] - v[:,2].min()) / (v[:,2].max() - v[:,2].min() )    
    
    return v
    
def minpoints(x,n):
    # find minimum n-points (& their indices) from a numpy array
    
    dx = x.copy()
    I  = np.array(range(n))
    V  = np.array(range(n))
    M  = dx.max() + 1
    for i in range(n):
        I[i],V[i] = min(enumerate(dx),key=operator.itemgetter(1))
        dx[I[i]] = M
    return V,I

def maxpoints(x,n):
    # find minimum n-points (& their indices) from a numpy array
    
    dx = x.copy()
    I  = np.array(range(n))
    V  = np.array(range(n))
    M  = dx.min() - 1
    for i in range(n):
        I[i],V[i] = max(enumerate(dx),key=operator.itemgetter(1))
        dx[I[i]] = M
    return V,I
    

def definesources(sourcemodel):
    # this function defines the source locations so that it doesn't have to be AAL data.
    # by inlcuding this, the software is no longer AAL specific but generalises.
    # Sources locations may now be specified from an array, or from a nifti  or
    # gifti file
    if (type(sourcemodel) == str and sourcemodel == 'aal'):
        v,l = GetAAL()
        return v
    elif (type(sourcemodel) == str and sourcemodel[-3:] == 'nii'):
        v,f,vals = ReadNiftiOverlay(sourcemodel)    
        return v, f, vals
    elif (type(sourcemodel) == str and sourcemodel[-3:] == 'gii'):
        v,f = GetMesh(sourcemodel)
        return v
    else:
        v = sourcemodel
        return v

def ComputeROIParcels(v,i,vals):
    # v is a long list of nx3 vertices, and i is nx1, where i[k] identifies which
    # 'parcel' v[k] belongs to. Vals is the same length as unique(i) and is the 
    # corrsponding functional value of that parcel
    new = np.zeros([len(i),1])
    for k in range(len(vals)):
        these = i==k
        new[these] = vals[k]
    return new
    
    

def alignoverlay_icp(mv,f,o,sv):
    # Iterative closest point search to align / project 90-element overlay (o) 
    # matched to AAL source vertices onto a mesh brain - defined  by the 
    # vertices and faces mv & f.
    
    # read AAL sources - not any more!
    #v,l = GetAAL()
    v = sv
    
    # initiate stuff
    OL = np.empty([len(o),len(mv)]) # brain vert * AAL vert matrix
    r  = 1200                       # number of closest points
    w  = np.linspace(.1,1,r)        # 'weights'
    w  = np.flipud(w)
    #M  = np.empty([len(o),len(mv)])
    
    S  = np.array([o.min(),o.max()]) # ensure we can retain min/max vals
    x  = v[:,0]
    
    str = "Performing iterative closest point (ICP) alignment routine...\n"
    print(str)
    
    for i in range(len(x)):
        dist = cdist(mv,v[i])
        [junk,ind] = minpoints(dist,r)
        OL[i,ind]  = w*o[i]
        #M[i,ind]   = w
    
    numnan = 0;
    L = np.zeros([OL.shape[1],1]) # use zeros (instead empty) as speed not issue
    for i in range(OL.shape[1]):
        Col = np.argwhere(OL[:,i] != 0)
        L[i] = sum(OL[:,i] / len(Col))
        if np.isnan(L[i]):
            numnan = numnan + 1
            if numnan == 1:
                print("Instability warning: killing NaNs")
            L[i] = 0

    
    # normalise and rescale
    y = S[0] + [S[1]-S[0]] * (L - L.min()) / (L.max() - L.min() )
    
    return mv,f,y
    
def CentreVerts(v):
    nv = v.shape
    v  = v - np.repeat(spherefit(v).T,nv[0],axis=0)
    dv = v - np.repeat(spherefit(v).T,nv[0],axis=0)
    return dv
    
def alignoverlay_raycast(mv,f,o,sv):
    # Iterative closest point search to align / project 90-element overlay (o) 
    # matched to AAL source vertices onto a mesh brain - defined  by the 
    # vertices and faces mv & f.
    
    str = "Performing ray-casting alignment routine...\n"
    print(str)
    
    # Get upper and lower boundaries of overlay
    S  = np.array([o.min(),o.max()])
    
    # read AAL sources - not any more! use supplied source vertices
    #v,l = GetAAL()
    v = sv
    y = mv[:,1].copy()*0
    
    
    # centre both vertex lists
    mv = CentreVerts(mv)
    v  = CentreVerts(v)
    
    # ensure commmon box boundaries - i.e. we're in the same ballpark
    b = v.min(axis=0)
    B = v.max(axis=0)
    
    mv[:,0] = b[0] + [B[0]-b[0]] * (mv[:,0] - mv[:,0].min()) / (mv[:,0].max() - mv[:,0].min() )
    mv[:,1] = b[1] + [B[1]-b[1]] * (mv[:,1] - mv[:,1].min()) / (mv[:,1].max() - mv[:,1].min() )
    mv[:,2] = b[2] + [B[2]-b[2]] * (mv[:,2] - mv[:,2].min()) / (mv[:,2].max() - mv[:,2].min() )
    
    mv = np.around(mv,decimals=1)
    v  = np.around(v,decimals=1)
    
    # compute vertex normals for source volume
    #from scipy.spatial import Delaunay
    tri = Delaunay(v)
    sf  = tri.simplices.copy()
    sf  = sf[:,0:3]
    sn = meshnormals(v,sf)
    nin = np.isnan(sn)
    sn[nin] = 1
    
    step = np.linspace(-1,1,5)
    
    for i in step:
        points = np.around(v + i*sn,decimals=1)
        
        for k in range(len(points)):
            #hits = points[k,:] == mv
            hits = mv == points[k,:]
            ind  = np.where( sum(hits.T) == 3 )
            y[ind] = y[ind] + o[k]
            
            
    # normalise and rescale
    y = S[0] + [S[1]-S[0]] * (y - y.min()) / (y.max() - y.min() )        
    return mv,f,y
        
        
    
    
#    # initiate stuff
#    OL = np.empty([len(o),len(mv)]) # brain vert * AAL vert matrix
#    r  = 1200                       # number of closest points
#    w  = np.linspace(.1,1,r)        # 'weights'
#    w  = np.flipud(w)
#    #M  = np.empty([len(o),len(mv)])
#    
#    S  = np.array([o.min(),o.max()]) # ensure we can retain min/max vals
#    x  = v[:,0]
#    
#    str = "Performing Raycasting alignment routine...\n"
#    print(str)
#    
#    for i in range(len(x)):
#        dist = cdist(mv,v[i])
#        [junk,ind] = minpoints(dist,r)
#        OL[i,ind]  = w*o[i]
#        #M[i,ind]   = w
#    
#    numnan = 0;
#    L = np.zeros([OL.shape[1],1]) # use zeros (instead empty) as speed not issue
#    for i in range(OL.shape[1]):
#        Col = np.argwhere(OL[:,i] != 0)
#        L[i] = sum(OL[:,i] / len(Col))
#        if np.isnan(L[i]):
#            numnan = numnan + 1
#            if numnan == 1:
#                print("Instability warning: killing NaNs")
#            L[i] = 0

    
    # normalise and rescale
    y = S[0] + [S[1]-S[0]] * (L - L.min()) / (L.max() - L.min() )
    
    return y

def PlotMeshwOverlay(v,f,y,a):
    # plot a surface (vert & faces) as a 3D patch (trisurf) with overlay

    
    fig = mlab.figure(1, bgcolor=(0, 0, 0))
    pts = mlab.triangular_mesh(v[:,0], v[:,1], v[:,2], f,scalars=y[:,0],opacity=a)
    mlab.get_engine().scenes[0].scene.x_plus_view()
    mlab.view(0., 0.)
    mlab.colorbar(title="overlay")
#    limits = [v.min(), v.max()]
#    #cmap   = 'coolwarm'
#    fig = plt.figure()
#    ax = fig.add_subplot(111, projection='3d', xlim=limits, ylim=limits)
#    ax.view_init(elev=0, azim=0)
#    ax.set_axis_off()
#    
#    #cmp = cm.jet(y)
#    #cmp=np.squeeze(cmp)
#            
#    collec = ax.plot_trisurf(v[:, 0], v[:, 1], v[:, 2],
#                                triangles=f, linewidth=0.,
#                                antialiased=False,
#                                cmap=cm.jet,alpha=a)
#    #shade=False,
#    #yf = y[f]
#    #colors = np.amax(yf,axis=1)
#    
#    colors = np.mean(y[f], axis=1) # map vertex cols to face cols!
#    newy=colors[:,0]
#    collec.set_array(newy)
#    ax.add_collection(collec)
#    fig.colorbar(collec, ax=ax)

    return pts, fig #ax, collec, fig
    
    
def GenTestNet():
    A=np.random.randint(2, size=(90, 90))
    A=A*np.random.randint(2, size=(90, 90))
    A=A*np.random.randint(2, size=(90, 90))
    A=A*np.random.randint(2, size=(90, 90))
    A=A*np.random.randint(2, size=(90, 90))
    return A

def GenTestOverlay():
    o = np.random.randint(10,size=(90,1))
    return o
    
def ReadNiftiOverlay(y):
    y = nib.load(y)
    data = np.squeeze(y.get_data())
    #Sv,Sf,N,vals  = measure.marching_cubes(data,0.5)
    Sv,Sf,N,vals  = measure.marching_cubes_lewiner(data,step_size=1)
    #vals = vals[np.newaxis].T
    return Sv,Sf,vals


#    from skimage import restoration
#    from skimage import img_as_float
#    im_float = img_as_float(data)
#    im_denoised = restoration.nl_means_denoising(data, h=0.001)
#    
#    from scipy import ndimage
#    from skimage import morphology
#    # Black tophat transformation (see https://en.wikipedia.org/wiki/Top-hat_transform)
#    hat = ndimage.black_tophat(im_denoised, 7)
#    # Combine with denoised image
##    hat -= 0.3 * im_denoised
##    # Morphological dilation to try to remove some holes in hat image
##    hat = morphology.dilation(hat)
#
#    R = [np.min(hat.flatten()),np.max(hat.flatten())]
#    d = (R[1]-R[0])*0.2
    