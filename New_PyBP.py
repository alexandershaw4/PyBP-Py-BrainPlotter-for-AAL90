#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
TritsUI qt interface for PyBP



AS
"""

# Enthought imports.
#import pyface.qt
import sip
sip.setapi('QDate', 2)
sip.setapi('QDateTime', 2)
sip.setapi('QString', 2)
sip.setapi('QTextStream', 2)
sip.setapi('QTime', 2)
sip.setapi('QUrl', 2)
sip.setapi('QVariant', 2)


from traits.api import HasTraits, Instance, Button, File
from traitsui.api import View, Item, MenuBar, Menu, Action, Separator, \
                         CloseAction, HGroup, VGroup
from traitsui.file_dialog import open_file, TextInfo
import traitsui.toolkit

from pyface.api import FileDialog, OK
from tvtk.pyface.scene_editor import SceneEditor

from mayavi.tools.mlab_scene_model import MlabSceneModel
from mayavi.core.ui.mayavi_scene import MayaviScene

from matplotlib import cm

from mayavi.mlab import clf, draw

import time
import numpy as np

# My actual tools
from PyBP import *



class ActorViewer(HasTraits):

    # The scene model.
    scene = Instance(MlabSceneModel, ())

    
    ######################
    # Using 'scene_class=MayaviScene' adds a Mayavi icon to the toolbar,
    # to pop up a dialog editing the pipeline.
    view = View(Item(name='scene',
                     editor=SceneEditor(scene_class=MayaviScene),
                     show_label=False,
                     resizable=True,
                     height=600,
                     width=1000),
                menubar=MenuBar(
                    Menu(Action(name="Load Gifti", action="opengifti"), # see Controller for
                         Action(name="Inflate Gii", action="inflategii"),
                         Action(name="Template", action="DoTemplate"),
                         Action(name="Load Overlay", action="loadoverlay"), # these callbacks
                         Action(name="Load Network", action="loadnetwork"),
                         Separator(),
                         CloseAction,
                         name="File"),
                         ),
                title="ByBP: AAL90 Brain Plotter",
                resizable=True
                )
    
    

    def __init__(self, **traits):
        HasTraits.__init__(self, **traits)
        #self.DoTemplate()
        
    def DoTemplate(self):
        v,f = template()
        self.DoPlot(v,f)
        
    def DoPlot(self,v,f):
        
        clf()
        self.pts = self.scene.mlab.triangular_mesh(v[:,0], v[:,1], v[:,2], f,color=(1,1,1),opacity=0.3)
        self.scene.mlab.get_engine().scenes[0].scene.x_plus_view()
        self.scene.mlab.draw()
        self.scene.mlab.view(0., 0.)
        self.v = v
        self.f = f
        ActorViewer.v = v
        ActorViewer.f = f
        ActorViewer.plot = self
        return self
        
    def opengifti(self):
        G = GetGifti()       
        G.configure_traits()

        
    def inflategii(self):
        iG = GetGiftiInflate()
        iG.configure_traits()
    
    def loadoverlay(self):
        
        o = LoadOverlay90()
        o.configure_traits()
        
    def alignoverlaydraw(self,o):
        #o = self.o
       
        y   = alignoverlay(self.v,self.f,o)
        v   = self.v # get these from store in ActorViewer
        f   = self.f
        ActorViewer.y = y
        a = 0.3
        
        #fig = mlab.figure(1, bgcolor=(0, 0, 0))
        #ActorViewer.plot.pts.mlab_source.set(x = v[:,0], y = v[:,1], z = v[:,2], triangles=f, scalars=y[:,0],opacity=a)
        ActorViewer.plot.pts = self.scene.mlab.triangular_mesh(v[:,0], v[:,1], v[:,2], f,scalars=y[:,0],opacity=a)
        #pts = self.scene.mlab.triangular_mesh(v[:,0], v[:,1], v[:,2], f,scalars=y[:,0],opacity=a)
        ActorViewer.plot.scene.mlab.get_engine().scenes[0].scene.x_plus_view()
        ActorViewer.plot.scene.mlab.view(0., 0.)
        ActorViewer.plot.scene.mlab.colorbar(title="overlay")
        ActorViewer.plot.scene.mlab.draw()
        
    def loadnetwork(self):
        
        n = LoadNetwork90()
        n.configure_traits()
        
    def PlotNet(self):
        
        xx = self.xx
        yy = self.yy
        vv = self.vv
        
        jet = cm.get_cmap('jet') 
        cNorm  = cm.colors.Normalize(vmin=vv.min(), vmax=vv.max())
        scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)
    
        for i in range(len(xx)):
            colorVal = scalarMap.to_rgba(vv[i])
            colorVal = colorVal[0:3]
            ActorViewer.plot.scene.mlab.plot3d([xx[i][0], yy[i][0]],[xx[i][1], yy[i][1]],[xx[i][2], yy[i][2]],color=colorVal,line_width=10,tube_radius=2)
            ActorViewer.plot.scene.mlab.points3d(xx[i][0],xx[i][1],xx[i][2],color=(1,0,0),scale_factor=5)
            ActorViewer.plot.scene.mlab.points3d(yy[i][0],yy[i][1],yy[i][2],color=(1,0,0),scale_factor=5)
        ActorViewer.plot.scene.mlab.colorbar(title="Network")
        

# Load gifti by popup dialog
#------------------------------------------------------------------
class GetGifti(HasTraits):
    
    Gifti_Name = File
    open = Button('Open...')
    traits_view = View( 
        VGroup( 
            HGroup(
              Item( 'open', show_label = False ),
              '_',
              Item( 'Gifti_Name', style = 'readonly', width = 200 ),
            ),
        ),
        title="Select Gifti",
        )
    def _open_fired(self):
        """ Handles the user clicking the 'Open...' button.
        """
        #inflate = False
        self.docontinue = False
        extns = ['*.gii',]#seems to handle only one extension...
        wildcard='|'.join(extns)

        self.dialog = FileDialog(title='Select gifti surface',
            action='open', wildcard=wildcard,
             default_path = self.Gifti_Name)

    
    def _open_changed(self):
        """ Handles the user clicking the 'Open...' button.
        """
        file_name = open_file()
        if file_name != '':
            self.file_name = file_name
            print("File selected: " + file_name)     
            v,f = GetMesh(file_name)
            self.v = v
            self.f = f
            return ActorViewer().DoPlot(v,f)

        
# Load gifti by popup dialog AND INFLATE![same code]
#------------------------------------------------------------------
class GetGiftiInflate(HasTraits):
        
    Gifti_Inflate_Name = File
    open = Button('Open...')
    traits_view = View( 
        VGroup( 
            HGroup(
              Item( 'open', show_label = False ),
              '_',
              Item( 'Gifti_Inflate_Name', style = 'readonly', width = 200 ),
            ),
        ),
        title="Select Gifti",
        )
    def _open_fired(self):
        """ Handles the user clicking the 'Open...' button.
        """
        #inflate = False
        self.docontinue = False
        extns = ['*.gii',]#seems to handle only one extension...
        wildcard='|'.join(extns)

        self.dialog = FileDialog(title='Select gifti surface',
            action='open', wildcard=wildcard,
             default_path = self.Gifti_Inflate_Name)

    
    def _open_changed(self):
        """ Handles the user clicking the 'Open...' button.
        """
        file_name = open_file()
        if file_name != '':
            self.file_name = file_name
            print("File selected: " + file_name)     
            v,f = GetMesh(file_name)
            v   = inflate(v,f)
            self.v = v
            self.f = f
            return ActorViewer().DoPlot(v,f)

# Load OVERLAY (textfile)
#------------------------------------------------------------------
class LoadOverlay90(HasTraits):
        
    Overlay_Name = File
    open = Button('Open...')
    traits_view = View( 
        VGroup( 
            HGroup(
              Item( 'open', show_label = False ),
              '_',
              Item( 'Overlay_Name', style = 'readonly', width = 200 ),
            ),
        ),
        title="Select Gifti",
        )
    def _open_fired(self):
        """ Handles the user clicking the 'Open...' button.
        """
        #inflate = False
        self.docontinue = False
        extns = ['*.gii',]#seems to handle only one extension...
        wildcard='|'.join(extns)

        self.dialog = FileDialog(title='Select gifti surface',
            action='open', wildcard=wildcard,
             default_path = self.Overlay_Name)

    
    def _open_changed(self):
        """ Handles the user clicking the 'Open...' button.
        """
        file_name = open_file()
        if file_name != '':
            self.file_name = file_name
            print("File selected: " + file_name)     
            o = np.loadtxt(file_name)
            self.o = o
            #y = alignoverlay(self.v,self.f,o) 
            #v = self.v
            #f = self.f
            #self.o = o

            return ActorViewer().alignoverlaydraw(o)
            
# Load NETWORK (node & edge files)
#------------------------------------------------------------------
class LoadNetwork90(HasTraits):
        
    Edge_Name = File
    open = Button('Open...')
    traits_view = View( 
        VGroup( 
            HGroup(
              Item( 'open', show_label = False ),
              '_',
              Item( 'Edge_Name', style = 'readonly', width = 200 ),
            ),
        ),
        title="Select Gifti",
        )
    def _open_fired(self):
        """ Handles the user clicking the 'Open...' button.
        """
        #inflate = False
        self.docontinue = False
        extns = ['*.gii',]#seems to handle only one extension...
        wildcard='|'.join(extns)

        self.dialog = FileDialog(title='Select edge file',
            action='open', wildcard=wildcard,
             default_path = self.Edge_Name)

    
    def _open_changed(self):
        """ Handles the user clicking the 'Open...' button.
        """
        file_name = open_file()
        if file_name != '':
            self.file_name = file_name
            print("File selected: " + file_name)  
            
            AALv,L   = GetAAL()               # AAL sourcemodel
            A        = np.loadtxt(file_name)  # load .edge file (90x90 matrix)
            xx,yy,vv = CtoN(A,AALv)
            ActorViewer.xx = xx # Store
            ActorViewer.yy = yy
            ActorViewer.vv = vv
            ActorViewer.AALv = AALv
            
            return ActorViewer().PlotNet()
            
            
if __name__ == '__main__':
    a = ActorViewer()
    a.configure_traits()
    a.DoTemplate()