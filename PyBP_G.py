#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
tkinter based GUI for PyBrainPlot AAL90 tools

AS
"""


import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure

from matplotlib import pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
from matplotlib import cm
import numpy as np
import FileDialog

# get PyBP functions from temp dir for now
#import sys
#sys.path.append("/Users/Alex/code/PyBP/")
from PyBP import *

import os
import sys
if sys.version_info < (3, 0):
    # Python 2
    from Tkinter import *
    import Tkinter as tk
    import Tkinter
else:
    # Python 3
    from tkinter import *
    import tkinter as tk
    import tkinter


import ttk
from tkFileDialog import askopenfile
from tkFileDialog import asksaveasfile


# Here, we are creating our class, Window, and inheriting from the Frame
# class. Frame is a class from the tkinter module. (see Lib/tkinter/__init__)
class Window(Frame):

    # Define settings upon initialization. Here you can specify
    def __init__(self, master=None):
        
        # parameters that you want to send through the Frame class. 
        Frame.__init__(self, master)    

        #reference to the master widget, which is the tk window                 
        self.master = master

        #with that, we want to then run init_window, which doesn't yet exist
        self.init_window()

    #Creation of init_window
    def init_window(self):

        self.v = np.array([0])
        self.f = np.array([0])
        
        # changing the title of our master widget      
        self.master.title("SurfacePlot: AAL90")

        # allowing the widget to take the full space of the root window
        self.pack(fill=BOTH, expand=1)

        
        # creating a menu instance
        menu = Menu(self.master)
        self.master.config(menu=menu)

        # create the file object)
        file = Menu(menu)

        # File menu items
        file.add_command(label="Open Gifti Surface (.gii)", command=self.opengifti)
        file.add_command(label="Inflate Gifti Surface (.gii)", command=self.inflategifti)
        file.add_command(label="Add AAL90 Overlay (.txt)", command=self.addoverlay)
        file.add_command(label="Add AAL90 Network (.edge)", command=self.addnet)
        file.add_command(label="Generate Example Overlay", command=self.genoverlay)
        file.add_command(label="Generate Example Network", command=self.gennet)
        file.add_command(label="Save Image", command=self.saveim)     
        file.add_command(label="Clear", command=self.clearfig)    
        file.add_command(label="Exit", command=self.client_exit)
        
        #added "file" to our menu
        menu.add_cascade(label="File", menu=file)
        
        # if command line included 'overlay' use the root template mesh?
        if len(sys.argv) == 2:
            inpt = sys.argv[1]
            if inpt == "template":
                print("Loading template")
                self.v,self.f = template()
                self.plotgii()
            elif os.path.isfile(inpt):
                pmsg = "Loading mesh: %s\n" % inpt 
                print(pmsg)
                self.v,self.f = GetMesh(inpt)
                self.plotgii()


        
    def client_exit(self):
        exit()
    
    def clearfig(self):
        self.canvas.get_tk_widget().delete("all")
        self.ax = None
        self.fig = None
        self.canvas = None
        self.collec = None
        root.update()        
        return self
        #self.canvas.draw()

                
    def opengifti(self):
        gifti = askopenfile(title="Select file:")
        self.gifti = gifti.name
        str = "\nUsing surface gifti: %s\n" % gifti.name
        print(str)
        self.v,self.f = GetMesh(self.gifti)
        self.plotgii()
        return self
    
    def saveim(self):
        thefile = asksaveasfile(title="Save name: ")
        fname = thefile.name
        fname = fname + ".png"
        self.fig.savefig(fname, format='png', dpi=800)
        
    def inflategifti(self):
        gifti = askopenfile(title="Select file:")
        self.gifti = gifti.name
        str = "\nUsing surface gifti: %s\n" % gifti.name
        print(str)
        self.v,self.f = GetMesh(self.gifti)
        
        self.v = inflate(self.v,self.f)
        self.plotgii()
        return self

    def genoverlay(self):
        overlay = GenTestOverlay()
        self.overlay = overlay
        self.plotoverlay()

    
    def addoverlay(self):
        # if a mesh already opened, add overlay
        
        overlayfile = askopenfile(title="Select overlay (90 value .txt file)")
        overlay = overlayfile.name
    
        strg = "Using overlay file: %s\n" % overlay
        msg  = "Interpolating AAL overlay onto this surface"
        root.update()
        # add a label ?

        print(strg)    
        o = np.loadtxt(overlay)
        self.overlay = o
        self.plotoverlay()
        
    def plotoverlay(self):
        v = self.v
        o = self.overlay
        y = alignoverlay(v,self.f,o)
        
        # Existing fig handle info
        ax = self.ax
        collec = self.collec
        collec.cmap = cm.jet
        
        colors = np.mean(y[self.f], axis=1) # map vertex cols to face cols!
        newy=colors[:,0]
        collec.set_array(newy)
        ax.add_collection(collec)
        self.fig.colorbar(collec, ax=ax)                
        
        canvas = self.canvas
        canvas.show()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)        
        
    
    def gennet(self):
        A = GenTestNet()
        self.A = A
        self.plotnet()


    def addnet(self):
        networkfile = askopenfile(title="Select network file (90x90 value .edge file)")
        network = networkfile.name
        str = "Using edge file: %s\n" % network
        print(str)           
    
        A = np.loadtxt(network)  # load .edge file (90x90 matrix)
        self.A = A
        self.plotnet()

    def plotnet(self):
        A = self.A
        AALv,L   = GetAAL()             # template vertices
        xx,yy,vv = CtoN(A,AALv)         # convert to vertices        
        
        # unpack existing figure info
        canvas = self.canvas
        fig = self.fig
        
        jet = cm.get_cmap('jet') 
        cNorm  = cm.colors.Normalize(vmin=vv.min(), vmax=vv.max())
        scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)
    
    
        for i in range(len(xx)):
            colorVal = scalarMap.to_rgba(vv[i])
            self.ax.plot([xx[i][0], yy[i][0]],[xx[i][1], yy[i][1]],[xx[i][2], yy[i][2]],c=colorVal)
            self.ax.scatter(xx[i][0],xx[i][1],xx[i][2])
            self.ax.scatter(yy[i][0],yy[i][1],yy[i][2])
    
        fig.subplots_adjust(bottom=0.25)
        ax1 = fig.add_axes([0.05, 0.10, 0.9, 0.1])
        norm = mpl.colors.Normalize(vmin=vv.min(), vmax=vv.max())
        cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cm.jet,
                               norm=norm,
                               orientation='horizontal')
        cb1.set_label('Edge strength')                   
        canvas.show()        
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True) 
        
    def orthogonal_proj(self, zfront, zback):
        a = (zfront+zback)/(zfront-zback)
        b = -2*(zfront*zback)/(zfront-zback)
        # -0.0001 added for numerical stability as suggested in:
            # http://stackoverflow.com/questions/23840756
        return numpy.array([[1,0,0,0],
                                [0,1,0,0],
                                [0,0,a,b],
                                [0,0,-0.0001,zback]])
    
    def plotgii(self):
        v = self.v
        # make figure - 
        limits = [self.v.min(), self.v.max()]
        self.fig = plt.figure()
        fig = self.fig
        f = fig
        
        canvas = FigureCanvasTkAgg(f, self)
        self.ax = fig.add_subplot(111, projection='3d',xlim=limits, ylim=limits)
        self.ax.view_init(elev=0, azim=0)
        self.ax.set_axis_off()        
  
        collec = self.ax.plot_trisurf(v[:, 0], v[:, 1], v[:, 2],
                                         triangles=self.f, linewidth=0.,
                                         antialiased=False,
                                         alpha=0.05,
                                         color='white')
        proj3d.persp_transformation = self.orthogonal_proj
        self.collec = collec
        canvas.show()
        canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2TkAgg(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True) 
        self.canvas = canvas
        
        
        
        
        
        
        
        
        
#    def openoverlay(self):
#        gifti = askopenfile(title="Select Surface (.gii):")
#        mesh  = gifti.name
#        str = "\nUsing surface gifti: %s\n" % mesh
#        print(str)
#        
#        overlayfile = askopenfile(title="Select overlay (90 value .txt file)")
#        overlay = overlayfile.name
#        str = "Using overlay file: %s\n" % overlay
#        print(str)        
#        
#        self.v,self.f = GetMesh(mesh)
#        v = self.v
#        o = np.loadtxt(overlay)
#        y = alignoverlay(v,self.f,o)
#        
#        limits = [v.min(), v.max()]
#        #cmap   = 'coolwarm'
#        fig = plt.figure()
#        f = fig
#        canvas = FigureCanvasTkAgg(f, self)
#        
#        self.ax = fig.add_subplot(111, projection='3d', xlim=limits, ylim=limits)
#        self.ax.view_init(elev=0, azim=0)
#        self.ax.set_axis_off()
#    
#        collec = self.ax.plot_trisurf(v[:, 0], v[:, 1], v[:, 2],
#                                 triangles=self.f, linewidth=0.,
#                                 antialiased=False,
#                                 cmap=cm.jet,alpha=.5)
#        #shade=False,
#        #yf = y[f]
#        #colors = np.amax(yf,axis=1)
#    
#        ax = self.ax
#        colors = np.mean(y[self.f], axis=1) # map vertex cols to face cols!
#        newy=colors[:,0]
#        collec.set_array(newy)
#        ax.add_collection(collec)
#        fig.colorbar(collec, ax=ax)
#        
#        canvas.show()
#        canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
#
#        toolbar = NavigationToolbar2TkAgg(canvas, self)
#        toolbar.update()
#        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)         
#        
#        
#    def openall(self):
#        gifti = askopenfile(title="Select Surface (.gii):")
#        mesh  = gifti.name
#        str = "\nUsing surface gifti: %s\n" % mesh
#        print(str)
#        
#        overlayfile = askopenfile(title="Select overlay (90 value .txt file)")
#        overlay = overlayfile.name
#        str = "Using overlay file: %s\n" % overlay
#        print(str)        
#        
#        networkfile = askopenfile(title="Select network file (90x90 value .edge file)")
#        network = networkfile.name
#        str = "Using edge file: %s\n" % network
#        print(str)  
#
#        
#        # MESH
#        self.v,self.f = GetMesh(mesh)
#        v = self.v
#        
#        # OVERLAY
#        o = np.loadtxt(overlay)
#        y = alignoverlay(v,self.f,o)    
#        
#        # NET
#        A        = np.loadtxt(network)  # load .edge file (90x90 matrix)
#        AALv,L   = GetAAL()
#        xx,yy,vv = CtoN(A,AALv)         # convert to vertices
#        
#        # plot mesh and overlay
#        limits = [v.min(), v.max()]
#        #cmap   = 'coolwarm'
#        fig = plt.figure()
#        f = fig
#        canvas = FigureCanvasTkAgg(f, self)
#        
#        self.ax = fig.add_subplot(111, projection='3d', xlim=limits, ylim=limits)
#        self.ax.view_init(elev=0, azim=0)
#        self.ax.set_axis_off()
#    
#        collec = self.ax.plot_trisurf(v[:, 0], v[:, 1], v[:, 2],
#                                 triangles=self.f, linewidth=0.,
#                                 antialiased=False,
#                                 cmap=cm.jet,alpha=.05)
#    
#        ax = self.ax
#        colors = np.mean(y[self.f], axis=1) # map vertex cols to face cols!
#        newy=colors[:,0]
#        collec.set_array(newy)
#        ax.add_collection(collec)
#        fig.colorbar(collec, ax=ax)        
#        
#        # Add net
#        jet = cm.get_cmap('jet') 
#        cNorm  = cm.colors.Normalize(vmin=vv.min(), vmax=vv.max())
#        scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)
#    
#    
#        for i in range(len(xx)):
#            colorVal = scalarMap.to_rgba(vv[i])
#            self.ax.plot([xx[i][0], yy[i][0]],[xx[i][1], yy[i][1]],[xx[i][2], yy[i][2]],c=colorVal)
#            self.ax.scatter(xx[i][0],xx[i][1],xx[i][2])
#            self.ax.scatter(yy[i][0],yy[i][1],yy[i][2])
#    
#        fig.subplots_adjust(bottom=0.25)
#        ax1 = fig.add_axes([0.05, 0.10, 0.9, 0.1])
#        norm = mpl.colors.Normalize(vmin=vv.min(), vmax=vv.max())
#        cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cm.jet,
#                               norm=norm,
#                               orientation='horizontal')
#        cb1.set_label('Edge strength')       
#
#        canvas.show()
#        canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
#
#        toolbar = NavigationToolbar2TkAgg(canvas, self)
#        toolbar.update()
#        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)                  
#        
#        
#        
#    
#    def opennetwork(self):
#        gifti = askopenfile(title="Select Surface (.gii):")
#        mesh  = gifti.name
#        str = "\nUsing surface gifti: %s\n" % mesh
#        print(str)
#        
#        networkfile = askopenfile(title="Select network file (90x90 value .edge file)")
#        network = networkfile.name
#        str = "Using edge file: %s\n" % network
#        print(str)        
#        
#        self.v,self.f = GetMesh(mesh)
#        v = self.v  
#        
#        A        = np.loadtxt(network)  # load .edge file (90x90 matrix)
#        AALv,L   = GetAAL()
#        xx,yy,vv = CtoN(A,AALv)         # convert to vertices
#
#        # make figure - 
#        limits = [self.v.min(), self.v.max()]
#        fig = plt.figure()
#        f = fig
#        # ax = self.ax
#        canvas = FigureCanvasTkAgg(f, self)
#        self.ax = fig.add_subplot(111, projection='3d',xlim=limits, ylim=limits)
#        self.ax.view_init(elev=0, azim=0)
#        self.ax.set_axis_off()        
#  
#        collec = self.ax.plot_trisurf(v[:, 0], v[:, 1], v[:, 2],
#                                         triangles=self.f, linewidth=0.,
#                                         antialiased=False,
#                                         color='white',alpha=0.1)
#        # add the network
#        jet = cm.get_cmap('jet') 
#        cNorm  = cm.colors.Normalize(vmin=vv.min(), vmax=vv.max())
#        scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)
#    
#    
#        for i in range(len(xx)):
#            colorVal = scalarMap.to_rgba(vv[i])
#            self.ax.plot([xx[i][0], yy[i][0]],[xx[i][1], yy[i][1]],[xx[i][2], yy[i][2]],c=colorVal)
#            self.ax.scatter(xx[i][0],xx[i][1],xx[i][2])
#            self.ax.scatter(yy[i][0],yy[i][1],yy[i][2])
#    
#        fig.subplots_adjust(bottom=0.25)
#        ax1 = fig.add_axes([0.05, 0.10, 0.9, 0.1])
#        norm = mpl.colors.Normalize(vmin=vv.min(), vmax=vv.max())
#        cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cm.jet,
#                               norm=norm,
#                               orientation='horizontal')
#        cb1.set_label('Edge strength')       
#            
#        canvas.show()
#        canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
#
#        toolbar = NavigationToolbar2TkAgg(canvas, self)
#        toolbar.update()
#        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)         
#        self.canvas = canvas
#        
        
        
        


        
# root window created. Here, that would be the only window, but
# you can later have windows within windows.
root = Tk()

root.geometry("1100x850")

#creation of an instance
app = Window(root)

#mainloop 
root.mainloop() 