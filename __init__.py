#! /usr/bin/env python

'''
__init__

for PyBP

AS2017
'''
from __future__ import absolute_import, division, print_function
from .PyBP import *  
from .nibabel import gifti
#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib as mpl
from matplotlib import cm

from traits.api import HasTraits, Instance, Button, File
from traitsui.api import View, Item, MenuBar, Menu, Action, Separator, \
                         CloseAction, HGroup, VGroup
from traitsui.file_dialog import open_file, TextInfo
import traitsui.toolkit

from pyface.api import FileDialog, OK
from tvtk.pyface.scene_editor import SceneEditor

from mayavi.tools.mlab_scene_model import MlabSceneModel
from mayavi.core.ui.mayavi_scene import MayaviScene
from mayavi.mlab import clf, draw