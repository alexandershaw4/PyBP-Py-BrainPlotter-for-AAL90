# PyBP: (Py) Brain Plotter for AAL90

The idea of this project is two-fold: first, I needed to be able to interpolate data from the 90-node AAL atlas onto subjects cortical meshs. Second, I wanted to be able to plot both a functional overlay and a network (nodes & edges) simultaneously. 

There are two deployment options: 

For mac, download the app (.dmg contains .app) from [here!](https://www.dropbox.com/s/iahvx7m6xtyfzp1/PyBP_G.dmg?dl=0)

Or download this repo and take a look at the example UserScript.py.


# The App.

Download the .dmg and open it. You should see the app inside:

![App Image](app_logo.png)

Open the app, and you're presented with a blank window. Go to File and select 'Open Gifti Surface'. 
Once the brain mesh appears, you can load either an overlay or network.

For a network, select your '.edge' file which is a textfile which contains a 90x90 matrix. 
For an overlay, select a '.txt' file which is a textfile containing a column vector of 90 values.

![GUI_Image](PyBPGUI.png)


# Examples:

Overlay & Network

![both](both.png)

Network

![test net fig](testfig.png)

Overlay

![test overlay fig](fig1.png)

# Rotate (vis Script only at the moment)

![BothRotate](rotation1.gif)
