# PyBP: (Py) Brain Plotter for AAL90

The idea of this project is two-fold: first, I needed to be able to interpolate data from the 90-node AAL atlas onto subjects cortical meshs. Second, I wanted to be able to plot both a functional overlay and a network (nodes & edges) simultaneously. 

There are three deployment options: 

1. For mac, download the app (.dmg contains .app) from [here!](https://www.dropbox.com/s/iahvx7m6xtyfzp1/PyBP_G.dmg?dl=0)
2. Download this repo, navigate to it and launch the gui using: $python PyBP_G.py
3. Download this repo, open up your python ide (e.g. spyder) and take a look at the example UserScript.py

4. (Get a template mesh & example overlay and network files, [here](https://www.dropbox.com/sh/w35j02u45602u4g/AACjzoSq-H7uskskiKBois3Ba?dl=0))

# The App.

Download the .dmg and open it. You should see the app inside:

![App Image](app_logo.png)

Open the app, and you're presented with a blank window. Go to File and select 'Open Gifti Surface'. 
Once the brain mesh appears, you can load either an overlay or network.

* For a network, select your '.edge' file which is a textfile which contains a 90x90 matrix. 
* For an overlay, select a '.txt' file which is a textfile containing a column vector of 90 values.

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
