#!/bin/bash

# Run the python code to make the images with the various options
./region_network_maps.py

# Make merge from one figure to the next
convert -delay 20 region_network_plot_01.png region_network_plot_02.png -morph 20 rnetwork1.gif
convert -delay 20 region_network_plot_02.png region_network_plot_03.png -morph 20 rnetwork2.gif
convert -delay 20 region_network_plot_03.png region_network_plot_04.png -morph 20 rnetwork3.gif
convert -delay 20 region_network_plot_04.png region_network_plot_05.png -morph 20 rnetwork4.gif
convert -delay 20 region_network_plot_05.png region_network_plot_01.png -morph 20 rnetwork5.gif

# Put it all together
convert -loop 0 -delay 18 region_network_plot_01.png region_network_plot_01.png region_network_plot_01.png region_network_plot_01.png region_network_plot_01.png region_network_plot_01.png region_network_plot_01.png region_network_plot_01.png region_network_plot_01.png region_network_plot_01.png rnetwork1.gif region_network_plot_02.png region_network_plot_02.png region_network_plot_02.png region_network_plot_02.png region_network_plot_02.png region_network_plot_02.png region_network_plot_02.png region_network_plot_02.png region_network_plot_02.png region_network_plot_02.png rnetwork2.gif rnetwork3.gif rnetwork4.gif rnetwork5.gif cpdn_network.gif
