# wdss2
python scripts to display radar imagery from netcdf files. The netcdf files first must be created with a separate 
application called wdss-ii (also known as wdss2):

      http://www.wdssii.org/


Workflow generally goes with the following sequence of scripts...


case_data.py
---------------------
Dictionaries of case metadata that get imported into other scripts to inform what/how/where they should be processing.
This is the best first step since it will ease workflow when subsequent scripts get run.

wdss_create_netcdfs.py
---------------------
Constructs commands that are passed to wdss-ii so it can generate netcdfs from WSR-88D level 2 archive radar files.
Creating netcdfs is usually only a one time deal for each case, unless of course you later decide 
you want to process additional radar data.

Note: importing case metadata from case_data.py doesn't seem to work here (because of using os.system ??).
      Therefore am hard-wiring required arguments into script instead.


wdss_create_figures.py
---------------------
Takes the staged files and uses <a href="https://matplotlib.org/" target="_blank">matplotlib</a> to create multi-pane images. Also leverages 'case_data.py' to control workflow and to define plotting domains.
Imports custom_maps.py to get the color tables needed to render the maps.
Imports several functions from my_functions.py
