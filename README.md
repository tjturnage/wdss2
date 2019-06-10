# wdss2
python scripts to display radar imagery from netcdf files. The netcdf files first must be created with a separate 
application called wdss-ii (also known as wdss2):

      http://www.wdssii.org/


Workflow generally goes with the following sequence of scripts...


case_data.py
---------------------
Dictionaries of case metadata that get imported into other scripts to inform what/how they should be processing.
This is the best first step since it will ease workflow when subsequent scripts get run.

wdss_create_netcdfs.py
---------------------
Constructs commands that are passed to wdss-ii so it can generate netcdfs from WSR-88D level 2 archive radar files.
Creating netcdfs is usually only a one time deal for each case, unless of course you later decide 
you want to process additional radar data.

Note: importing case metadata from case_data.py doesn't seem to work here (because of using os.cmd ??).
      Therefore am hard-wiring required arguments into script instead.

wdss_stage_files.py
---------------------
Moves the desired netcdf files to a staging directory and gives each one a unique name
to eliminate duplicates and any overwriting issues. Case information imported from 'case_data.py'
informs the file paths and selections. The new files are sorted in chronological order, which is 
what 'wdss_create_figures.py' is expecting.

wdss_create_figures.py
---------------------
Takes the staged files and creates multi-pane images. Also leverages 'case_data.py' to control workflow.
Imports custom_maps.py to get the color tables needed to render the maps.

custom_maps.py
---------------------
Color maps that are created and then imported into wdss_create_figures.py
