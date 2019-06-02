# wdss2
python scripts to display radar imagery from wdss-ii (also known as wdss2) netcdf files.
http://www.wdssii.org/

Workflow generally goes in the following sequence 

case_data.py
---------------------
Dictionaries of case metadata that get imported into other scripts to inform what/how they should be processing.
This is the best first step since it will ease workflow when subsequent scripts get run.

wdss_create_netcdfs.py
---------------------
Constructs commands to be passed to wdss-ii so it can generate netcdfs from radar files.
Creating netcdfs is usually only a one time deal for each case, unless of course it's later
determined that additional radar data files are desired.

wdss_stage_files.py
---------------------
Moves the desired netcdf files to a staging directory and gives each one a unique name
to eliminate duplicates and any overwriting issues. Case information imported from 'case_data.py'
is what informs the file paths and selections. The new files are sorted in chronological order, which is 
what 'wdss_create_figures.py' wants.

wdss_create_figures.py
---------------------
Takes the staged files and creates multi-pane images. Also leverages 'case_data.py' to control workflow.
Imports custom_maps.py to get the color tables needed to render the maps.

custom_maps.py
---------------------
Color maps that are created and then imported into wdss_create_figures.py

