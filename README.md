# wdss2
python scripts to display radar imagery from wdss-ii (also known as wdss2) netcdf files.
http://www.wdssii.org/

Description of each file

wdss_file_process.py
---------------------
Moves the desired netcdf files to a staging directory and gives each one a unique name
to eliminate duplicates and any overwriting issues. The new filenames are kept in
chronological order to satisfy the next srcipt...

wdss_image_process.py
---------------------
Takes staged files from the wdss_file_process.py script and creates multi-pane images
