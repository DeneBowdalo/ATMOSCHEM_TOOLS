import os

basedir = '/work/home/db876/observations/ozonesonde/WOUDC'

for fn in os.listdir(basedir):
    print fn
    full_filename = basedir+'/'+fn
    for fx in os.listdir(full_filename):
        print fx
        current_name = full_filename+'/'+fx
        new_name = full_filename+'/'+fx.upper()
        print new_name
        os.rename(current_name,new_name)
