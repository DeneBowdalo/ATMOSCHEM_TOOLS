import os
import glob

present_dir = os.getcwd()

if present_dir.split('/')[-4] == 'xy':
    from xy_plotting_files import files
if present_dir.split('/')[-4] == 'grid':
    from grid_plotting_files import files

#try:
#    from grid_plotting_files import files
#except:
#    pass

#get files to symboliclly link in directories
set_files = files()

all_files = glob.glob('*')
valid_dirs = []

for i in range(len(all_files)):
    if ('obs' not in all_files[i]) & ('ALL' not in all_files[i]) & ('ACCMIP' not in all_files[i]) & ('CCMI' not in all_files[i]) & ('hourly' not in all_files[i]) & ('daily' not in all_files[i]) & ('monthly' not in all_files[i]) & ('clean' not in all_files[i]) & ('plotting' not in all_files[i]) & ('set_model' not in all_files[i]):
        valid_dirs.append(all_files[i])

for d in valid_dirs:
    for f in set_files:
        filename = f.split('/')[-1]
        try:
            os.symlink("%s"%(f),"%s/%s"%(d,filename))
        except:
            continue
