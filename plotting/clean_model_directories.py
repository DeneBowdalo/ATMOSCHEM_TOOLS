import os
import glob

all_files = glob.glob('*')
valid_dirs = []

for i in range(len(all_files)):
    if ('obs' not in all_files[i]) & ('ALL' not in all_files[i]) & ('ACCMIP' not in all_files[i]) & ('CCMI' not in all_files[i]) & ('hourly' not in all_files[i]) & ('daily' not in all_files[i]) & ('monthly' not in all_files[i]) & ('clean' not in all_files[i]) & ('set_model' not in all_files[i]) & ('plotting' not in all_files[i]):
        valid_dirs.append(all_files[i])

for d in valid_dirs:
    #get all files in directory
    valid_files = glob.glob('%s/*'%(d))

    for f in valid_files:
        os.remove("%s"%(f))
