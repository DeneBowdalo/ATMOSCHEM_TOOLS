import os
import glob

all_files = glob.glob('*')
valid_files = []

for i in range(len(all_files)):
    try:
        timeres = all_files[i].split('_')[-2]
        if ('obs' not in all_files[i]) & ('ALL' not in all_files[i]) & ('ACCMIP' not in all_files[i]) & ('CCMI' not in all_files[i]) & ('M' == timeres) & ('hourly' not in all_files[i]) & ('daily' not in all_files[i]) & ('monthly' not in all_files[i]) & ('clean' not in all_files[i]) & ('plotting' not in all_files[i]) & ('set_model' not in all_files[i]):
            valid_files.append(all_files[i])
    except:
        continue

#process model lsp specific and bulge_calc
for f in valid_files:
    os.chdir('%s'%(f))
    execfile("model_mag_phase_calc_specific_fullgrid.py")
    os.chdir('../')
