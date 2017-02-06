import numpy as np
import os
import glob

#read in all files want to convert to binary .npy format
all_files = glob.glob('plane.log.2*')
all_files.sort()

print len(all_files)
if len(all_files) > 0:

    test_input = np.genfromtxt(all_files[0],names=True)
    len_var = len(test_input.dtype.names)

    for one_file in all_files:
        print one_file

    #get timstamp from file processing
        timestamp = one_file[10:]

    #set format and read file
        float_array = (len_var-4)*',f4'
        read = np.genfromtxt(one_file,names=True,dtype='i10,S3,S8,S4%s'%(float_array))

    #save file, using timestamp in name
        np.save('%s'%(timestamp), read)

        os.remove(one_file)
