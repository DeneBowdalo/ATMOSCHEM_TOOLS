import modules
import os
import numpy as np

raw_species = os.getcwd().split('/')[-2]

if raw_species == 'ISOP':
    species = 'C5H8'
elif raw_species == 'NO2-MOLYBDENUM':
    species = 'NO2-M'
elif raw_species == 'NO2-PHOTOLYTIC':
    species = 'NO2-O'
else:
    species = raw_species

#do gridding
modules.grid_data_badc(raw_species,species)
