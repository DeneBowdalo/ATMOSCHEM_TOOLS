import glob
import numpy as np

files = glob.glob('*.npy')

files = sorted(files)

species_2_rm = ['N2O5','PPN','PMN','R4N2','H2O2','MP','CH2O','MO2','ETO2','PRPE','ALK4','ACET','ALD2','MEK','RCHO','MVK','SO2','DMS','SO4','REA_N2O5','REA_307',
                'REA_323','REA_324','REA_325','REA_326','MSA','TRA_1','TRA_2','TRA_3','TRA_4','TRA_5','TRA_6','TRA_7','TRA_8','TRA_9','TRA_10','TRA_11',
                'TRA_12','TRA_13','TRA_14','TRA_15','TRA_16','TRA_17','TRA_18','TRA_19','TRA_20','TRA_21','TRA_22','TRA_23','TRA_24',
                'TRA_25','TRA_26','TRA_27','TRA_28','TRA_29','TRA_30','TRA_31','TRA_32','TRA_33','TRA_34','TRA_35','TRA_36','TRA_37','TRA_38','TRA_39',
                'TRA_40','TRA_41','TRA_42','TRA_43','REA_327','REA_328','REA_329','POINT','TYPE','LAT','LON','PRESS','NO3','HNO4','HNO3','HNO2',
                'HO2','TRA_44','TRA_45','TRA_46','TRA_47','TRA_48','TRA_49','TRA_50','TRA_51','TRA_52','TRA_53','GMAO_ABSH','GMAO_SURF','GMAO_HFLUX']

for f in files:
    valid_list = []
    print f
    data = np.load(f)
    name_list = data.dtype.names
    for spec in name_list:
        if spec not in species_2_rm:
            valid_list.append(name_list.index(spec))
    name_list = np.array(name_list)
    name_list = name_list[valid_list]
    valids = [spec for spec in name_list]
    print valids
    data = data[valids]
    np.save(f,data)
