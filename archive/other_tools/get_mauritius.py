import numpy as np
import modules
import matplotlib.pyplot as plt

model_version = '2x2.5'


species_dict = {'O3':'/work/home/db876/plotting_tools/binary_logs/v90103_2x2.5_GRID_O3.npy','CO':'/work/home/db876/plotting_tools/binary_logs/v90103_2x2.5_GRID_CO.npy','NO':'/work/home/db876/plotting_tools/binary_logs/v90103_2x2.5_GRID_NO.npy','NO2':'/work/home/db876/plotting_tools/binary_logs/v90103_2x2.5_GRID_NO2.npy', \
                'PAN':'/work/home/db876/plotting_tools/binary_logs/v90103_2x2.5_GRID_PAN.npy','C2H6':'/work/home/db876/plotting_tools/binary_logs/v90103_2x2.5_GRID_C2H6.npy','C3H8':'/work/home/db876/plotting_tools/binary_logs/v90103_2x2.5_GRID_C3H8.npy','ISOP':'/work/home/db876/plotting_tools/binary_logs/v90103_2x2.5_GRID_ISOP.npy',}
#species = raw_input('Choose species\n%s\n'%('   '.join(i for i in species_dict)))


for species in species_dict:

    model_f = species_dict[species]

    #get model grid dims. for sim. type
    lat_c,lat_e,lon_c,lon_e = modules.model_grids(model_version)
    gridbox_count = len(lat_c)*len(lon_c)
 
    #get model gridbox for obs site
    gridbox_n = modules.obs_model_gridbox(lat_e,lon_e,-20.2, 57.5)

    model_data  = np.load(model_f)
    model_time = np.arange(0,2191,1./24)

    #get model grid dims. for sim. type
    lat_c,lat_e,lon_c,lon_e = modules.model_grids(model_version)                                                                                                                                                                              
    gridbox_count = len(lat_c)*len(lon_c)
 
    model_var = model_data[gridbox_n::gridbox_count]
    model_var = model_var*1e9

    #set plotting area & background to white
    fig=plt.figure(figsize=(20,12))
    fig.patch.set_facecolor('white')
    ax = fig.add_subplot(1,1,1)


    plt.plot(model_time,model_var)
    plt.xlabel('Days Since 2006',fontsize=16)
    plt.ylabel('Concentration (ppb)',fontsize=16)
    plt.title('GEOS-Chem 2x2.5 %s at Mauritius'%(species),fontsize=20)
    plt.grid(True)
    plt.savefig('mauritius_GC_2x2.5/%s'%(species))
#plt.show()


