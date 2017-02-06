import glob

def files():

    stem = '/work/home/db876/plotting_tools/AtmosChem_Tools/plotting/xy/*.py'

    f = glob.glob('%s'%(stem))

    f.append('/work/home/db876/plotting_tools/AtmosChem_Tools/spectral_analysis/obsmodel_mag_phase_calc_specific.py')
    f.append('/work/home/db876/plotting_tools/AtmosChem_Tools/spectral_analysis/obsmodel_specific_process.pbs')
    f.append('/work/home/db876/plotting_tools/AtmosChem_Tools/spectral_analysis/obsmodel_mag_phase_calc_spectra.py')
    f.append('/work/home/db876/plotting_tools/AtmosChem_Tools/spectral_analysis/obsmodel_spectra_process.pbs')

    return f                                                                                                                  
