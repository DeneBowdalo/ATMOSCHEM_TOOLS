import glob

def files():

    stem = '/work/home/db876/plotting_tools/AtmosChem_Tools/plotting/grid/*.py'

    f = glob.glob('%s'%(stem))

    f.append('/work/home/db876/plotting_tools/AtmosChem_Tools/spectral_analysis/model_mag_phase_calc_fullgrid_specific.py')
    f.append('/work/home/db876/plotting_tools/AtmosChem_Tools/spectral_analysis/model_fullgrid_specific_process.pbs')
    f.append('/work/home/db876/plotting_tools/AtmosChem_Tools/spectral_analysis/model_mag_phase_calc_fullgrid_spectra.py')
    f.append('/work/home/db876/plotting_tools/AtmosChem_Tools/spectral_analysis/model_spectra_fullgrid_process.pbs')

    return f                                                                                                                  
