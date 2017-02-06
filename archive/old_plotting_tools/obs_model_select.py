import obs_model_lomb_plotter

#List of plotable species:

#O3
#CO
#NO
#NO2
#C2H6
#C3H8
#DMS
#ISOPRENE
#ACETONE
#TEMP
#SURFACE_PRES
#WINDSPEED
#SURFACE_SOLAR_RADIATION
#ABS_HUMIDITY
#REL_HUMIDITY

#Input species below from list above

species = 'O3'


#Gives species exact model tags for convenience
if species == 'ISOPRENE':
    species = 'TRA_6'

elif species == 'ACETONE':
    species = 'ACET'

elif species == 'TEMP':
    species = 'GMAO_TEMP'

elif species == 'SURFACE_PRES':
    species = 'GMAO_PSFC'

elif species == 'WINDSPEED':
    species = 'GMAO_WIND'

elif species == 'SURFACE_SOLAR_RADIATION':
    species = 'GMAO_RADSW'

elif species == 'ABS_HUMIDITY':
    species = 'GMAO_ABSH'

elif species == 'REL_HUMIDITY':
    species = 'GMAO_RHUM'


obs_model_lomb_plotter.plot(species)

