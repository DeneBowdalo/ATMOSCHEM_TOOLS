import model_model_lomb_plotter

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

#Input species below from list above

species = 'C2H6'


#Gives species exact model tags for convenience
if species == 'ISOPRENE':
    species = 'TRA_6'

elif species == 'ACETONE':
    species = 'ACET'

elif species == 'TEMP':
    species = 'GMAO_TEMP'

elif species == 'SURFACE_PRES':
    species = 'GMAO_PSFC'


model_model_lomb_plotter.plot(species)

