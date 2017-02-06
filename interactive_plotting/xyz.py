import matplotlib
matplotlib.use("TkAgg")
import sys
from Tkinter import *
import glob
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
import modules
from matplotlib.widgets import Lasso
from matplotlib import path
import matplotlib.pyplot as plt
from numpy import nonzero
from numpy.random import rand
import collections
import inspect
import colormaps as cmaps

groups = ['standard','spring','summer','autumn','winter','day','night']
param_dict = {'standard':['diurnal_amp','diurnal_ph','seasonal_amp','seasonal_ph','mean','p1','p5','p25','p50','p75','p95','p99','pc_var_daily','pc_var_seasonal','pc_var_full','pc_var_noise','total_var'],'spring':['diurnal_amp','diurnal_ph','mean','p1','p5','p25','p50','p75','p95','p99'],'summer':['diurnal_amp','diurnal_ph','mean','p1','p5','p25','p50','p75','p95','p99'],'autumn':['diurnal_amp','diurnal_ph','mean','p1','p5','p25','p50','p75','p95','p99'],'winter':['diurnal_amp','diurnal_ph','mean','p1','p5','p25','p50','p75','p95','p99'],'day':['seasonal_amp','seasonal_ph','mean','p1','p5','p25','p50','p75','p95','p99'],'night':['seasonal_amp','seasonal_ph','mean','p1','p5','p25','p50','p75','p95','p99']}
params_str = ['diurnal_amp','diurnal_ph','seasonal_amp','seasonal_ph','mean','p1','p5','p25','p50','p75','p95','p99','diurnal_ave_waveform','seasonal_waveform','full_ave_waveform','pc_var_daily','pc_var_seasonal','pc_var_full','pc_var_noise','total_var',
        'diurnal_amp_spring','diurnal_ph_spring','mean_spring','p1_spring','p5_spring','p25_spring','p50_spring','p75_spring','p95_spring','p99_spring','diurnal_waveform_spring',
        'diurnal_amp_summer','diurnal_ph_summer','mean_summer','p1_summer','p5_summer','p25_summer','p50_summer','p75_summer','p95_summer','p99_summer','diurnal_waveform_summer',
        'diurnal_amp_autumn','diurnal_ph_autumn','mean_autumn','p1_autumn','p5_autumn','p25_autumn','p50_autumn','p75_autumn','p95_autumn','p99_autumn','diurnal_waveform_autumn',
        'diurnal_amp_winter','diurnal_ph_winter','mean_winter','p1_winter','p5_winter','p25_winter','p50_winter','p75_winter','p95_winter','p99_winter','diurnal_waveform_winter',
        'seasonal_amp_day','seasonal_ph_day','mean_day','p1_day','p5_day','p25_day','p50_day','p75_day','p95_day','p99_day','seasonal_waveform_day',
        'seasonal_amp_night','seasonal_ph_night','mean_night','p1_night','p5_night','p25_night','p50_night','p75_night','p95_night','p99_night','seasonal_waveform_night',
        'daily_h3_mag','daily_h2_mag','daily_h1_mag','daily_mag','annual_h3_mag','annual_h2_mag','annual_h1_mag','annual_mag']

#Get all valid year ranges and associated species

set_type_desc = []
all_xy = []
all_xy_1 = glob.glob('/work/home/db876/xy/*/*/*/*/LSP_stats.nc')
valid_xy_obs = glob.glob('/work/home/db876/xy/*/*/*/obs_*/LSP_stats.nc')
#remove emission alt sims 
for i in range(len(all_xy_1)):
    if '_NOX' not in all_xy_1[i]:
        all_xy.append(all_xy_1[i])

valid_xy_model = [f for f in all_xy if 'obs_' not in f.split('/')[-2]]
obs_spec = [i.split('/')[5] for i in valid_xy_obs]
obs_range = [i.split('/')[7] for i in valid_xy_obs]
obs_set = [i.split('/')[-2] for i in valid_xy_obs]
set_type_desc=np.append(set_type_desc,['obs']*len(obs_spec))

modelxy_spec = [i.split('/')[5] for i in valid_xy_model]
modelxy_range = [i.split('/')[7] for i in valid_xy_model]
modelxy_set = [i.split('/')[-2] for i in valid_xy_model]
set_type_desc=np.append(set_type_desc,['modelxy']*len(modelxy_spec))

valid_grid = glob.glob('/work/home/db876/grid/*/*/*/*/LSP_stats.nc')
modelgrid_spec = [i.split('/')[5] for i in valid_grid]
modelgrid_range = [i.split('/')[7] for i in valid_grid]
modelgrid_set = [i.split('/')[-2] for i in valid_grid]
set_type_desc=np.append(set_type_desc,['modelgrid']*len(modelgrid_spec))

all_filenames = np.array(valid_xy_obs+valid_xy_model+valid_grid)
all_year_ranges = np.array(obs_range+modelxy_range+modelgrid_range)
all_spec = np.array(obs_spec+modelxy_spec+modelgrid_spec)
all_set = np.array(obs_set+modelxy_set+modelgrid_set)
unique_year_ranges = np.unique(all_year_ranges)

year_range_dict = {}
for yr in unique_year_ranges:
    spec_cut = np.unique(all_spec[all_year_ranges == yr])
    year_range_dict[yr] = list(spec_cut) 
 
#------------------------------
#all_year_ranges = np.array(['2005_2010','2005_2010','2005_2010','2005_2010','2006_2012','2006_2012','2006_2012','2006_2012','2006_2012','2006_2012'])
#all_spec = np.array(['O3','O3','O3','NO','NO','NO','NO','CO','NO2','NO2'])
#all_set = np.array(['obs_H','GEOS-Chem','GEOS-Chem','obs_H','GEOS-Chem','GEOS-Chem-FP','GFDLAM3','obs_H','obs_H','GFDLAM3'])
#set_type_desc = np.array(['obs','modelxy','modelgrid','obs','modelxy','modelxy','modelgrid','obs','obs','modelgrid'])
#unique_year_ranges = np.array(['2005_2010','2006_2012'])
#year_range_dict = {'2005_2010': ['O3','NO'],'2006_2012': ['NO','CO','NO2']}
    
#-------------------------------------------------------------------
#SETUP TKINTER MENU SYSTEM
class App(Frame):

    def __init__(self, master):
        Frame.__init__(self, master)
        
        #setup initial default menu strings
        self.yearrange_str = StringVar(self)
        self.speciesx_str = StringVar(self)
        self.speciesy_str = StringVar(self)
        self.speciesz_str = StringVar(self)
        self.typex_str = StringVar(self)
        self.typey_str = StringVar(self)
        self.typez_str = StringVar(self)
        self.setx_str = StringVar(self)
        self.sety_str = StringVar(self)
        self.setz_str = StringVar(self)
        self.groupx_str = StringVar(self)
        self.groupy_str = StringVar(self)
        self.groupz_str = StringVar(self)
        self.paramx_str = StringVar(self)
        self.paramy_str = StringVar(self)
        self.paramz_str = StringVar(self)

        #setup updates to menu on changes
        self.yearrange_str.trace('w', self.change_year)
        self.speciesx_str.trace('w', self.change_speciesx)
        self.speciesy_str.trace('w', self.change_speciesy)
        self.speciesz_str.trace('w', self.change_speciesz)
        self.typex_str.trace('w', self.change_typex)
        self.typey_str.trace('w', self.change_typey)
        self.typez_str.trace('w', self.change_typez)
        self.groupx_str.trace('w', self.change_groupx)
        self.groupy_str.trace('w', self.change_groupy)
        self.groupz_str.trace('w', self.change_groupz)

        #setup menus with labels
        self.optionmenu_yearrange = OptionMenu(self, self.yearrange_str, *year_range_dict.keys())
        self.optionmenu_yearrange_label = Label(self)
        self.optionmenu_yearrange_label['text'] = 'Year Range'
        self.optionmenu_speciesx = OptionMenu(self, self.speciesx_str, '')
        self.optionmenu_speciesx_label = Label(self)
        self.optionmenu_speciesx_label['text'] = 'X'
        self.optionmenu_speciesy = OptionMenu(self, self.speciesy_str, '')
        self.optionmenu_speciesy_label = Label(self)
        self.optionmenu_speciesy_label['text'] = 'Y'
        self.optionmenu_speciesz = OptionMenu(self, self.speciesz_str, '')
        self.optionmenu_speciesz_label = Label(self)
        self.optionmenu_speciesz_label['text'] = 'Z'
        self.optionmenu_typex = OptionMenu(self, self.typex_str, '')
        self.optionmenu_typex_label = Label(self)
        self.optionmenu_typex_label['text'] = 'X'
        self.optionmenu_setx = OptionMenu(self, self.setx_str, '')
        self.optionmenu_groupx = OptionMenu(self, self.groupx_str, *groups)
        self.optionmenu_paramx = OptionMenu(self, self.paramx_str, '')
        self.optionmenu_typey = OptionMenu(self, self.typey_str, '')
        self.optionmenu_typey_label = Label(self)
        self.optionmenu_typey_label['text'] = 'Y'
        self.optionmenu_sety = OptionMenu(self, self.sety_str, '')
        self.optionmenu_groupy = OptionMenu(self, self.groupy_str, *groups)
        self.optionmenu_paramy = OptionMenu(self, self.paramy_str, '')
        self.optionmenu_typez = OptionMenu(self, self.typez_str, '')
        self.optionmenu_typez_label = Label(self)
        self.optionmenu_typez_label['text'] = 'Z'
        self.optionmenu_setz = OptionMenu(self, self.setz_str, '')
        self.optionmenu_groupz = OptionMenu(self, self.groupz_str, *groups)
        self.optionmenu_paramz = OptionMenu(self, self.paramz_str, '')
    
        #set initial year range
        self.yearrange_str.set(unique_year_ranges[0])
    
        #create button to update plot
        self.button = Button(self,text="Update", command=self.update)
    
        #setup layout of menus
        self.optionmenu_yearrange_label.grid(row=0,column=0)
        self.optionmenu_yearrange.grid(row=0,column=1)
        self.optionmenu_speciesx_label.grid(row=0,column=2)
        self.optionmenu_speciesx.grid(row=0,column=3)
        self.optionmenu_speciesy_label.grid(row=0,column=4)
        self.optionmenu_speciesy.grid(row=0,column=5)
        self.optionmenu_speciesz_label.grid(row=0,column=6)
        self.optionmenu_speciesz.grid(row=0,column=7)
        self.optionmenu_typex_label.grid(row=0,column=8)
        self.optionmenu_typex.grid(row=0,column=9)
        self.optionmenu_setx.grid(row=0,column=10)
        self.optionmenu_groupx.grid(row=0,column=11)
        self.optionmenu_paramx.grid(row=0,column=12)
        self.optionmenu_typey_label.grid(row=1,column=8)
        self.button.grid(row=1,column=1)
        self.optionmenu_typey.grid(row=1,column=9)
        self.optionmenu_sety.grid(row=1,column=10)
        self.optionmenu_groupy.grid(row=1,column=11)
        self.optionmenu_paramy.grid(row=1,column=12)
        self.optionmenu_typez_label.grid(row=2,column=8)
        self.optionmenu_typez.grid(row=2,column=9)
        self.optionmenu_setz.grid(row=2,column=10)
        self.optionmenu_groupz.grid(row=2,column=11)
        self.optionmenu_paramz.grid(row=2,column=12)
        self.pack()

    def change_year(self, *args):
        print 'update year'
        
        self.update_year_flag = True
        
        #get valid species for year range and set up species lists for x,y,z variables
        species = year_range_dict[self.yearrange_str.get()]
        self.speciesx_str.set(species[0])
        self.speciesy_str.set(species[0])
        self.speciesz_str.set(species[0])

        menux = self.optionmenu_speciesx['menu']
        menuy = self.optionmenu_speciesy['menu']
        menuz = self.optionmenu_speciesz['menu']
        menux.delete(0, 'end')
        menuy.delete(0, 'end')
        menuz.delete(0, 'end')

        for spec in species:
            menux.add_command(label=spec, command=lambda s=spec: self.speciesx_str.set(s))
            menuy.add_command(label=spec, command=lambda s=spec: self.speciesy_str.set(s))
            menuz.add_command(label=spec, command=lambda s=spec: self.speciesz_str.set(s))
            
        self.update_year_flag = False
        
    def change_speciesx(self, *args):
        print 'update species x'
        
        self.update_species_flag = True
        
        current_yearrange = self.yearrange_str.get()
        
        #get current x_types
        x_types,ind = np.unique(set_type_desc[(all_year_ranges == current_yearrange) & (all_spec == self.speciesx_str.get())],return_index=True)
        x_types = list(x_types[np.argsort(ind)])
        
        #if updating year range all params are same, so set same as x
        if self.update_year_flag == True:
            current_x_species = self.speciesx_str.get()
            current_x_type = x_types[0]
            current_y_species = current_x_species
            current_y_type = x_types[0]
            current_z_species = current_x_species
            current_z_type = x_types[0]
        #else if updating species get params individually
        else:
            current_x_species = self.speciesx_str.get()
            current_x_type = self.typex_str.get()
            current_y_species = self.speciesy_str.get()
            current_y_type = self.typey_str.get()
            current_z_species = self.speciesz_str.get()
            current_z_type = self.typez_str.get()
        
        print current_x_species,current_y_species,current_z_species
        print current_x_type,current_y_type,current_z_type
        
        #update types

        #deal with rules in updating types when changing species
        #3 same species
        if (current_x_species == current_y_species) & (current_x_species == current_z_species):
            print '3 same species'
            self.species_flag = '3 same'
            #Cannot be 3 modelgrids but all other options are fine
            if (current_y_type != 'modelgrid') or (current_z_type != 'modelgrid'):
                #remove modelgrid type
                try:
                    x_types.remove('modelgrid')
                except:
                    pass
        #3 different species
        elif (current_x_species != current_y_species) & (current_x_species != current_z_species) & (current_y_species != current_z_species):
            print '3 diff species'
            self.species_flag = '3 diff'
            #Must be only 1 obs or modelxy type
            #if y or z are obs or modelxy, x must be grid
            if (current_y_type != 'modelgrid') or (current_z_type != 'modelgrid'):
                #if modelgrid in x_types, otherwise set blank
                if 'modelgrid' in x_types:
                    x_types = ['modelgrid']
                else:
                    x_types = ['']
            #if y & z both modelgrid then x must be obs or modelxy
            elif all([current_y_type == 'modelgrid', current_z_type == 'modelgrid']):
                #remove modelgrid type
                try:
                    x_types.remove('modelgrid')
                except:
                    pass
        #1,2 species
        else:
            print '1,2 species'
            self.species_flag = '1,2'
            #Must be only 1 obs or modelxy type
            #if y or z are obs or modelxy, x must be grid
            if (current_y_type != 'modelgrid') or (current_z_type != 'modelgrid'):
                #if modelgrid in x_types, otherwise set blank
                if 'modelgrid' in x_types:
                    x_types = ['modelgrid']
                else:
                    x_types = ['']
        
            #if y or z both modelgrid then x must be obs or modelxy
            elif all([current_y_type == 'modelgrid', current_z_type == 'modelgrid']):
                #remove modelgrid type
                try:
                    x_types.remove('modelgrid')
                except:
                    pass
        
        #set x_types
        if len(x_types) == 0:
            x_types = [''] 
        self.typex_str.set(x_types[0])
        menu = self.optionmenu_typex['menu']
        menu.delete(0, 'end')
        for x_type in x_types:
            menu.add_command(label=x_type, command=lambda s=x_type: self.typex_str.set(s))
           
        #refresh groups
        self.groupx_str.set(groups[0])
        
        self.update_species_flag = False
        
    def change_speciesy(self, *args):
        print 'update species y'

        self.update_species_flag = True

        current_yearrange = self.yearrange_str.get()
        
        #get current y_types
        y_types,ind = np.unique(set_type_desc[(all_year_ranges == current_yearrange) & (all_spec == self.speciesy_str.get())],return_index=True)
        y_types = list(y_types[np.argsort(ind)])
        
        #if updating year range all params are same, so set same as x
        if self.update_year_flag == True:
            current_x_species = self.speciesx_str.get()
            current_x_type = y_types[0]
            current_y_species = current_x_species
            current_y_type = y_types[0]
            current_z_species = current_x_species
            current_z_type = y_types[0]
        #else if updating species get params individually
        else:
            current_x_species = self.speciesx_str.get()
            current_x_type = self.typex_str.get()
            current_y_species = self.speciesy_str.get()
            current_y_type = self.typey_str.get()
            current_z_species = self.speciesz_str.get()
            current_z_type = self.typez_str.get()
        
        #update types
        
        #deal with rules in updating types when changing species
        #3 same species
        if (current_x_species == current_y_species) & (current_x_species == current_z_species):
            print '3 same species'
            self.species_flag = '3 same'
            #Cannot be 3 modelgrids but all other options are fine
            if (current_x_type != 'modelgrid') or (current_z_type != 'modelgrid'):
                #remove modelgrid type
                try:
                    y_types.remove('modelgrid')
                except:
                    pass
        #3 different species
        elif (current_x_species != current_y_species) & (current_x_species != current_z_species) & (current_y_species != current_z_species):
            print '3 diff species'
            self.species_flag = '3 diff'
            #Must be only 1 obs or modelxy type
            #if x or z are obs or modelxy, y must be grid
            if (current_x_type != 'modelgrid') or (current_z_type != 'modelgrid'):
                #if modelgrid in y_types, otherwise set blank
                if 'modelgrid' in y_types:
                    y_types = ['modelgrid']
                else:
                    y_types = ['']
            #if x & z both modelgrid then y must be obs or modelxy
            elif all([current_x_type == 'modelgrid', current_z_type == 'modelgrid']):
                #remove modelgrid type
                try:
                    y_types.remove('modelgrid')
                except:
                    pass
        #1,2 species
        else:
            print '1,2 species'
            self.species_flag = '1,2'
            #Must be only 1 obs or modelxy type
            #if x or z are obs or modelxy, y must be grid
            if (current_x_type != 'modelgrid') or (current_z_type != 'modelgrid'):
                #if modelgrid in y_types, otherwise set blank
                if 'modelgrid' in y_types:
                    y_types = ['modelgrid']
                else:
                    y_types = ['']
        
            #if x or z both modelgrid then y must be obs or modelxy
            elif all([current_x_type == 'modelgrid', current_z_type == 'modelgrid']):
                #remove modelgrid type
                try:
                    y_types.remove('modelgrid')
                except:
                    pass
        
        #set y_types
        if len(y_types) == 0:
            y_types = [''] 
        self.typey_str.set(y_types[0])
        menu = self.optionmenu_typey['menu']
        menu.delete(0, 'end')
        for y_type in y_types:
            menu.add_command(label=y_type, command=lambda s=y_type: self.typey_str.set(s))
           
        #refresh groups
        self.groupy_str.set(groups[0])
        
        self.update_species_flag = False
            
    def change_speciesz(self, *args):
        print 'update species z'

        self.update_species_flag = True

        current_yearrange = self.yearrange_str.get()
        
        #get current z_types
        z_types,ind = np.unique(set_type_desc[(all_year_ranges == current_yearrange) & (all_spec == self.speciesz_str.get())],return_index=True)
        z_types = list(z_types[np.argsort(ind)])
        
        #if updating year range all params are same, so set same as x
        if self.update_year_flag == True:
            current_x_species = self.speciesx_str.get()
            current_x_type = z_types[0]
            current_y_species = current_x_species
            current_y_type = z_types[0]
            current_z_species = current_x_species
            current_z_type = z_types[0]
        #else if updating species get params individually
        else:
            current_x_species = self.speciesx_str.get()
            current_x_type = self.typex_str.get()
            current_y_species = self.speciesy_str.get()
            current_y_type = self.typey_str.get()
            current_z_species = self.speciesz_str.get()
            current_z_type = self.typez_str.get()
        
        #update types
        
        #deal with rules in updating types when changing species
        #3 same species
        if (current_x_species == current_y_species) & (current_x_species == current_z_species):
            print '3 same species'
            self.species_flag = '3 same'
            #Cannot be 3 modelgrids but all other options are fine
            if (current_x_type != 'modelgrid') or (current_y_type != 'modelgrid'):
                #remove modelgrid type
                try:
                    z_types.remove('modelgrid')
                except:
                    pass
        #3 different species
        elif (current_x_species != current_y_species) & (current_x_species != current_z_species) & (current_y_species != current_z_species):
            print '3 diff species'
            self.species_flag = '3 diff'
            #Must be only 1 obs or modelxy type
            #if x or y are obs or modelxy, z must be grid
            if (current_x_type != 'modelgrid') or (current_y_type != 'modelgrid'):
                #if modelgrid in z_types, otherwise set blank
                if 'modelgrid' in z_types:
                    z_types = ['modelgrid']
                else:
                    z_types = ['']
            #if x & y both modelgrid then z must be obs or modelxy
            elif all([current_x_type == 'modelgrid', current_y_type == 'modelgrid']):
                #remove modelgrid type
                try:
                    z_types.remove('modelgrid')
                except:
                    pass
        #1,2 species
        else:
            print '1,2 species'
            self.species_flag = '1,2'
            #Must be only 1 obs or modelxy type
            #if x or y are obs or modelxy, z must be grid
            if (current_x_type != 'modelgrid') or (current_y_type != 'modelgrid'):
                #if modelgrid in z_types, otherwise set blank
                if 'modelgrid' in z_types:
                    z_types = ['modelgrid']
                else:
                    z_types = ['']
        
            #if x or y both modelgrid then z must be obs or modelxy
            elif all([current_x_type == 'modelgrid', current_y_type == 'modelgrid']):
                #remove modelgrid type
                try:
                    z_types.remove('modelgrid')
                except:
                    pass
        
        #set y_types
        if len(z_types) == 0:
            z_types = [''] 
        self.typez_str.set(z_types[0])
        menu = self.optionmenu_typez['menu']
        menu.delete(0, 'end')
        for z_type in z_types:
            menu.add_command(label=z_type, command=lambda s=z_type: self.typez_str.set(s))
           
        #refresh groups
        self.groupz_str.set(groups[0])
        
        self.update_species_flag = False
    
    def change_typex(self, *args):
        print 'update type x'
        
        current_yearrange = self.yearrange_str.get()
        current_x_species = self.speciesx_str.get()
        current_x_type = self.typex_str.get()
        
        #update x sets 
        if current_x_type == '':
            sets = ['']
        else:
            sets = list(all_set[(all_year_ranges == current_yearrange) & (all_spec == current_x_species) & (set_type_desc == current_x_type)])        
        self.setx_str.set(sets[0])
        menu = self.optionmenu_setx['menu']
        menu.delete(0, 'end')
        for set in sets:
            menu.add_command(label=set, command=lambda s=set: self.setx_str.set(s))
        
        #update y and z types when changing x type
        current_y_species = self.speciesy_str.get()
        current_y_type = self.typey_str.get()
        current_y_types = self.optionmenu_typey.keys()
        current_z_species = self.speciesz_str.get()
        current_z_type = self.typez_str.get()
        current_z_types = self.optionmenu_typez.keys()
        y_types,ind = np.unique(set_type_desc[(all_year_ranges == current_yearrange) & (all_spec == self.speciesy_str.get())],return_index=True)
        y_types = list(y_types[np.argsort(ind)])
        z_types,ind = np.unique(set_type_desc[(all_year_ranges == current_yearrange) & (all_spec == self.speciesz_str.get())],return_index=True)
        z_types = list(z_types[np.argsort(ind)])
        
        
        #update y and z types - if possible
        #only update if called from species update
        if (self.update_species_flag == True) & (self.update_year_flag == False):
            #update y types - if can
            if y_types != current_y_types:
                reset_flag = False
                #3 same species
                if self.species_flag == '3 same':
                     #Cannot be 3 modelgrids but all other options are fine
                     if (current_x_type != 'modelgrid') or (current_z_type != 'modelgrid'):
                        reset_flag = True
                #3 diff species
                elif self.species_flag == '3 diff':
                     #Must be only 1 obs or modelxy type, #if x or z are obs or modelxy, y must be grid
                     if all([current_x_type != 'obs', current_x_type == 'modelxy']) & all([current_z_type != 'obs', current_z_type == 'modelxy']):
                        reset_flag = True  
                #1,2 species
                else:
                     #Must be only 1 obs or modelxy type, #if x or z are obs or modelxy, y must be grid
                     if all([current_x_type != 'obs', current_x_type == 'modelxy']) & all([current_z_type != 'obs', current_z_type == 'modelxy']):
                        reset_flag = True  
                    
                #reset y types if can
                if reset_flag == True:
                    #reset y types
                     self.typey_str.set(y_types[0])
                     menu = self.optionmenu_typey['menu']
                     menu.delete(0, 'end')
                     for type in y_types:
                         menu.add_command(label=type, command=lambda s=type: self.typey_str.set(s))
                     #update y sets
                     sets = list(all_set[(all_year_ranges == current_yearrange) & (all_spec == current_y_species) & (set_type_desc == y_types[0])])        
                     self.sety_str.set(sets[0])
                     menu = self.optionmenu_sety['menu']
                     menu.delete(0, 'end')
                     for set in sets:
                        menu.add_command(label=set, command=lambda s=set: self.sety_str.set(s))
                        
            #update z types - if can
            if z_types != current_z_types:
                reset_flag = False
                #3 same species
                if self.species_flag == '3 same':
                     #Cannot be 3 modelgrids but all other options are fine
                     if (current_x_type != 'modelgrid') or (current_y_type != 'modelgrid'):
                        reset_flag = True
                #3 diff species
                elif self.species_flag == '3 diff':
                     #Must be only 1 obs or modelxy type, #if x or y are obs or modelxy, z must be grid
                     if all([current_x_type != 'obs', current_x_type == 'modelxy']) & all([current_y_type != 'obs', current_y_type == 'modelxy']):
                        reset_flag = True  
                #1,2 species
                else:
                     #Must be only 1 obs or modelxy type, #if x or y are obs or modelxy, z must be grid
                     if all([current_x_type != 'obs', current_x_type == 'modelxy']) & all([current_y_type != 'obs', current_y_type == 'modelxy']):
                        reset_flag = True  
                    
                #reset z types if can
                if reset_flag == True:
                    #reset y types
                     self.typez_str.set(z_types[0])
                     menu = self.optionmenu_typez['menu']
                     menu.delete(0, 'end')
                     for type in z_types:
                         menu.add_command(label=type, command=lambda s=type: self.typez_str.set(s))
                     #update z sets
                     sets = list(all_set[(all_year_ranges == current_yearrange) & (all_spec == current_z_species) & (set_type_desc == z_types[0])])        
                     self.setz_str.set(sets[0])
                     menu = self.optionmenu_setz['menu']
                     menu.delete(0, 'end')
                     for set in sets:
                        menu.add_command(label=set, command=lambda s=set: self.setz_str.set(s))
        

            
    def change_typey(self, *args):
        print 'update type y'
        
        current_yearrange = self.yearrange_str.get()
        current_y_species = self.speciesy_str.get()
        current_y_type = self.typey_str.get()
        
        #update y sets 
        if current_y_type == '':
            sets = ['']
        else:
            sets = list(all_set[(all_year_ranges == current_yearrange) & (all_spec == current_y_species) & (set_type_desc == current_y_type)])        
        self.sety_str.set(sets[0])
        menu = self.optionmenu_sety['menu']
        menu.delete(0, 'end')
        for set in sets:
            menu.add_command(label=set, command=lambda s=set: self.sety_str.set(s))
        
        #update x and z types when changing y type
        current_x_species = self.speciesx_str.get()
        current_x_type = self.typex_str.get()
        current_x_types = self.optionmenu_typex.keys()
        current_z_species = self.speciesz_str.get()
        current_z_type = self.typez_str.get()
        current_z_types = self.optionmenu_typez.keys()
        x_types,ind = np.unique(set_type_desc[(all_year_ranges == current_yearrange) & (all_spec == self.speciesx_str.get())],return_index=True)
        x_types = list(x_types[np.argsort(ind)])
        z_types,ind = np.unique(set_type_desc[(all_year_ranges == current_yearrange) & (all_spec == self.speciesz_str.get())],return_index=True)
        z_types = list(z_types[np.argsort(ind)])
        
        
        #update x and z types - if possible
        #only update if called from species update
        if (self.update_species_flag == True) & (self.update_year_flag == False):
            #update x types - if can
            if x_types != current_x_types:
                reset_flag = False
                #3 same species
                if self.species_flag == '3 same':
                     #Cannot be 3 modelgrids but all other options are fine
                     if (current_y_type != 'modelgrid') or (current_z_type != 'modelgrid'):
                        reset_flag = True
                #3 diff species
                elif self.species_flag == '3 diff':
                     #Must be only 1 obs or modelxy type, #if y or z are obs or modelxy, x must be grid
                     if all([current_y_type != 'obs', current_y_type == 'modelxy']) & all([current_z_type != 'obs', current_z_type == 'modelxy']):
                        reset_flag = True  
                #1,2 species
                else:
                     #Must be only 1 obs or modelxy type, #if y or z are obs or modelxy, x must be grid
                     if all([current_y_type != 'obs', current_y_type == 'modelxy']) & all([current_z_type != 'obs', current_z_type == 'modelxy']):
                        reset_flag = True  
                    
                #reset x types if can
                if reset_flag == True:
                    #reset y types
                     self.typex_str.set(x_types[0])
                     menu = self.optionmenu_typex['menu']
                     menu.delete(0, 'end')
                     for type in x_types:
                         menu.add_command(label=type, command=lambda s=type: self.typex_str.set(s))
                     #update x sets
                     sets = list(all_set[(all_year_ranges == current_yearrange) & (all_spec == current_x_species) & (set_type_desc == x_types[0])])        
                     self.setx_str.set(sets[0])
                     menu = self.optionmenu_setx['menu']
                     menu.delete(0, 'end')
                     for set in sets:
                        menu.add_command(label=set, command=lambda s=set: self.setx_str.set(s))
                        
            #update z types - if can
            if z_types != current_z_types:
                reset_flag = False
                #3 same species
                if self.species_flag == '3 same':
                     #Cannot be 3 modelgrids but all other options are fine
                     if (current_x_type != 'modelgrid') or (current_y_type != 'modelgrid'):
                        reset_flag = True
                #3 diff species
                elif self.species_flag == '3 diff':
                     #Must be only 1 obs or modelxy type, #if x or y are obs or modelxy, z must be grid
                     if all([current_x_type != 'obs', current_x_type == 'modelxy']) & all([current_y_type != 'obs', current_y_type == 'modelxy']):
                        reset_flag = True  
                #1,2 species
                else:
                     #Must be only 1 obs or modelxy type, #if x or y are obs or modelxy, z must be grid
                     if all([current_x_type != 'obs', current_x_type == 'modelxy']) & all([current_y_type != 'obs', current_y_type == 'modelxy']):
                        reset_flag = True  
                    
                #reset z types if can
                if reset_flag == True:
                    #reset y types
                     self.typez_str.set(z_types[0])
                     menu = self.optionmenu_typez['menu']
                     menu.delete(0, 'end')
                     for type in z_types:
                         menu.add_command(label=type, command=lambda s=type: self.typez_str.set(s))
                     #update z sets
                     sets = list(all_set[(all_year_ranges == current_yearrange) & (all_spec == current_z_species) & (set_type_desc == z_types[0])])        
                     self.setz_str.set(sets[0])
                     menu = self.optionmenu_setz['menu']
                     menu.delete(0, 'end')
                     for set in sets:
                        menu.add_command(label=set, command=lambda s=set: self.setz_str.set(s))

            
    def change_typez(self, *args):
        print 'update type z'
        
        current_yearrange = self.yearrange_str.get()
        current_z_species = self.speciesz_str.get()
        current_z_type = self.typez_str.get()
        
        #update z sets 
        if current_z_type == '':
            sets = ['']
        else:
            sets = list(all_set[(all_year_ranges == current_yearrange) & (all_spec == current_z_species) & (set_type_desc == current_z_type)])        
        self.setz_str.set(sets[0])
        menu = self.optionmenu_setz['menu']
        menu.delete(0, 'end')
        for set in sets:
            menu.add_command(label=set, command=lambda s=set: self.setz_str.set(s))
        
        #update x and y types when changing z type
        current_x_species = self.speciesx_str.get()
        current_x_type = self.typex_str.get()
        current_x_types = self.optionmenu_typex.keys()
        current_y_species = self.speciesy_str.get()
        current_y_type = self.typey_str.get()
        current_y_types = self.optionmenu_typey.keys()
        x_types,ind = np.unique(set_type_desc[(all_year_ranges == current_yearrange) & (all_spec == self.speciesx_str.get())],return_index=True)
        x_types = list(x_types[np.argsort(ind)])
        y_types,ind = np.unique(set_type_desc[(all_year_ranges == current_yearrange) & (all_spec == self.speciesy_str.get())],return_index=True)
        y_types = list(y_types[np.argsort(ind)])


        #update x and y types - if possible
        #only update if called from species update
        if (self.update_species_flag == True) & (self.update_year_flag == False):
            #update x types - if can
            if x_types != current_x_types:
                reset_flag = False
                #3 same species
                if self.species_flag == '3 same':
                     #Cannot be 3 modelgrids but all other options are fine
                     if (current_y_type != 'modelgrid') or (current_z_type != 'modelgrid'):
                        reset_flag = True
                #3 diff species
                elif self.species_flag == '3 diff':
                     #Must be only 1 obs or modelxy type, #if y or z are obs or modelxy, x must be grid
                     if all([current_y_type != 'obs', current_y_type == 'modelxy']) & all([current_z_type != 'obs', current_z_type == 'modelxy']):
                        reset_flag = True  
                #1,2 species
                else:
                     #Must be only 1 obs or modelxy type, #if y or z are obs or modelxy, x must be grid
                     if all([current_y_type != 'obs', current_y_type == 'modelxy']) & all([current_z_type != 'obs', current_z_type == 'modelxy']):
                        reset_flag = True  
                    
                #reset x types if can
                if reset_flag == True:
                    #reset y types
                     self.typex_str.set(x_types[0])
                     menu = self.optionmenu_typex['menu']
                     menu.delete(0, 'end')
                     for type in x_types:
                         menu.add_command(label=type, command=lambda s=type: self.typex_str.set(s))
                     #update x sets
                     sets = list(all_set[(all_year_ranges == current_yearrange) & (all_spec == current_x_species) & (set_type_desc == x_types[0])])        
                     self.setx_str.set(sets[0])
                     menu = self.optionmenu_setx['menu']
                     menu.delete(0, 'end')
                     for set in sets:
                        menu.add_command(label=set, command=lambda s=set: self.setx_str.set(s))
                        
            #update y types - if can
            if y_types != current_y_types:
                reset_flag = False
                #3 same species
                if self.species_flag == '3 same':
                     #Cannot be 3 modelgrids but all other options are fine
                     if (current_x_type != 'modelgrid') or (current_z_type != 'modelgrid'):
                        reset_flag = True
                #3 diff species
                elif self.species_flag == '3 diff':
                     #Must be only 1 obs or modelxy type, #if x or z are obs or modelxy, y must be grid
                     if all([current_x_type != 'obs', current_x_type == 'modelxy']) & all([current_z_type != 'obs', current_z_type == 'modelxy']):
                        reset_flag = True  
                #1,2 species
                else:
                     #Must be only 1 obs or modelxy type, #if x or z are obs or modelxy, y must be grid
                     if all([current_x_type != 'obs', current_x_type == 'modelxy']) & all([current_z_type != 'obs', current_z_type == 'modelxy']):
                        reset_flag = True  
                    
                #reset y types if can
                if reset_flag == True:
                    #reset y types
                     self.typey_str.set(y_types[0])
                     menu = self.optionmenu_typey['menu']
                     menu.delete(0, 'end')
                     for type in y_types:
                         menu.add_command(label=type, command=lambda s=type: self.typey_str.set(s))
                     #update y sets
                     sets = list(all_set[(all_year_ranges == current_yearrange) & (all_spec == current_y_species) & (set_type_desc == y_types[0])])        
                     self.sety_str.set(sets[0])
                     menu = self.optionmenu_sety['menu']
                     menu.delete(0, 'end')
                     for set in sets:
                        menu.add_command(label=set, command=lambda s=set: self.sety_str.set(s))

    
    def change_groupx(self, *args):
        #update params
        params = param_dict[self.groupx_str.get()]
        
        self.paramx_str.set(params[0])
        menu = self.optionmenu_paramx['menu']
        menu.delete(0, 'end')
        for p in params:
            menu.add_command(label=p, command=lambda s=p: self.paramx_str.set(s))
    
    def change_groupy(self, *args):
        #update params
        params = param_dict[self.groupy_str.get()]
        self.paramy_str.set(params[0])
        menu = self.optionmenu_paramy['menu']
        menu.delete(0, 'end')
        for p in params:
            menu.add_command(label=p, command=lambda s=p: self.paramy_str.set(s))
            
    def change_groupz(self, *args):
        #update params
        params = param_dict[self.groupz_str.get()]
        self.paramz_str.set(params[0])
        menu = self.optionmenu_paramz['menu']
        menu.delete(0, 'end')
        for p in params:
            menu.add_command(label=p, command=lambda s=p: self.paramz_str.set(s))
    
    
    def update(self):
        #update plot
        plot = update_data(self)

class update_data(object):
    def __init__(self,app): 
        self.app_class = app
        self.new_x_fname = all_filenames[(all_year_ranges == self.app_class.yearrange_str.get()) & (all_spec == self.app_class.speciesx_str.get()) & (all_set == self.app_class.setx_str.get()) & (set_type_desc == self.app_class.typex_str.get())][0]
        self.new_y_fname = all_filenames[(all_year_ranges == self.app_class.yearrange_str.get()) & (all_spec == self.app_class.speciesy_str.get()) & (all_set == self.app_class.sety_str.get()) & (set_type_desc == self.app_class.typey_str.get())][0]
        self.new_z_fname = all_filenames[(all_year_ranges == self.app_class.yearrange_str.get()) & (all_spec == self.app_class.speciesz_str.get()) & (all_set == self.app_class.setz_str.get()) & (set_type_desc == self.app_class.typez_str.get())][0]    
        self.process_read()
        self.key_vals()
        
    def process_read(self):
        #initialise old filenames if first time updating plots
        try:
            if old_x_fname:
                self.first_pass_flag = False
        except:
            self.first_pass_flag = True
            old_x_fname = ''
            old_y_fname = ''
            old_z_fname = ''
            global old_x_fname
            global old_y_fname
            global old_z_fname
            
        #setup instance variables that determine if need to update global variables at end
        self.set_new_meta = False
        self.set_new_x = False
        self.set_new_y = False
        self.set_new_z = False
        self.x_model_flag = False
        self.y_model_flag = False
        self.z_model_flag = False
        
        #determine if need to read x
        if (old_x_fname == self.new_x_fname):
            print 'dont read x - have old x'
            self.x_params = x_params[:] 
            self.x_ts = x_ts[:]
        elif (old_y_fname == self.new_x_fname):
            print 'dont read x - have old y'
            self.x_params = y_params[:]
            self.x_ts = y_ts[:]
        elif (old_z_fname == self.new_x_fname):
            print 'dont read x - have old z'
            self.x_params = z_params[:]
            self.x_ts = z_ts[:]
        else:
            print 'read new x'
            self.set_new_x = True
            if self.app_class.typex_str.get() == 'obs':
                self.set_new_meta = True
                get_meta_file = self.new_x_fname[:-13]
                obs_fname,species,start_year,end_year,vres,timeres,run_type = modules.get_obs_info(get_meta_file)
                self.refs,na,na,self.datetime_time,self.x_ts,self.lats,self.lons,self.alt,self.process_groups,self.raw_class,self.anthrome_class,na = modules.read_obs_all(obs_fname,species,start_year,end_year) 
            elif self.app_class.typex_str.get() == 'modelxy':   
                self.set_new_meta = True
                self.x_model_flag = True
                get_meta_file = ''
                for i in self.new_x_fname.split('/')[:-2]:
                    get_meta_file+=(i+'/')
                get_meta_file = get_meta_file+'obs_SURFACE_H'   
                obs_fname,species,start_year,end_year,vres,timeres,run_type = modules.get_obs_info(get_meta_file)
                self.refs,na,na,self.datetime_time,na,self.lats,self.lons,self.alt,self.process_groups,self.raw_class,self.anthrome_class,na = modules.read_obs_all(obs_fname,species,start_year,end_year) 
                model_fname,species,start_year,end_year,vres,timeres = modules.get_model_info(self.new_x_fname[:-13])
                na,na,self.datetime_time,self.x_ts,self.lat_e,self.lon_e,self.lat_c,self.lon_c,na,na = modules.read_model_all(model_fname,species,start_year,end_year)    
                na,na,na,lat_i,lon_i = modules.obs_model_gridboxes(self.lat_e,self.lon_e,self.lats,self.lons)
                self.x_ts = self.x_ts[:,lat_i,lon_i]
            elif self.app_class.typex_str.get() == 'modelgrid':
                self.x_model_flag = True
                model_fname,species,start_year,end_year,vres,timeres = modules.get_model_info(self.new_x_fname[:-13])
                na,na,self.datetime_time,self.x_ts,self.lat_e,self.lon_e,self.lat_c,self.lon_c,na,na = modules.read_model_all(model_fname,species,start_year,end_year)
                try:
                    na,na,na,lat_i,lon_i = modules.obs_model_gridboxes(self.lat_e,self.lon_e,self.lats,self.lons)
                except:
                    na,na,na,lat_i,lon_i = modules.obs_model_gridboxes(self.lat_e,self.lon_e,lats,lons)
                self.x_ts = self.x_ts[:,lat_i,lon_i]
            self.x_params = self.read_data(self.new_x_fname,self.app_class.typex_str.get())

        #determine if need to read y
        if (self.new_x_fname) == (self.new_y_fname):
            print 'dont read y - have new x'
            self.y_params = self.x_params[:] 
            self.y_ts = self.x_ts[:]
        elif (old_x_fname == self.new_y_fname):
            print 'dont read y - have old x'
            self.y_params = x_params[:] 
            self.y_ts = x_ts[:]
        elif (old_y_fname == self.new_y_fname):
            print 'dont read y - have old y'
            self.y_params = y_params[:]
            self.y_ts = y_ts[:]
        elif (old_z_fname == self.new_y_fname):
            print 'dont read y - have old z'
            self.y_params = z_params[:]
            self.y_ts = z_ts[:]
        else:
            print 'read new y'
            self.set_new_y = True
            if self.app_class.typey_str.get() == 'obs':
                self.set_new_meta = True
                get_meta_file = self.new_y_fname[:-13]
                obs_fname,species,start_year,end_year,vres,timeres,run_type = modules.get_obs_info(get_meta_file)
                self.refs,na,na,self.datetime_time,self.y_ts,self.lats,self.lons,self.alt,self.process_groups,self.raw_class,self.anthrome_class,na = modules.read_obs_all(obs_fname,species,start_year,end_year) 
            elif self.app_class.typey_str.get() == 'modelxy':   
                self.set_new_meta = True
                self.y_model_flag = True
                get_meta_file = ''
                for i in self.new_y_fname.split('/')[:-2]:
                    get_meta_file+=(i+'/')
                get_meta_file = get_meta_file+'obs_SURFACE_H'   
                obs_fname,species,start_year,end_year,vres,timeres,run_type = modules.get_obs_info(get_meta_file)
                self.refs,na,na,self.datetime_time,na,self.lats,self.lons,self.alt,self.process_groups,self.raw_class,self.anthrome_class,na = modules.read_obs_all(obs_fname,species,start_year,end_year) 
                model_fname,species,start_year,end_year,vres,timeres = modules.get_model_info(self.new_y_fname[:-13])
                na,na,self.datetime_time,self.y_ts,self.lat_e,self.lon_e,self.lat_c,self.lon_c,na,na = modules.read_model_all(model_fname,species,start_year,end_year)
                na,na,na,lat_i,lon_i = modules.obs_model_gridboxes(self.lat_e,self.lon_e,self.lats,self.lons)
                self.y_ts = self.y_ts[:,lat_i,lon_i]
            elif self.app_class.typey_str.get() == 'modelgrid':
                self.y_model_flag = True
                model_fname,species,start_year,end_year,vres,timeres = modules.get_model_info(self.new_y_fname[:-13])
                na,na,self.datetime_time,self.y_ts,self.lat_e,self.lon_e,self.lat_c,self.lon_c,na,na = modules.read_model_all(model_fname,species,start_year,end_year)
                try:
                    na,na,na,lat_i,lon_i = modules.obs_model_gridboxes(self.lat_e,self.lon_e,self.lats,self.lons)
                except:
                    na,na,na,lat_i,lon_i = modules.obs_model_gridboxes(self.lat_e,self.lon_e,lats,lons)
                self.y_ts = self.y_ts[:,lat_i,lon_i]
            self.y_params = self.read_data(self.new_y_fname,self.app_class.typey_str.get())
            
        #determine if need to read z
        if (self.new_x_fname) == (self.new_z_fname):
            print 'dont read z - have new x'
            self.z_params = self.x_params[:]
            self.z_ts = self.x_ts[:]
        elif (self.new_y_fname) == (self.new_z_fname):
            print 'dont read z - have new y'
            self.z_params = self.y_params[:]
            self.z_ts = self.y_ts[:]
        elif (old_x_fname == self.new_z_fname):
            print 'dont read z - have old x'
            self.z_params = x_params[:]
            self.z_ts = x_ts[:]
        elif (old_y_fname == self.new_z_fname):
            print 'dont read z - have old y'
            self.z_params = y_params[:]
            self.z_ts = y_ts[:]
        elif (old_z_fname == self.new_z_fname):
            print 'dont read z - have old z'
            self.z_params = z_params[:]
            self.z_ts = z_ts[:]
        else:
            print 'read new z'
            self.set_new_z = True            
            if self.app_class.typez_str.get() == 'obs':
                self.set_new_meta = True
                get_meta_file = self.new_z_fname[:-13]
                obs_fname,species,start_year,end_year,vres,timeres,run_type = modules.get_obs_info(get_meta_file)
                self.refs,na,na,self.datetime_time,self.z_ts,self.lats,self.lons,self.alt,self.process_groups,self.raw_class,self.anthrome_class,na = modules.read_obs_all(obs_fname,species,start_year,end_year) 
            elif self.app_class.typez_str.get() == 'modelxy':   
                self.set_new_meta = True
                self.z_model_flag = True
                get_meta_file = ''
                for i in self.new_z_fname.split('/')[:-2]:
                    get_meta_file+=(i+'/')
                get_meta_file = get_meta_file+'obs_SURFACE_H'   
                obs_fname,species,start_year,end_year,vres,timeres,run_type = modules.get_obs_info(get_meta_file)
                self.refs,na,na,self.datetime_time,na,self.lats,self.lons,self.alt,self.process_groups,self.raw_class,self.anthrome_class,na = modules.read_obs_all(obs_fname,species,start_year,end_year) 
                model_fname,species,start_year,end_year,vres,timeres = modules.get_model_info(self.new_z_fname[:-13])
                na,na,self.datetime_time,self.z_ts,self.lat_e,self.lon_e,self.lat_c,self.lon_c,na,na = modules.read_model_all(model_fname,species,start_year,end_year)
                na,na,na,lat_i,lon_i = modules.obs_model_gridboxes(self.lat_e,self.lon_e,self.lats,self.lons)
                self.z_ts = self.z_ts[:,lat_i,lon_i]
            elif self.app_class.typez_str.get() == 'modelgrid':
                self.z_model_flag = True
                model_fname,species,start_year,end_year,vres,timeres = modules.get_model_info(self.new_z_fname[:-13])
                na,na,self.datetime_time,self.z_ts,self.lat_e,self.lon_e,self.lat_c,self.lon_c,na,na = modules.read_model_all(model_fname,species,start_year,end_year)
                try:
                    na,na,na,lat_i,lon_i = modules.obs_model_gridboxes(self.lat_e,self.lon_e,self.lats,self.lons)
                except:
                    na,na,na,lat_i,lon_i = modules.obs_model_gridboxes(self.lat_e,self.lon_e,lats,lons)
                self.z_ts = self.z_ts[:,lat_i,lon_i]
            self.z_params = self.read_data(self.new_z_fname,self.app_class.typez_str.get())

        #get specific x param
        if self.app_class.groupx_str.get() == 'standard':
            self.x_param_name = self.app_class.paramx_str.get()
        else:
            self.x_param_name = self.app_class.paramx_str.get()+'_'+self.app_class.groupx_str.get()
        x_ind = params_str.index(self.x_param_name)
        self.x = self.x_params[x_ind]
        
        #get specific y param
        if self.app_class.groupy_str.get() == 'standard':
            self.y_param_name = self.app_class.paramy_str.get()
        else:
            self.y_param_name = self.app_class.paramy_str.get()+'_'+self.app_class.groupy_str.get()
        y_ind = params_str.index(self.y_param_name)
        self.y = self.y_params[y_ind]
        
        #get specific z param
        if self.app_class.groupz_str.get() == 'standard':
            self.z_param_name = self.app_class.paramz_str.get()
        else:
            self.z_param_name = self.app_class.paramz_str.get()+'_'+self.app_class.groupz_str.get()
        z_ind = params_str.index(self.z_param_name)
        self.z = self.z_params[z_ind]
        
        #print self.refs[np.isnan(self.z)]
        #print len((self.lats))
        #print self.x[:10],np.min(self.z),np.max(self.z),len(self.x)
        #print self.y[:10],np.min(self.z),np.max(self.z),len(self.y)
        #print self.z[:10],np.min(self.z),np.max(self.z),len(self.z)
        
        #get diurnal/seasonal waveforms
        #x
        if (self.app_class.groupx_str.get() == 'standard'):
            self.x_diurnal_waveform = self.x_params[params_str.index('diurnal_ave_waveform')]
            self.x_seasonal_waveform = self.x_params[params_str.index('seasonal_waveform')]
        elif self.app_class.groupx_str.get() == 'spring':
            self.x_diurnal_waveform = self.x_params[params_str.index('diurnal_waveform_spring')]
            self.x_seasonal_waveform = self.x_params[params_str.index('seasonal_waveform')]
        elif self.app_class.groupx_str.get() == 'summer':
            self.x_diurnal_waveform = self.x_params[params_str.index('diurnal_waveform_summer')]
            self.x_seasonal_waveform = self.x_params[params_str.index('seasonal_waveform')]
        elif self.app_class.groupx_str.get() == 'autumn':
            self.x_diurnal_waveform = self.x_params[params_str.index('diurnal_waveform_autumn')]
            self.x_seasonal_waveform = self.x_params[params_str.index('seasonal_waveform')]
        elif self.app_class.groupx_str.get() == 'winter':
            self.x_diurnal_waveform = self.x_params[params_str.index('diurnal_waveform_winter')]
            self.x_seasonal_waveform = self.x_params[params_str.index('seasonal_waveform')]
        elif (self.app_class.groupx_str.get() == 'day'):
            self.x_diurnal_waveform = self.x_params[params_str.index('diurnal_ave_waveform')]
            self.x_seasonal_waveform = self.x_params[params_str.index('seasonal_waveform_day')]
        elif (self.app_class.groupx_str.get() == 'night'):
            self.x_diurnal_waveform = self.x_params[params_str.index('diurnal_ave_waveform')]
            self.x_seasonal_waveform = self.x_params[params_str.index('seasonal_waveform_night')]
            
        #y
        if (self.app_class.groupy_str.get() == 'standard'):
            self.y_diurnal_waveform = self.y_params[params_str.index('diurnal_ave_waveform')]
            self.y_seasonal_waveform = self.y_params[params_str.index('seasonal_waveform')]
        elif self.app_class.groupy_str.get() == 'spring':
            self.y_diurnal_waveform = self.y_params[params_str.index('diurnal_waveform_spring')]
            self.y_seasonal_waveform = self.y_params[params_str.index('seasonal_waveform')]
        elif self.app_class.groupy_str.get() == 'summer':
            self.y_diurnal_waveform = self.y_params[params_str.index('diurnal_waveform_summer')]
            self.y_seasonal_waveform = self.y_params[params_str.index('seasonal_waveform')]
        elif self.app_class.groupy_str.get() == 'autumn':
            self.y_diurnal_waveform = self.y_params[params_str.index('diurnal_waveform_autumn')]
            self.y_seasonal_waveform = self.y_params[params_str.index('seasonal_waveform')]
        elif self.app_class.groupy_str.get() == 'winter':
            self.y_diurnal_waveform = self.y_params[params_str.index('diurnal_waveform_winter')]
            self.y_seasonal_waveform = self.y_params[params_str.index('seasonal_waveform')]
        elif (self.app_class.groupy_str.get() == 'day'):
            self.y_diurnal_waveform = self.y_params[params_str.index('diurnal_ave_waveform')]
            self.y_seasonal_waveform = self.y_params[params_str.index('seasonal_waveform_day')]
        elif (self.app_class.groupy_str.get() == 'night'):
            self.y_diurnal_waveform = self.y_params[params_str.index('diurnal_ave_waveform')]
            self.y_seasonal_waveform = self.y_params[params_str.index('seasonal_waveform_night')]
            
        #z
        if (self.app_class.groupz_str.get() == 'standard'):
            self.z_diurnal_waveform = self.z_params[params_str.index('diurnal_ave_waveform')]
            self.z_seasonal_waveform = self.z_params[params_str.index('seasonal_waveform')]
        elif self.app_class.groupz_str.get() == 'spring':
            self.z_diurnal_waveform = self.z_params[params_str.index('diurnal_waveform_spring')]
            self.z_seasonal_waveform = self.z_params[params_str.index('seasonal_waveform')]
        elif self.app_class.groupz_str.get() == 'summer':
            self.z_diurnal_waveform = self.z_params[params_str.index('diurnal_waveform_summer')]
            self.z_seasonal_waveform = self.z_params[params_str.index('seasonal_waveform')]
        elif self.app_class.groupz_str.get() == 'autumn':
            self.z_diurnal_waveform = self.z_params[params_str.index('diurnal_waveform_autumn')]
            self.z_seasonal_waveform = self.z_params[params_str.index('seasonal_waveform')]
        elif self.app_class.groupz_str.get() == 'winter':
            self.z_diurnal_waveform = self.z_params[params_str.index('diurnal_waveform_winter')]
            self.z_seasonal_waveform = self.z_params[params_str.index('seasonal_waveform')]
        elif (self.app_class.groupz_str.get() == 'day'):
            self.z_diurnal_waveform = self.z_params[params_str.index('diurnal_ave_waveform')]
            self.z_seasonal_waveform = self.z_params[params_str.index('seasonal_waveform_day')]
        elif (self.app_class.groupz_str.get() == 'night'):
            self.z_diurnal_waveform = self.z_params[params_str.index('diurnal_ave_waveform')]
            self.z_seasonal_waveform = self.z_params[params_str.index('seasonal_waveform_night')]
            
        #update scatter plot
        self.plot()

    def read_data(self,filename,readtype):
        cut_filename = filename[:-13]
        
        if readtype == 'obs':        
            (diurnal_mag,diurnal_ph,seasonal_mag,seasonal_ph,mean,p1,p5,p25,p50,p75,p95,p99,diurnal_ave_waveform,seasonal_waveform,full_ave_waveform,pc_var_daily,pc_var_seasonal,pc_var_full,pc_var_noise,total_var,
            diurnal_mag_spring,diurnal_ph_spring,mean_spring,p1_spring,p5_spring,p25_spring,p50_spring,p75_spring,p95_spring,p99_spring,diurnal_waveform_spring,
            diurnal_mag_summer,diurnal_ph_summer,mean_summer,p1_summer,p5_summer,p25_summer,p50_summer,p75_summer,p95_summer,p99_summer,diurnal_waveform_summer,
            diurnal_mag_autumn,diurnal_ph_autumn,mean_autumn,p1_autumn,p5_autumn,p25_autumn,p50_autumn,p75_autumn,p95_autumn,p99_autumn,diurnal_waveform_autumn,
            diurnal_mag_winter,diurnal_ph_winter,mean_winter,p1_winter,p5_winter,p25_winter,p50_winter,p75_winter,p95_winter,p99_winter,diurnal_waveform_winter,
            seasonal_mag_day,seasonal_ph_day,mean_day,p1_day,p5_day,p25_day,p50_day,p75_day,p95_day,p99_day,seasonal_waveform_day,
            seasonal_mag_night,seasonal_ph_night,mean_night,p1_night,p5_night,p25_night,p50_night,p75_night,p95_night,p99_night,seasonal_waveform_night,
            daily_h3_mag,daily_h2_mag,daily_h1_mag,daily_mag,annual_h3_mag,annual_h2_mag,annual_h1_mag,annual_mag) = modules.get_periodic_specific(filename,'na')

        elif (readtype == 'modelxy') or (readtype == 'modelgrid'):        
            (diurnal_mag,diurnal_ph,seasonal_mag,seasonal_ph,mean,p1,p5,p25,p50,p75,p95,p99,diurnal_ave_waveform,seasonal_waveform,full_ave_waveform,pc_var_daily,pc_var_seasonal,pc_var_full,pc_var_noise,total_var,
            diurnal_mag_spring,diurnal_ph_spring,mean_spring,p1_spring,p5_spring,p25_spring,p50_spring,p75_spring,p95_spring,p99_spring,diurnal_waveform_spring,
            diurnal_mag_summer,diurnal_ph_summer,mean_summer,p1_summer,p5_summer,p25_summer,p50_summer,p75_summer,p95_summer,p99_summer,diurnal_waveform_summer,
            diurnal_mag_autumn,diurnal_ph_autumn,mean_autumn,p1_autumn,p5_autumn,p25_autumn,p50_autumn,p75_autumn,p95_autumn,p99_autumn,diurnal_waveform_autumn,
            diurnal_mag_winter,diurnal_ph_winter,mean_winter,p1_winter,p5_winter,p25_winter,p50_winter,p75_winter,p95_winter,p99_winter,diurnal_waveform_winter,
            seasonal_mag_day,seasonal_ph_day,mean_day,p1_day,p5_day,p25_day,p50_day,p75_day,p95_day,p99_day,seasonal_waveform_day,
            seasonal_mag_night,seasonal_ph_night,mean_night,p1_night,p5_night,p25_night,p50_night,p75_night,p95_night,p99_night,seasonal_waveform_night,
            daily_h3_mag,daily_h2_mag,daily_h1_mag,daily_mag,annual_h3_mag,annual_h2_mag,annual_h1_mag,annual_mag) = modules.get_periodic_specific('na',filename)  

        #need to cut modelgrid params to obs lat/lons
        if readtype == 'modelgrid':
            try:
                na,na,na,lat_i,lon_i = modules.obs_model_gridboxes(self.lat_e,self.lon_e,self.lats,self.lons)
            except:
                na,na,na,lat_i,lon_i = modules.obs_model_gridboxes(self.lat_e,self.lon_e,lats,lons)

            return [diurnal_mag[lat_i,lon_i],diurnal_ph[lat_i,lon_i],seasonal_mag[lat_i,lon_i],seasonal_ph[lat_i,lon_i],mean[lat_i,lon_i],p1[lat_i,lon_i],p5[lat_i,lon_i],p25[lat_i,lon_i],p50[lat_i,lon_i],p75[lat_i,lon_i],p95[lat_i,lon_i],p99[lat_i,lon_i],diurnal_ave_waveform[:,lat_i,lon_i],seasonal_waveform[:,lat_i,lon_i],full_ave_waveform[:,lat_i,lon_i],pc_var_daily[lat_i,lon_i],pc_var_seasonal[lat_i,lon_i],pc_var_full[lat_i,lon_i],pc_var_noise[lat_i,lon_i],total_var[lat_i,lon_i],
            diurnal_mag_spring[lat_i,lon_i],diurnal_ph_spring[lat_i,lon_i],mean_spring[lat_i,lon_i],p1_spring[lat_i,lon_i],p5_spring[lat_i,lon_i],p25_spring[lat_i,lon_i],p50_spring[lat_i,lon_i],p75_spring[lat_i,lon_i],p95_spring[lat_i,lon_i],p99_spring[lat_i,lon_i],diurnal_waveform_spring[:,lat_i,lon_i],
            diurnal_mag_summer[lat_i,lon_i],diurnal_ph_summer[lat_i,lon_i],mean_summer[lat_i,lon_i],p1_summer[lat_i,lon_i],p5_summer[lat_i,lon_i],p25_summer[lat_i,lon_i],p50_summer[lat_i,lon_i],p75_summer[lat_i,lon_i],p95_summer[lat_i,lon_i],p99_summer[lat_i,lon_i],diurnal_waveform_summer[:,lat_i,lon_i],
            diurnal_mag_autumn[lat_i,lon_i],diurnal_ph_autumn[lat_i,lon_i],mean_autumn[lat_i,lon_i],p1_autumn[lat_i,lon_i],p5_autumn[lat_i,lon_i],p25_autumn[lat_i,lon_i],p50_autumn[lat_i,lon_i],p75_autumn[lat_i,lon_i],p95_autumn[lat_i,lon_i],p99_autumn[lat_i,lon_i],diurnal_waveform_autumn[:,lat_i,lon_i],
            diurnal_mag_winter[lat_i,lon_i],diurnal_ph_winter[lat_i,lon_i],mean_winter[lat_i,lon_i],p1_winter[lat_i,lon_i],p5_winter[lat_i,lon_i],p25_winter[lat_i,lon_i],p50_winter[lat_i,lon_i],p75_winter[lat_i,lon_i],p95_winter[lat_i,lon_i],p99_winter[lat_i,lon_i],diurnal_waveform_winter[:,lat_i,lon_i],
            seasonal_mag_day[lat_i,lon_i],seasonal_ph_day[lat_i,lon_i],mean_day[lat_i,lon_i],p1_day[lat_i,lon_i],p5_day[lat_i,lon_i],p25_day[lat_i,lon_i],p50_day[lat_i,lon_i],p75_day[lat_i,lon_i],p95_day[lat_i,lon_i],p99_day[lat_i,lon_i],seasonal_waveform_day[:,lat_i,lon_i],
            seasonal_mag_night[lat_i,lon_i],seasonal_ph_night[lat_i,lon_i],mean_night[lat_i,lon_i],p1_night[lat_i,lon_i],p5_night[lat_i,lon_i],p25_night[lat_i,lon_i],p50_night[lat_i,lon_i],p75_night[lat_i,lon_i],p95_night[lat_i,lon_i],p99_night[lat_i,lon_i],seasonal_waveform_night[:,lat_i,lon_i],
            daily_h3_mag[lat_i,lon_i],daily_h2_mag[lat_i,lon_i],daily_h1_mag[lat_i,lon_i],daily_mag[lat_i,lon_i],annual_h3_mag[lat_i,lon_i],annual_h2_mag[lat_i,lon_i],annual_h1_mag[lat_i,lon_i],annual_mag[lat_i,lon_i]]

        else:
            return [diurnal_mag,diurnal_ph,seasonal_mag,seasonal_ph,mean,p1,p5,p25,p50,p75,p95,p99,diurnal_ave_waveform,seasonal_waveform,full_ave_waveform,pc_var_daily,pc_var_seasonal,pc_var_full,pc_var_noise,total_var,
            diurnal_mag_spring,diurnal_ph_spring,mean_spring,p1_spring,p5_spring,p25_spring,p50_spring,p75_spring,p95_spring,p99_spring,diurnal_waveform_spring,
            diurnal_mag_summer,diurnal_ph_summer,mean_summer,p1_summer,p5_summer,p25_summer,p50_summer,p75_summer,p95_summer,p99_summer,diurnal_waveform_summer,
            diurnal_mag_autumn,diurnal_ph_autumn,mean_autumn,p1_autumn,p5_autumn,p25_autumn,p50_autumn,p75_autumn,p95_autumn,p99_autumn,diurnal_waveform_autumn,
            diurnal_mag_winter,diurnal_ph_winter,mean_winter,p1_winter,p5_winter,p25_winter,p50_winter,p75_winter,p95_winter,p99_winter,diurnal_waveform_winter,
            seasonal_mag_day,seasonal_ph_day,mean_day,p1_day,p5_day,p25_day,p50_day,p75_day,p95_day,p99_day,seasonal_waveform_day,
            seasonal_mag_night,seasonal_ph_night,mean_night,p1_night,p5_night,p25_night,p50_night,p75_night,p95_night,p99_night,seasonal_waveform_night,
            daily_h3_mag,daily_h2_mag,daily_h1_mag,daily_mag,annual_h3_mag,annual_h2_mag,annual_h1_mag,annual_mag]

    def plot(self):
        #clear scatter and map plot
        ax1.cla()
        ax5.cla()
        m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=-180,urcrnrlon=180,resolution='c',ax=ax5)
        m.drawcoastlines()
        m.drawmapboundary()
        
        #plot new scatter and map plot
        if '_ph' not in self.z_param_name:
            ax1.scatter(self.x,self.y,s=100,c=self.z,picker=1,vmin=np.min(self.z),vmax=np.max(self.z),cmap=cmaps.viridis,edgecolor='none')
            #have exception for cases when lats/lons not read and no self instance, so revert to global variables of previous read 
            try:
                ax5.scatter(self.lons,self.lats,s=25,c=self.z,picker=1,vmin=np.min(self.z),vmax=np.max(self.z),cmap=cmaps.viridis,edgecolor='none')
            except:
                ax5.scatter(lons,lats,s=25,c=self.z,picker=1,vmin=np.min(self.z),vmax=np.max(self.z),cmap=cmaps.viridis,edgecolor='none')
            
            #change colorbar range
            cb.set_clim(vmin=np.min(self.z),vmax=np.max(self.z)) 
            #set colormap
            cb.set_cmap(cmaps.viridis)
        else:
            z_min = 0
            if 'diurnal' in self.z_param_name.lower():
                z_max = 24
                #change colorbar range
                cb.set_clim(vmin=z_min,vmax=z_max)
            elif 'seasonal' in self.z_param_name.lower():
                z_max = 12
                #change colorbar range
                cb.set_clim(vmin=z_min,vmax=z_max)
            ax1.scatter(self.x,self.y,s=100,c=self.z,picker=1,vmin=z_min,vmax=z_max,cmap=plt.cm.hsv,edgecolor='none')
            #have exception for cases when lats/lons not read and no self instance, so revert to global variables of previous read
            try:
                ax5.scatter(self.lons,self.lats,s=25,c=self.z,picker=1,vmin=z_min,vmax=z_max,cmap=plt.cm.hsv,edgecolor='none')
            except:
                ax5.scatter(lons,lats,s=25,c=self.z,picker=1,vmin=z_min,vmax=z_max,cmap=plt.cm.hsv,edgecolor='none')
        
            #set colormap
            cb.set_cmap(plt.cm.hsv)                    
    
        #set xy axis limits
        ax1.set_xlim(np.min(self.x)-1,np.max(self.x)+1)
        ax1.set_ylim(np.min(self.y)-1,np.max(self.y)+1)
        
        #plot one to one line if x and y param names are same
        if self.x_param_name == self.y_param_name:
            one_to_one_line = np.arange(np.min(np.append(self.x,self.y))-1,np.max(np.append(self.x,self.y))+1,1.)
            ax1.plot(one_to_one_line,one_to_one_line,color='black',linestyle='--')
    
        #set colorbar label
        cb.set_label('%s %s\n%s %s'%(self.app_class.speciesz_str.get(),self.app_class.setz_str.get(),self.app_class.groupz_str.get(),self.app_class.paramz_str.get()))
        #set axis labels
        ax1.set_xlabel('%s %s\n%s %s'%(self.app_class.speciesx_str.get(),self.app_class.setx_str.get(),self.app_class.groupx_str.get(),self.app_class.paramx_str.get()))
        ax1.set_ylabel('%s %s\n%s %s'%(self.app_class.speciesy_str.get(),self.app_class.sety_str.get(),self.app_class.groupy_str.get(),self.app_class.paramy_str.get()))
        #redraw
        cb.draw_all() 
        fig.canvas.draw() 
     
    def key_vals(self):
        
        #set global file dependent variables if first run
        if self.first_pass_flag == True:
            datetime_time = self.datetime_time
            refs = self.refs
            lats = self.lats
            lons = self.lons
            alt = self.alt
            global datetime_time
            global refs
            global lats
            global lons
            global alt
        
        #set global file dependent variables if file changed
        else:
            if self.set_new_x == True:
                datetime_time = self.datetime_time
                global datetime_time
                
            if self.set_new_y == True:
                datetime_time = self.datetime_time
                global datetime_time

            if self.set_new_z == True:
                datetime_time = self.datetime_time
                global datetime_time
                
            if self.set_new_meta == True:
                refs = self.refs
                lats = self.lats
                lons = self.lons
                alt = self.alt
                global refs
                global lats
                global lons
                global alt
        
        #these variables will always be processed so always update
        x = self.x
        x_diurnal_waveform = self.x_diurnal_waveform
        x_seasonal_waveform = self.x_seasonal_waveform
        x_ts = self.x_ts
        x_species = self.app_class.speciesx_str.get()
        x_type = self.app_class.typex_str.get()
        x_set = self.app_class.setx_str.get()
        x_group = self.app_class.groupx_str.get()
        x_param = self.app_class.paramx_str.get()
        x_params = self.x_params
        y = self.y
        y_diurnal_waveform = self.y_diurnal_waveform
        y_seasonal_waveform = self.y_seasonal_waveform
        y_ts = self.y_ts
        y_species = self.app_class.speciesy_str.get()
        y_type = self.app_class.typey_str.get()
        y_set = self.app_class.sety_str.get()
        y_group = self.app_class.groupy_str.get()
        y_param = self.app_class.paramy_str.get()
        y_params = self.y_params
        z = self.z
        z_diurnal_waveform = self.z_diurnal_waveform
        z_seasonal_waveform = self.z_seasonal_waveform
        z_ts = self.z_ts
        z_species = self.app_class.speciesz_str.get()
        z_type = self.app_class.typez_str.get()
        z_set = self.app_class.setz_str.get()
        z_group = self.app_class.groupz_str.get()
        z_param = self.app_class.paramz_str.get()
        z_params = self.z_params
        global x
        global x_diurnal_waveform
        global x_seasonal_waveform
        global x_ts
        global x_species
        global x_type
        global x_set
        global x_group
        global x_param
        global x_params
        global y 
        global y_diurnal_waveform
        global y_seasonal_waveform
        global y_ts
        global y_species
        global y_type
        global y_set
        global y_group
        global y_param
        global y_params
        global z 
        global z_diurnal_waveform
        global z_seasonal_waveform
        global z_ts
        global z_species
        global z_type
        global z_set
        global z_group
        global z_param
        global z_params
        old_x_fname = self.new_x_fname
        old_y_fname = self.new_y_fname
        old_z_fname = self.new_z_fname
        global old_x_fname
        global old_y_fname
        global old_z_fname
        
#---------------------------------------------------------
#SETUP CALLBACKS FOR POINT SELECT AND LASSO SELECT

class select_callbacks(object):
    def __init__(self):
        self.canvas = ax1.figure.canvas
        self.pickEvent = False
        self.bid = self.canvas.mpl_connect('pick_event', self.click_callback)
        self.aid = self.canvas.mpl_connect('button_press_event', self.onpress)

    def lasso_callback(self, verts):
        p = path.Path(verts)
        ind = p.contains_points(self.xys)
        try:
            self.selected.remove()
            self.selectedmap.remove()
        except:
            pass
        self.selected, =ax1.plot(x[ind], y[ind], 'o', ms=16, alpha=0.6,color='yellow', visible=True,zorder=10)
        self.selectedmap, =ax5.plot(lons[ind], lats[ind], 'o', ms=9, alpha=0.6,color='yellow', visible=True,zorder=10)
        
        ax2.cla()
        ax3.cla()
        ax4.cla()
        if 'obs' in x_type.lower():
            x_tsi = x_ts[ind,:]
            x_dw = x_diurnal_waveform[ind,:]
            x_sw = x_seasonal_waveform[ind,:]
            xaxis_ts_n = 0 
            xaxis_p_n = 0
        elif 'modelxy' == x_type.lower():
            x_tsi = x_ts[:,ind]
            x_dw = x_diurnal_waveform[ind,:]
            x_sw = x_seasonal_waveform[ind,:]
            xaxis_ts_n = 1
            xaxis_p_n = 0
        elif 'modelgrid' == x_type.lower():
            x_tsi = x_ts[ind,:]
            x_dw = x_diurnal_waveform[:,ind]
            x_sw = x_seasonal_waveform[:,ind]
            xaxis_ts_n = 1
            xaxis_p_n = 1
        if 'obs' in y_type.lower():
            y_tsi = y_ts[ind,:]
            y_dw = y_diurnal_waveform[ind,:]
            y_sw = y_seasonal_waveform[ind,:]
            yaxis_ts_n = 0 
            yaxis_p_n = 0
        elif 'modelxy' == y_type.lower():
            y_tsi = y_ts[:,ind]
            y_dw = y_diurnal_waveform[ind,:]
            y_sw = y_seasonal_waveform[ind,:]
            yaxis_ts_n = 1
            yaxis_p_n = 0
        elif 'modelgrid' == y_type.lower():
            y_tsi = y_ts[:,ind]
            y_dw = y_diurnal_waveform[:,ind]
            y_sw = y_seasonal_waveform[:,ind]
            yaxis_ts_n = 1
            yaxis_p_n = 1
            
        x_tsi[x_tsi < 0] = np.NaN 
        y_tsi[y_tsi < 0] = np.NaN 
        
        ax2.plot(datetime_time,np.nanmean(x_tsi,axis=xaxis_ts_n),color='black')
        ax2.plot(datetime_time,np.nanmean(y_tsi,axis=yaxis_ts_n),color='red',alpha=0.4)
        #ax2.set_title('Average time series of %s sites'%(),loc='right')
        
        ax3.plot(datetime_time[:24],np.nanmean(x_dw,axis=xaxis_p_n),color='black',label='%s %s %s'%(x_species,x_set,x_group))
        ax3.plot(datetime_time[:24],np.nanmean(y_dw,axis=yaxis_p_n),color='red',alpha=0.4,label='%s %s %s'%(y_species,y_set,y_group))
        leg = ax3.legend(loc=4,fontsize=10)
        leg.get_frame().set_alpha(0) 
        leg.get_frame().set_edgecolor('white')
        ax4.plot(datetime_time[:8766],np.nanmean(x_sw,axis=xaxis_p_n),color='black')
        ax4.plot(datetime_time[:8766],np.nanmean(y_sw,axis=yaxis_p_n),color='red',alpha=0.4)
        
        self.canvas.draw_idle()
        self.canvas.widgetlock.release(self.lasso)
        del self.lasso
    
    def click_callback(self, event):
        self.pickEvent = True
        if len(event.ind) > 0:
            ind = event.ind[0]
        else:
            return
            
        try:
            self.selected.remove()
            self.selectedmap.remove()
        except:
            pass
        self.selected, = ax1.plot(x[ind], y[ind], 'o', ms=16, alpha=0.6,color='yellow', visible=True,zorder=10)
        self.selectedmap, = ax5.plot(lons[ind], lats[ind], 'o', ms=9, alpha=0.6,color='yellow', visible=True,zorder=10)
        
        ax2.cla()
        ax3.cla()
        ax4.cla()
        if 'obs' in x_type.lower():
            x_tsi = x_ts[ind]
        else:
            x_tsi = x_ts[:,ind]
        if 'obs' in y_type.lower():
            y_tsi = y_ts[ind]
        else:
            y_tsi = y_ts[:,ind]
        x_tsi[x_tsi < 0] = np.NaN 
        y_tsi[y_tsi < 0] = np.NaN 
        ax2.plot(datetime_time,x_tsi,color='black')
        ax2.plot(datetime_time,y_tsi,color='red',alpha=0.4)
        ax2.set_title('%s Lat:%s Lon:%s Alt:%s'%(refs[ind],np.round(lats[ind],2),np.round(lons[ind],2),np.round(alt[ind],2)),loc='right')
        
        if 'modelgrid' == x_type.lower(): 
            x_dw = x_diurnal_waveform[:,ind]
            x_sw = x_seasonal_waveform[:,ind]
        else:
            x_dw = x_diurnal_waveform[ind]
            x_sw = x_seasonal_waveform[ind]
        if 'modelgrid' == y_type.lower():
            y_dw = y_diurnal_waveform[:,ind]
            y_sw = y_seasonal_waveform[:,ind]
        else:
            y_dw = y_diurnal_waveform[ind]
            y_sw = y_seasonal_waveform[ind]
        
        ax3.plot(datetime_time[:24],x_dw,color='black',label='%s %s %s'%(x_species,x_set,x_group))
        ax3.plot(datetime_time[:24],y_dw,color='red',alpha=0.4,label='%s %s %s'%(y_species,y_set,y_group))
        leg = ax3.legend(loc=4,fontsize=10)
        leg.get_frame().set_alpha(0) 
        leg.get_frame().set_edgecolor('white')
        ax4.plot(datetime_time[:8766],x_sw,color='black')
        ax4.plot(datetime_time[:8766],y_sw,color='red',alpha=0.4)
        
        self.canvas.draw_idle()

    def onpress(self, event):
        self.merge_scatter = [(xx, yy) for xx,yy in zip(x,y)]
        self.merge_map = [(xx, yy) for xx,yy in zip(lons,lats)]
        if event.inaxes == ax1:
            self.xys = self.merge_scatter
        elif event.inaxes == ax5:
            self.xys = self.merge_map
        self.canvas = event.canvas
        if self.pickEvent:
            self.pickEvent = False
            return
        if self.canvas.widgetlock.locked():
            return
        if event.inaxes is None:
            return
        self.lasso = Lasso(event.inaxes,(event.xdata, event.ydata),self.lasso_callback)
        # acquire a lock on the widget drawing
        self.canvas.widgetlock(self.lasso)

if __name__ == "__main__":
    root = Tk()
        
    fig =plt.figure(figsize=(17,8.5))
    fig.patch.set_facecolor('white')
    gs1 = gridspec.GridSpec(4, 4)
    gs1.update(left=0.06,right=0.50,top=0.98,bottom=0.38,wspace=0.03)
    ax1 = plt.subplot(gs1[:,:])
    gs2 = gridspec.GridSpec(4, 4)
    gs2.update(left=0.06,right=0.50,top=0.28,bottom=0.06,wspace=0.03)
    ax2 = plt.subplot(gs2[:,:]) 
    gs3 = gridspec.GridSpec(4, 4)
    gs3.update(left=0.55,right=0.975,bottom=0.26,top=0.42,wspace=0.03)
    ax3 = plt.subplot(gs3[:,:])         
    gs4 = gridspec.GridSpec(4, 4)
    gs4.update(left=0.55,right=0.975,top=0.23,bottom=0.06,wspace=0.03)
    ax4 = plt.subplot(gs4[:,:]) 
    gs5 = gridspec.GridSpec(4, 4)
    gs5.update(left=0.51,right=0.98,top=0.98,bottom=0.55,wspace=0.03)
    ax5 = plt.subplot(gs5[:,:])
    m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=-180,urcrnrlon=180,resolution='c',ax=ax5)
    m.drawcoastlines()
    m.drawmapboundary()
    cbar_ax = fig.add_axes([0.515, 0.49, 0.46, 0.05])

    #initialise plots with synthetic data
    x = np.arange(0.0, 3.0, 0.01)
    y = np.sin(2*np.pi*x)
    c = y*10.
    lats = np.random.randint(low=-90, high=90, size=len(x))
    lons = np.random.randint(low=-180, high=180, size=len(x))
    p1 = ax1.scatter(x,y,s=100,c=c,picker=1)
    p5 = ax5.scatter(lons,lats,s=50,c=c,picker=1)
    cb = fig.colorbar(p1, orientation = 'horizontal', cax=cbar_ax,label='')
    
    app = App(root)

    # a tk.DrawingArea
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.get_tk_widget().pack(side=BOTTOM, fill=BOTH, expand=1)
    canvas.show()
    toolbar = NavigationToolbar2TkAgg(canvas, root)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
    
    #update plots with actual data
    plot = update_data(app)
    #setup scatter and map callbacks
    callback = select_callbacks()
    #callback.setup_callbacks

    app.mainloop()
