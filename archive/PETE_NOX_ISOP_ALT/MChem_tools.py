# ------------  MChem Tools - tms -------------------------------------
# --------------  
# ---- Section 0 ----- Modules required

# I/O
import os
import glob
import csv
import sys
import pygchem
if pygchem.__version__ == '0.2.0':
    import pygchem.diagnostics as gdiag
else:
    # standard code or dev branch?
    try:
        from pygchem import datasets
    except:
        import pygchem.datafields as datasets


# Plotting/Analysis
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import datetime as datetime 
import numpy as np

# --------------- ------------- ------------- ------------- ------------- 
# ---- Section 1 ----- Core Programmes
# 1.01 - Open ctm.bpch  ( pygchem version 0.2.0)
# 1.02 - Get np array (4D) of ctm.bpch ( lon, lat, alt, time)
# 1.03 - Get Air mass 
# 1.04 - Plot up GEOS-Chem slice
# 1.06 - Process time/date to CV days equivilent - mje
# 1.07 - incremental increase datetime by given months - credit: Dave Webb
# 1.08 - processes provided files to extract data/names
# 1.09 - date specific (Year,Month,Day) planeflight output reader 
# 1.10 - Get var data for model run ( pygchem version > 0.3.0)

# --------------- ------------- ------------- ------------- ------------- 
# ---- Section 2 ----- Generic Tools
# 2.01 - find nearest values in array
# 2.02 - Get GEOS-Chem longitude index for a given latitude
# 2.03 - Get GEOS-Chem latitude index for a given latitude
# 2.04 - Get GC datetime
# 2.05 - Get GEOS-Chem altitude index for a given latitude

# --------------- ------------- ------------- ------------- ------------- 
# ---- Section 3 ----- Prod-Loss processing
# 3.01 - extract reactions tracked by prod loss diag for a given p/l family
# 3.02- Get tags for reactions
# 3.03 - Get prod loss reactions for a given family.
# 3.04 - Extract reactions to form a dictionary of active reactions

# --------------- ------------- ------------- ------------- ------------- 
# ---- Section 4 ----- Species and Reference GEOS-Chem data
# 4.01 - Class for holding extra Geos-Chem (GC) species data.
# 4.02 - GEOS-Chem species in latex form
# 4.03 - What GEOS-Chem (GC) Species am i? takes TRA_## & returns GC ID - setup for tms iodine tracer # > 53 (v9-01-03/v9-2 - iodine branch)
# 4.04 - Retrieve lat, lon and alt for a given resolution and region.
# 4.05 - Get CTM (GEOS-Chem) array dimension for a given resolution
# 4.06 - Convert gamap category/species name to Iris/bpch name
# 4.07 - Returns tracers unit and scale (if requested)
# 4.99 - Reference data (from GChem - credit: Gerrit Kuhlmann)

# --------------- ------------- ------------- ------------- ------------- 
# ---- Section 1 ----- Core Programmes

# --------------                                                                                 
# 1.01 - open ctm.bpch using PyGChem <= REDUNDENT
# --------------                                                                                
def open_ctm_bpch(wd, bpch_fname='ctm.bpch'):
    ctm_f = gdiag.CTMFile.fromfile(os.path.join(wd, bpch_fname))
    return ctm_f

# --------------                                                                                 
# 1.02 - get np array (4D) of ctm.bpch ( lon,lat , alt,time) <= REDUNDENT                         
# --------------                                                                                
def get_gc_data_np(ctm_f, species,category="IJ-AVG-$", debug=False):
    if (debug):
        print 'called get_np_gc_4D_diags'
    diagnostics = ctm_f.filter(name=species, category=category)
    for diag in diagnostics:
        scalar = (diag.values[:,:,:])[:,:,:,np.newaxis]
        if (debug):
            print diag.name ,'len(scalar)',len(scalar), 'type(scalar)' , \
                type(scalar) , 'diag.scale', diag.scale, 'scalar.shape', \
                scalar.shape,'diag.unit',diag.unit
        try:
            np_scalar = np.concatenate( (np_scalar,scalar), axis=3 )
        except NameError:
            np_scalar = scalar
        if (debug):
            print 'np_scalar' , type(np_scalar), len(np_scalar), \
                np_scalar.shape, 'scalar', type(scalar), len(scalar), \
                scalar.shape
    return np_scalar

# --------------
# 1.03 - get air mass (4D) numpy array <= REDUNDENT 
# -------------       
def get_air_mass_np(ctm_f, times=None, debug=False):
    if (debug):
        print 'called get air mass'

    diagnostics = ctm_f.filter(name='AD', category="BXHGHT-$",time=times)
    for diag in diagnostics:
        scalar = np.array( diag.values[:,:,:] )[:,:,:,np.newaxis]              # Grab data                                                                                                                          
        if (debug):
            print diag.name ,'len(scalar)',len(scalar), 'type(scalar)' , \
                type(scalar) , 'diag.scale', diag.scale, 'scalar.shape', \
                scalar.shape,'diag.unit',diag.unit
        try:
            np_scalar = np.concatenate( (np_scalar, scalar), axis=3 )
        except NameError:
            np_scalar = scalar
        if (debug):
            print 'np_scalar' , type(np_scalar), len(np_scalar), \
                np_scalar.shape, 'scalar', type(scalar), len(scalar),\
                scalar.shape
    return np_scalar

# -----                                                                                                            
# 1.04 - plot geos slice                                                                                           
# -----                                                                                                            
def map_plot(scalar, species=None, unit=None, res='4x5', **Kwargs):

    # Setup slices                                                                                                 
    # Grid/Mesh values for Lat, lon, & alt                                                                         
    lon, lat, alt = get_latlonalt4res( res=res )

    # Setup mesh grids                                                                                             
    x, y = np.meshgrid(lon,lat)
    print len(x), len(y)

    plt.ylabel('Latitude', fontsize = 20)
    plt.xlabel('Longitude',fontsize = 20)
    plt.title( '{} / {} '.format( latex_spec_name( species ), unit ) )

    # Setup map ("m") using Basemap                                                                                
    m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
                    llcrnrlon=-182.5,\
                    urcrnrlon=177.5,\
                    resolution='c')
    m.drawcoastlines()

    parallels = np.arange(-90,91,15)
    meridians = np.arange(-180,151,30)

    plt.xticks(meridians) # draw meridian lines
    plt.yticks(parallels) # draw parrelel lines
#    m.drawparallels(parallels) # add to map
                                                                                                                   
    # Create meshgrid to plot onto                                                                             
    x, y = np.meshgrid(*m(lon, lat))
    print len(x), len(y)

    plt.xlim(-180,175)
    plt.ylim(-89,89)

    poly = m.pcolor(lon, lat, scalar, cmap = plt.cm.Blues)

    # Add labels/annotations                                                                                       
    cb = plt.colorbar(poly, ax = m.ax,shrink=0.4)

    return plt , cb 

# --------------
# 1.06 - Process time/date to CV days equivilent - mje
# -------------
# translate year to "since2006" function
def year_to_since_2006(model):
    year=(model[:,0]//10000)
    month=((model[:,0]-year*10000)//100)
    day=(model[:,0]-year*10000-month*100)
    hour=model[:,1]//100
    min=(model[:,1]-hour*100)
    doy=[ datetime.datetime(np.int(year[i]),np.int(month[i]),np.int(day[i]),\
                                np.int(hour[i]),np.int(min[i]),0)- \
              datetime.datetime(2006,1,1,0,0,0) \
              for i in range(len(year))]
    since2006=[doy[i].days+doy[i].seconds/(24.*60.*60.) for i in range(len(doy))]
    return since2006


# -------------
# 1.07 -  incremental increase datetime by given months - credit: Dave Webb
# -------------
def add_months(sourcedate,months):
    month = sourcedate.month - 1 + months
    year = sourcedate.year + month / 12
    month = month % 12 + 1
    day = min(sourcedate.day,calendar.monthrange(year,month)[1])
    return datetime.datetime(year,month,day)
    
# --------------
# 1.08 - processes provided files to extract data/names
# -------------    
def process_files_to_read(files, location, big, names, debug=True):
    if (debug):
        print 'process_files_to_read called'
    print files
    reader=csv.reader(open(files,'rb'), delimiter=' ', skipinitialspace = True)
    for row in reader:
#        print location , (row[0] == 'POINT'), (row[1] == location) ,  len(row) , len(big), len(names)#, big.shape#, (row[2:])[0], (row[-1]) , l

#        print row
#        sys.exit( 0 )

        if row[1] == location: 
            new=row[2:]
            try:    
                big.append(new)
            except:
                big=[new]
        if row[0] == 'POINT':
            names = row[2:]
    return big, names


# --------------
# 1.09 - date specific (Year,Month,Day) planeflight output reader 
# -------------
def readfile(filename, location,  years_to_use, months_to_use, days_to_use, plot_all_data=False,debug=True, **kwargs):
    print 'readfile called'
    big, names = [],[]
    # sort for choosen years/months
    for files in filename:
             # loop through list on choosen years
        if (not plot_all_data):
            lll = 0
            for year in range(len(years_to_use)):
                if (("{0}".format(years_to_use[year])) in (files)) :
                    # is it not the last year specificied?
                    if (debug):
                        print 'years_to_use[year]', years_to_use[year], \
                            'years_to_use[-1]', years_to_use[-1]
                    if (not (years_to_use[year] == years_to_use[-1])):
                        # just read all years upto point uptyo final year
                        big, names=process_files_to_read(files, location,\
                            big, names)
                        print 'i got to line 91'
                    # If last year selected, then retrieve only 
                    # given months & days
                    if (years_to_use[year] == years_to_use[-1]):
                        # Plot months exceot last one
                        for month in range(len(months_to_use)):                                                                                                  
                            if (debug):
                                print 'months_to_use[month]', \
                                    months_to_use[month], 'months_to_use[-1]', \
                                    months_to_use[-1], 'months_to_use', \
                                    months_to_use, 'type(months_to_use)',\
                                     type(months_to_use)
                            if (("{0}{1}".format(years_to_use[year],\
                                months_to_use[month])) in files) :  

                                if (not (months_to_use[month] == \
                                     months_to_use[-1])):
                                    big, names=process_files_to_read(files,\
                                         location,big, names)
                                    print 'i got to line 100',  'month', month,\
                                     'in', len(months_to_use),  'year', year, \
                                     'in' , len(years_to_use)
                                if (months_to_use[month] == months_to_use[-1]):
                                    # For last month, plot days upto last day
                                    for day in range(len(days_to_use)):                                                                                          
                                        if (("{0}{1}{2}".format( \
                                            years_to_use[year], \
                                            months_to_use[month], \
                                            days_to_use[day])) in files) : 
                                            if (debug):
                                                print 'days_to_use[day]', \
                                                days_to_use[day], \
                                                'days_to_use[-1]', \
                                                days_to_use[-1]
                                            big, names = process_files_to_read(\
                                                files, location,big, names)
                                            if (debug):
                                                print 'i got to line 108'
                                                print 'readfile read big of'+\
                                                    ' size: ', len(big)
                                                
        if (plot_all_data):
            big, names=process_files_to_read(files, location,big, names)
            print 'reading all data'

    big=np.float64(big)             
    print 'readfile read big of size: ', len(big)
    if len(big) == 0:
        print 'ERROR - no output read in for:', location,  years_to_use, \
            months_to_use, days_to_use
        sys.exit( 0 )
        
    return big, names


# ----
# 1.10 - Get var data for model run
# ---
def get_GC_output( wd, vars=None, species=None, category=None, 
            r_cubes=False, r_res=False, restore_zero_scaling=True, 
            trop_limit=False, debug=False):
    """
        Data extractor for GEOS-Chem using pygchem (>= 0.3.0 )
        ( Credit: Ben Bovy -  https://github.com/benbovy/PyGChem )

        ctm.bpch files are extracted for a given directory, with only 
        the specific category and species returned. This can either be as
        a Iris Cube (retaining metadata) or as a numpy array.

            - To return Iris Cubes, set r_cubes to True
            - To get resolution of model run, set r_res to True

        For simplicity use variable names (species, category) as used in 
        the iris cube ( e.g. IJ_AVG_S__CO). if variable is unknown, just 
        print full dataset extracted to screen to see active diagnostics. 

        Species and category variables are maintained ( and translated ) 
        to allow for backwards compatibility with functions written for 
        pygchem version 0.2.0
        
        This replaces the now defunct functions: open_ctm_bpch 
        and get_gc_data_np.

        Examples and use of pygchem is discussed on Ben Bovy's GITHib 
        https://github.com/benbovy/PyGChem_examples/blob/master/Tutorials/Datasets.ipynb        
    """
    
    if debug:
        print 'Opening >{}<, for var: >{}<'.format( wd, ','.join(vars) ) +   \
            '(additional) gamap variables provided: >{}< + >{}<'.format( \
            category, species )

    # Error catch for older pygchem versions
    if pygchem.__version__ == '0.2.0':
        print 'WARNING: using {}, to extact ctm.bpch files'.format( pygchem.__version__) + \
            'using this dist. please use funcs: open_ctm_bpch and get_gc_data_np' 
        sys.exit( 0 )
                
    # Temporary fix for back compatibility: 
    # Convert gamap names ( species  + category ) to iris cube names
    # Just for use whilst updating functions written to use pygchem 0.2.0
    if not any( [ isinstance(i, type(None) ) for i in species, category ] ):
        # convert to Iris Cube name
        if (category == None) and ( vars ==  None ):
            category = "IJ-AVG-$"
        if (species == None) and ( vars ==  None ):
            species =  'O3'
    
        # remove scaling for 'IJ-AVG-$' - KLUDGE - needed for all ?
        if category == 'IJ-AVG-$':
            category = diagnosticname_gamap2iris( category )

        # setup var for Iris Cube, then remove known rm. chars.
        var = category+'__'+species
        vars=  [ var.replace('-', '_').replace('$', 'S') ]

        pass
    else:
        # Get default settings for reader
        if isinstance(vars, type(None)):
            vars = [ 'IJ_AVG_S__O3'  ]

    # ensure wd has a leading '/'
    if wd[-1] != '/':
        wd +=  '/'

    # Get files in dir ( more than one? )
    fns = sorted( glob.glob( wd+ '*ctm*' ) )
    if debug:
        print fns

    # Load files into Iris Cube
    cubes = datasets.load( fns, vars )

    # If no data extracted, print our variables
    try:
        [ i[0].data for i in cubes ]
    except:
        print datasets.load( fns[0] )
        print 'WARNING: no vars found for >{}<'.format( ','.join(vars) )
        sys.exit( 0 )
        
    # Temporary fix for back compatibility: 
    # Just extract as numpy 

    if not r_cubes:

        # Extract data
        arr = [ cubes[i].data for i in range( len(vars) ) ]

        # Limit to GEOS-Chem "chemical troposphere'
        if trop_limit:
            arr = [ i[...,:38] for i in arr]

        # convert to GC standard 4D fmt. - lon, lat, alt, time 
        if len((arr[0].shape)) == 4:
            if debug:
                print [i.shape for i in arr]
            arr = [np.rollaxis(i,0, 4) for i in arr]
            if debug:
                print [i.shape for i in arr]

        # For multiple vars, concatenate to var, lon, lat, lat, time 
        if len(vars) >1:
            arr = np.concatenate( [ i[None,...] for i in arr ], axis=0 )
        else:
            arr =arr[0]
            
        # Add altitude dimension to 2D (lon, lat)
        # a bug might occur for emission etc ( where lon, lat, time are dims )
        if len((arr.shape)) == 2:
            arr =arr[...,None]

    # Get res by comparing 1st 3 dims. against dict of GC dims.
    if r_res:
        res=get_dims4res( r_dims=True )[arr.shape[:3]]

    # Sort output - return Cubes?
    if r_cubes:
        output =  cube
    else:
        output = arr
#        return cube.data

    # Return model resolution? 
    if r_res:
        return output, res
    else:
        return output


# --------------- ------------- ------------- ------------- ------------- 
# ---- Section 2 ----- Generic Tools
# --------   
# 2.01 - Find nearest
# --------
def find_nearest(array,value):
    """
        Find nearest point. Adapted from HappyLeapSecond's Stackoverflow 
        answer. http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
    """
    idx = (np.abs(array-value)).argmin()
    return idx

# ----
# 2.02 - Get Longitube in GC grid box number
# ---- 
def get_gc_lon(lon, res='4x5',debug=False):
    lon_c, NIU, NIU = get_latlonalt4res( res=res )
    del NIU
    return find_nearest( lon_c, lon )

# ----
# 2.03 - Get Latitude in GC grid box number
# ---- 
def get_gc_lat(lat, res='4x5',debug=False):
    NIU, lat_c, NIU = get_latlonalt4res(res=res)
    del NIU
    return find_nearest( lat_c, lat )

# ----
# 2.04 - Get gc datetime
# -----
def get_gc_datetime(ctm_f, spec='O3', cat='IJ-AVG-$', debug=False):
    d  = ctm_f.filter(name=spec, category=cat)
    if debug:
        print '>'*30,d
        print '>.'*15, sorted(d)
    return [ i.times for i in d ]

# -------------- 
# 2.05 - get_gc_alt  ( KM => box num )   
# -------------
def get_gc_alt(alt):
    alt_c = gchemgrid('c_km_geos5_r')
    return find_nearest( alt_c, alt )


# --------------- ------------- ------------- ------------- ------------- 
# ---- Section 3 ----- Prod-Loss processing

# -------------
# 3.01 - extract reactions tracked by prod loss diag for a given p/l family
# ------------- 
def rxns_in_pl( wd, spec='LOX' ):
    fn =  'smv2.log'
    file_ =  open( wd+'/'+fn, 'rb' )
    readrxn  = False
    for row in file_:
        row = row.split()
        if all( [ i in row for i in 'Family','coefficient' , 'rxns', spec] ):
            readrxn=True
        if len(row) < 1 :
            readrxn=False
        if  readrxn:
            try:
                rxns.append( row )
            except:
                rxns = [ row ]          

    # -- remove 'Family' 
    rxns = [ i for i in rxns if (  'Family' not in i ) ]
    n = [int(rxn[1]) for rxn in rxns ]
    rxns = [rxn[2:] for rxn in rxns ]

    rdict = dict( zip(n, rxns) )
    return rdict

# -------------
# 3.02 - Get tags for reactions
# ------------- 
def get_p_l_tags( rxns ):

    # (PD??, RD??, LO3_??, PO3_??, LR??)
    for rxn in rxns:
#        print rxn
#        print [ i for i in rxn if  any( [x in i for x in 'PD', 'RD', 'PO3','LO3' , 'LR' ]) ]
        tags =  [i for i in rxn if any( [x in i for x in 'PD', 'RD', 'PO3','LO3' , 'LR' ]) ]
        try:
            tagsl.append( tags)
        except:
            tagsl = [tags]

    return tagsl

# -------------
# 3.03 - Get prod loss reactions for a given family.
# ------------- 
def prod_loss_4_spec( wd, spec ):
    # ---  Get Dict of all reactions, Keys = #s
    rdict = rxn_dict_from_smvlog(wd)

    # ---  Get reaction # tracked by p/l diag for spec and coefficient.
    rxns = rxns_in_pl(wd, spec)
    nums =  rxns.keys() 
    Coe = [ rxn[-1] for rxn in rxns.values() ]

    # --- get all details from full reaction dictionary
    rxns =  [ rdict[i] for i in nums ]

    # --- get tags for tracked reactions, state where reactions are un tracked
    tags = get_p_l_tags( rxns )

    return nums, rxns, tags, Coe

# -------------
# 3.04 - extract reactions to form a dictionary of active reactions
# ------------- 
def rxn_dict_from_smvlog( wd, spec='LOX' ):
    fn =  'smv2.log'
    file_ =  open( wd+'/'+fn, 'rb' )
    readrxn  = False
    for row in file_:
        row = row.split()
        if 'NMBR' in row:
            readrxn=True
        if len(row) < 1 :
            readrxn=False
        if  readrxn:
            try:
                rxns.append( row )
            except:
                rxns = [ row ]          

    # -- remove 'NMBR'
    rxns = [ i for i in rxns if (  'NMBR' not in i ) ]
    n = [int(rxn[0]) for rxn in rxns ]
    rxns = [rxn[1:] for rxn in rxns ]
    rdict = dict( zip(n, rxns) )
    return rdict
    


# --------------- ------------- ------------- ------------- ------------- 
# ---- Section 4 ----- Species and Reference GEOS-Chem data

# ------------
# 4.01 - Class for holding extra Geos-Chem (GC) species data. Input the string into the class init function and it can give several forms of data out.
# e.g. 
#  O3 = species('O3')
#  print O3.Latex
#  print O3.RMM
#  Output: 'O_3'   
# -------------

class species:
   def __init__(self, name):
      self.name = name
      self.help = ("""This is a class to get information on species from a local CSV folder
   It might contain the following information:
   self.RMM       = The Mean Mass of the species.
   self.latex     = The latex name of the species.
   self.smiles    = The smiles string of the species.
   self.InChI     = The InChI string of the species.
   """)
      species_filename = os.path.dirname(__file__) + "/Species.csv"

      try:
         species_file     = open(species_filename, 'rb')
      except IOError:
         print "Error: Species.csv does not appear to exist."
      species_csv      = csv.reader(species_file)

      if (name == 'OH'):
         self.group = 'CHEM-L=$'
      else:
         self.group = 'IJ-AVG-$'

      for row in species_csv:
         try:
            if (str(self.name) == row[0].strip()):
               self.formula   = row[1]
               self.InChI     = row[2]
               self.smiles    = row[3]
               self.RMM       = float(row[4])
               self.Latex     = row[5]
         except NameError:
            print "Species not found in CSV file"

# --------------
# 4.02 -  GEOS-Chem species in latex form
# --------------
def latex_spec_name(input_x, debug=False):
     spec_dict = {'OIO': 'OIO', 'C3H7I': 'C$_{3}$H$_{7}$I', 'IO': 'IO', 'I': 'I', 'I2': 'I$_{2}$', 'CH2ICl': 'CH$_{2}$ICl', 'HOI': 'HOI', 'CH2IBr': 'CH$_{2}$IBr', 'C2H5I': 'C$_{2}$H$_{5}$I', 'CH2I2': 'CH$_{2}$I$_{2}$', 'CH3IT': 'CH$_{3}$I', 'IONO': 'IONO','HIO3': 'HIO$_{3}$', 'ICl': 'ICl', 'I2O3': 'I$_{2}$O$_{3}$', 'I2O4': 'I$_{2}$O$_{4}$', 'I2O5': 'I$_{2}$O$_{5}$', 'INO': 'INO', 'I2O': 'I$_{2}$O', 'IBr': 'IBr','I2O2': 'I$_{2}$O$_{2}$', 'IONO2': 'IONO$_{2}$', 'HI':'HI', 'BrO':'BrO','Br':'Br','HOBr':'HOBr','Br2':'Br$_{2}$','CH3Br':'CH$_{3}$Br','CH2Br2':'CH$_{2}$Br$_{2}$', 'CHBr3':'CHBr$_{3}$','O3':'O$_{3}$', 'CO':'CO' , 'DMS':'DMS', 'NOx':'NOx', 'NO':'NO', 'NO2':'NO$_{2}$', 'NO3':'NO$_{3}$','HNO3':'HNO$_{3}$', 'HNO4':'HNO$_{4}$','PAN':'PAN', 'HNO2':'HNO$_{2}$', 'N2O5':'N$_{2}$O$_{5}$','ALK4':'>= C4 alkanes','ISOP':'CH$_{2}$=C(CH$_{3}$)CH=CH$_{2}$' ,'H2O2':'H$_{2}$O$_{2}$','ACET':'CH$_{3}$C(O)CH$_{3}$', 'MEK':'>C3 ketones', 'ALD2':'CH$_{3}$CHO', 'RCHO': 'CH$_{3}$CH$_{2}$CHO', 'MVK':'CH$_{2}$=CHC(O)CH$_{3}$', 'MACR':'CH$_{2}$=C(CH$_{3}$)CHO', 'PMN':'CH$_{2}$=C(CH$_{3}$)C(O)OONO$_{2}$', 'PPN':'CH$_{3}$CH$_{2}$C(O)OONO$_{2}$', 'R4N2':'>= C4 alkylnitrates','PRPE':'>= C4 alkenes', 'C3H8':'C$_{3}$H$_{8}$','CH2O':'CH$_{2}$O', 'C2H6':'C$_{2}$H$_{6}$', 'MP':'CH$_{3}$OOH', 'SO2':'SO$_{2}$', 'SO4':'SO$_{4}$','SO4s':'SO$_{4}$ on SSA', 'MSA':'CH$_{4}$SO$_{3}$','NH3':'NH$_{3}$', 'NH4': 'NH$_{4}$', 'NIT': 'InOrg N', 'NITs': 'InOrg N on SSA', 'BCPI':'BCPI', 'OCPI':'OCPI', 'BCPO':'BCPO','OCPO':'OCPO', 'DST1':'DST1', 'DST2':'DST2','DST3':'DST3','DST4':'DST4','SALA':'SALA', 'SALC':'SALC',  'HBr':'HBr', 'BrNO2': 'BrNO$_{2}$', 'BrNO3': 'BrNO$_{3}$', 'MPN':'CH$_{3}$ON$_{2}$', 'ISOPN':'ISOPN', 'MOBA':'MOBA', 'PROPNN':'PROPNN', 'HAC':'HAC', 'GLYC':'GLYC', 'MMN':'MMN', 'RIP':'RIP', 'IEPOX':'IEPOX','MAP':'MAP', 'AERI':'AERI' , 'Cl2':'Cl$_{2}$', 'Cl':'Cl','HOCl':'HOCl','ClO':'ClO','OClO':'OClO','BrCl':'BrCl' }

# --------------
# 4.03 - What GEOS-Chem (GC) Species am i? takes TRA_## & returns GC ID - setup for tms iodine tracer # > 53 (v9-01-03)
# -------------
def what_species_am_i(x) :
    TRA_lib={'O3':'O3','CO':'CO','NO':'NO','TRA_68': 'I2O', 'TRA_79': 'BrCl', 'TRA_71': 'I2O4', 'TRA_70': 'I2O3', 'TRA_59': 'IONO', 'TRA_58': 'HI', 'TRA_75': 'Cl', 'TRA_74': 'Cl2', 'TRA_77': 'ClO', 'TRA_76': 'HOCl', 'TRA_53': 'CH3Br', 'TRA_52': 'CH2Br2', 'TRA_51': 'CHBr3', 'TRA_50': 'BrNO3', 'TRA_57': 'OIO', 'TRA_56': 'IO', 'TRA_55': 'HOI', 'TRA_54': 'I2', 'TRA_69': 'INO', 'TRA_62': 'CH3I', 'TRA_63': 'CH2I2', 'TRA_60': 'IONO2', 'TRA_61': 'I2O2', 'TRA_48': 'HBr', 'TRA_49': 'BrNO2', 'TRA_64': 'IBr', 'TRA_65': 'ICl', 'TRA_44': 'Br2', 'TRA_45': 'Br', 'TRA_46': 'BrO', 'TRA_47': 'HOBr', 'TRA_73': 'AERI', 'TRA_72': 'I2O5', 'TRA_78': 'OClO', 'TRA_66': 'I', 'TRA_67': 'HIO3','CH3Br': 'REA_53', 'HOBr': 'REA_47', 'CHBr3': 'REA_51', 'Br2': 'REA_44', 'BrO': 'REA_46', 'Br': 'REA_45', 'CH2Br2': 'REA_52', 'BrNO2': 'REA_49', 'BrNO3': 'REA_50', 'HBr': 'REA_48','NO': 'NO', 'O3': 'O3', 'IO': 'TRA_56', 'CH3Br': 'TRA_53', 'HI': 'TRA_58', 'Br': 'TRA_45', 'IONO': 'TRA_59', 'Cl': 'TRA_75', 'BrO': 'TRA_46', 'HIO3': 'TRA_67', 'OClO': 'TRA_78', 'CH3I': 'TRA_62', 'CHBr3': 'TRA_51', 'ClO': 'TRA_77', 'I2O2': 'TRA_61', 'HOI': 'TRA_55', 'BrNO2': 'TRA_49', 'BrNO3': 'TRA_50', 'I': 'TRA_66', 'I2O': 'TRA_68', 'OIO': 'TRA_57', 'Cl2': 'TRA_74', 'BrCl': 'TRA_79', 'CH2Br2': 'TRA_52', 'ICl': 'TRA_65', 'CH2I2': 'TRA_63', 'IBr': 'TRA_64', 'I2O5': 'TRA_72', 'CO': 'CO', 'HBr': 'TRA_48', 'HOCl': 'TRA_76', 'HOBr': 'TRA_47', 'Br2': 'TRA_44', 'I2': 'TRA_54', 'I2O4': 'TRA_71', 'AERI': 'TRA_73', 'IONO2': 'TRA_60', 'I2O3': 'TRA_70', 'INO': 'TRA_69','REA_327': 'HO + CH3IT => H2O + I (CH2I)', 'REA_378': 'HI => .5I2', 'REA_325': 'IO + CH3O2 =(M)> I + HO2 + HCHO', 'REA_324': 'OIO + OH => HIO3', 'REA_373': 'I2O3 + I2O4 => 4AERI', 'REA_322': 'IO + IO (O2) =>I2O2', 'REA_363': 'I2O2 =(M)> IO + IO', 'REA_320': 'IO + IO (O2) => I + OIO', 'REA_309': 'I + O3 => IO + O2', 'REA_362': 'OIO + OIO => I2O4', 'REA_360': 'IO + OIO =(M)> I2O3', 'REA_367': 'I2O2 + O3 => I2O3 +O2', 'REA_369': 'I2O4 + O3 => I2O5 + O2', 'REA_365': 'I2O2 =(M)> OIO + I', 'REA_328': 'I + NO <=(M)> INO', 'REA_368': 'I2O3 + O3 => I2O4 +O2', 'REA_444': 'IONO =(hv)> I + NO2', 'REA_370': 'I2O4 =(M)> 2OIO', 'REA_446': 'I2O2 => IO + IO', 'REA_440': 'I2 =(hv)> 2I', 'REA_451': 'INO => I + NO', 'REA_447': 'CH3IT =>  <CH3 +I>', 'REA_441': 'HOI =(hv)> I + OH', 'REA_448': 'CH2I2 => <CH2 + I +I >', 'REA_332': 'I + NO2 <=(M)> IONO', 'REA_379': 'IONO2 => 0.5I2', 'REA_449': 'IBr =(hv)> I + Br', 'REA_445': 'IONO2 =(hv)> I + NO3', 'REA_375': 'I2O4 + I2O4 => 4AERI', 'REA_352': 'IO + BrO => Br + I + O2', 'REA_353': 'IO + BrO => Br +OIO', 'REA_350': 'IO + ClO => I + Cl + O2', 'REA_351': 'IO + ClO => I + OClO', 'REA_356': 'ICl + OH => HOCl  +I', 'REA_357': 'IBr + Br => I + Br2', 'REA_354': 'ICl + Cl => Cl2 + I', 'REA_355': 'ICL + Br => BrCl + I', 'REA_358': 'IBr + OH => HOI + Br', 'REA_359': 'IBr + OH => HOBr + I', 'REA_372': 'I2O3 + OIO => 3AERI', 'REA_334': 'IONO =(delta)> I + NO2', 'REA_335': 'IONO + IONO => I2 + 2NO2', 'REA_336': 'I2 + NO3 => I + IONO2', 'REA_337': 'IO + NO2 => IONO2', 'REA_330': 'INO =(delta)> NO + I', 'REA_331': 'INO + INO => I2 + 2NO', 'REA_318': 'IO + HO2 => HOI + O2', 'REA_319': 'IO + NO => I + NO2', 'REA_316': 'HI + OH => I + H2O', 'REA_317': 'HOI + OH => IO + H2O', 'REA_314': 'I + I2O => IO + I2', 'REA_315': 'I2 + OH => HOI + I', 'REA_374': 'I2O4 + OIO => 3AERI', 'REA_313': 'I + IO => I2O', 'REA_310': 'I + HO2 => HI + O2', 'REA_311': 'I + I =(M)> I2', 'REA_376': 'OIO + NO => NO2 + IO', 'REA_377': 'IO => .5I2', 'REA_442': 'IO =(hv)> I + O3', 'REA_450': 'ICl =(hv)> I + Cl', 'REA_339': 'IONO2 + I = I2 + NO3', 'REA_345': 'I2 + Br => I + IBr', 'REA_344': 'I2 + Cl => I + ICl', 'REA_347': 'IO + Cl => I + ClO', 'REA_346': 'I2 +BrO => IO + IBr', 'REA_340': 'IONO2 =(M)> IO + NO2', 'REA_343': 'I + BrO => IO + Br', 'REA_342': 'I + Br2 => Br + IBr', 'REA_381': '"HIO3 => AERI ( ""aerosol"")"', 'REA_380': '"OIO => AERI (""aerosol"")"', 'REA_383': 'HOI =>0.5I2', 'REA_382': 'I2O5 => 2AERI', 'REA_349': 'IO + ClO => ICl + O2', 'REA_348': 'IO + Br => I + BrO', 'REA_443': 'OIO =(hv)> I + O2(M)','PD59': 'IONO2 => 0.5I2', 'PD58': 'HI => .5I2', 'PD57': 'IO => .5I2', 'PD56': 'OIO + NO => NO2 + IO', 'PD55': 'I2O4 + I2O4 => 4AERI', 'PD54': 'I2O4 + OIO => 3AERI', 'PD53': 'I2O3 + I2O4 => 4AERI', 'PD52': 'I2O3 + OIO => 3AERI', 'PD51': 'I2O4 =(M)> 2OIO', 'PD50': 'I2O4 + O3 => I2O5 + O2', 'PD70': 'I2O2 => IO + IO', 'PD63': 'HOI =>0.5I2', 'PD69': 'IONO2 =(hv)> I + NO3', 'PD67': 'OIO =(hv)> I + O2(M)', 'PD62': 'I2O5 => 2AERI', 'PD71': 'CH3IT =>  <CH3 +I>', 'PD49': 'I2O3 + O3 => I2O4 +O2', 'PD39': 'ICL + Br => BrCl + I', 'PD38': 'ICl + Cl => Cl2 + I', 'PD60': '"OIO => AERI (""aerosol"")"', 'PD01': 'I + O3 => IO + O2', 'PD02': 'I + HO2 => HI + O2', 'PD03': 'I + I =(M)> I2', 'PD04': 'I + IO => I2O', 'PD05': 'I + I2O => IO + I2', 'PD06': 'I2 + OH => HOI + I', 'PD07': 'HI + OH => I + H2O', 'PD08': 'HOI + OH => IO + H2O', 'PD09': 'IO + HO2 => HOI + O2', 'PD24': 'IONO2 + I = I2 + NO3', 'PD25': 'IONO2 =(M)> IO + NO2', 'PD22': 'I2 + NO3 => I + IONO2', 'PD23': 'IO + NO2 => IONO2', 'PD20': 'IONO =(delta)> I + NO2', 'PD21': 'IONO + IONO => I2 + 2NO2', 'PD68': 'IONO =(hv)> I + NO2', 'PD47': 'I2O2 =(M)> OIO + I', 'PD48': 'I2O2 + O3 => I2O3 +O2', 'PD28': 'I2 + Cl => I + ICl', 'PD75': 'INO => I + NO', 'PD44': 'IO + OIO =(M)> I2O3', 'PD45': 'OIO + OIO => I2O4', 'PD46': 'I2O2 =(M)> IO + IO', 'PD29': 'I2 + Br => I + IBr', 'PD40': 'ICl + OH => HOCl  +I', 'PD41': 'IBr + Br => I + Br2', 'PD42': 'IBr + OH => HOI + Br', 'PD43': 'IBr + OH => HOBr + I', 'PD26': 'I + Br2 => Br + IBr', 'PD64': 'I2 =(hv)> 2I', 'PD27': 'I + BrO => IO + Br', 'PD65': 'HOI =(hv)> I + OH', 'PD74': 'ICl =(hv)> I + Cl', 'PD33': 'IO + ClO => ICl + O2', 'PD61': '"HIO3 => AERI ( ""aerosol"")"', 'PD72': 'CH2I2 => <CH2 + I +I >', 'PD32': 'IO + Br => I + BrO', 'PD66': 'IO =(hv)> I + O3', 'PD73': 'IBr =(hv)> I + Br', 'PD13': 'OIO + OH => HIO3', 'PD12': 'IO + IO (O2) =>I2O2', 'PD11': 'IO + IO (O2) => I + OIO', 'PD10': 'IO + NO => I + NO2', 'PD17': 'INO =(delta)> NO + I', 'PD16': 'I + NO <=(M)> INO', 'PD15': 'HO + CH3IT => H2O + I (CH2I)', 'PD14': 'IO + CH3O2 =(M)> I + HO2 + HCHO', 'PD31': 'IO + Cl => I + ClO', 'PD30': 'I2 +BrO => IO + IBr', 'PD19': 'I + NO2 <=(M)> IONO', 'PD18': 'INO + INO => I2 + 2NO', 'PD35': 'IO + ClO => I + OClO', 'PD34': 'IO + ClO => I + Cl + O2', 'PD37': 'IO + BrO => Br +OIO', 'PD36': 'IO + BrO => Br + I + O2','IONO2_hv': 'REA_445', 'ICl_hv': 'REA_450', 'CH2I2_hv': 'REA_448', 'IBr_hv': 'REA_449', 'IONO_hv': 'REA_444', 'OIO_hv': 'REA_443', 'IO_hv': 'REA_442', 'I2O2_hv': 'REA_446', 'INO_hv': 'REA_451', 'CH3IT_hv': 'REA_447', 'HOI_hv': 'REA_441', 'I2_hv': 'REA_440', 'NO2_hv': 'REA_385', 'BrNO2_hv': 'REA_438', 'NO3_hv_II': 'REA_394', 'BrNO3_hv_II': 'REA_437', 'O3_hv': 'REA_384', 'NO3_hv': 'REA_393', 'CH3Br_hv': 'REA_439', 'HOBr_hv': 'REA_435', 'HONO_hv': 'REA_391', 'Br2_hv': 'REA_433', 'BrNO3_hv': 'REA_436', 'BrO_hv': 'REA_434','REA_247': 'O3 emission', 'REA_391': 'HONO =(hv)>', 'REA_150': 'O3 + PRPE =>', 'REA_249': 'CH3Br emission', 'REA_152': 'O3 + PMN =>', 'REA_437': 'BrNO3 =(hv)>BrO +NO2', 'REA_436': 'BrNO3 =(hv)> Br + NO3', 'REA_435': 'HOBr =(hv)>', 'REA_434': 'BrO =(hv)>', 'REA_433': 'Br2 =(hv)>', 'REA_439': 'CH3Br =(hv)>', 'REA_438': 'BrNO2 =(hv)>Br + NO2', 'REA_394': 'NO3 =(hv)> ONO + O2', 'REA_393': 'NO3 =(hv)> NO2 + O3', 'REA_1': 'O3 + NO =>', 'REA_2': 'O3 + OH =>', 'REA_3': 'O3 +HO2 =>', 'REA_4': 'O3 + NO2 =>', 'REA_251': 'Br2 emission', 'REA_197': 'O3 + IALD =>', 'REA_169': 'O3 + MACR =>', 'REA_168': 'O3 + MVK =>', 'REA_254': 'I2 emission', 'REA_253': 'CH2I2 emission', 'REA_252': 'CH3IT emission', 'REA_277': 'Br + O3 =>', 'REA_250': 'CH2Br2 emission', 'REA_279': 'Br + OH2 =>', 'REA_278': 'Br + OH =>', 'REA_167': 'O3 + ISOP =>', 'REA_385': 'NO2 =(hv)>', 'REA_384': 'O3 =(hv)> OH + OH', 'REA_248': 'HNO3 emission', 'TRA_80':'CH2ICl', 'TRA_81': 'CH2IBr', 'TRA_83': 'C3H5I','TRA_82': 'C3H7I'}
    return TRA_lib[x]


# ----
#  4.04 - returns lon/lat/alt for res
# ----
def get_latlonalt4res(res='4x5', centre=True, hPa=False, nest=None, \
            debug=False):
    if debug:
        print res, centre, hPa, nest

    # Get dictionary variable name in Gerrit's GEOS-Chem dimensions list
    if centre:
        lon, lat = { 
        '4x5': ['c_lon_4x5' , 'c_lat_4x5'  ],
        '2x2.5':  ['c_lon_2x25' , 'c_lat_2x25' ], 
        '0.5x0.666' :[ 'c_lon_05x0667_EU', 'c_lat_05x0667_EU' ]  
        }[res]
    else:
        lon, lat = { 
        '4x5': ['e_lon_4x5' , 'e_lat_4x5' ],
        '2x2.5': ['e_lon_2x25' , 'e_lat_2x25' ], 
        '0.5x0.666' : ['e_lon_05x0667_EU', 'e_lat_05x0667_EU'] 
        }[res]
              
    # Override lat/lon variable name selection for nested grid (China)
    if nest == 'CH':
        if centre:
            lon, lat = 'c_lon_05x0667_CH', 'c_lat_05x0667_CH'
        else:
            lon, lat = 'e_lon_05x0667_CH', 'c_lat_05x0667_CH'
    if hPa:
        alt = 'c_hPa_geos5_r'
    else:
        alt='c_km_geos5_r'

    # Also provide high resolution grid if requested from this function all
    if nest =='high res global':
        lon, lat = np.arange(-180, 180, 0.25), np.arange(-90, 90, 0.25)
        return lon, lat, alt

    # Return values from gchemgrid
    d = gchemgrid(return_dict=True)
    if debug:
        print lon, lat, alt
    return  [ d[i] for i in lon, lat, alt ]

# --------   
# 4.05 - Get CTM (GEOS-Chem) array dimension for a given resolution
# --------
def get_dims4res(res=None, r_dims=False, invert=True, trop_limit=False):
    dims = {
    '4x5' :  (72,46,47), 
    '2x2.5':(144,91,47) , 
    '0.5x0.666':(121,81,47) 
    }
    if trop_limit:
        dims = list(dims)
        dims[-1]= 38
    if r_dims:
        if invert==True:
            return {v: k for k, v in dims.items()}
        else:
            return dims
    else:
        return dims[res ]

# --------   
# 4.06 - Convert gamap category/species name to Iris/bpch name
# --------
def diagnosticname_gamap2iris( x  ):
    d={
    "IJ-AVG-$": 'IJ_AVG_S'
    }
    return d[x]

# --------------
#  4.07 - Returns tracer's unit (and scale if requested)
# --------------
def tra_unit(x, scale=False, debug=False ):
    tra_unit = {
    'OCPI': 'ppbv', 'OCPO': 'ppbv', 'PPN': 'ppbv', 'HIO3': 'pptv', 'O3': 'ppbv', 'PAN': 'ppbv', 'ACET': 'ppbC', 'RIP': 'ppbv', 'BrNO3': 'pptv', 'Br': 'pptv', 'HBr': 'pptv', 'HAC': 'ppbv', 'ALD2': 'ppbC', 'HNO3': 'ppbv', 'HNO2': 'ppbv', 'C2H5I': 'pptv', 'HNO4': 'ppbv', 'OIO': 'pptv', 'MAP': 'ppbv', 'PRPE': 'ppbC', 'HI': 'pptv', 'CH2I2': 'pptv', 'IONO2': 'pptv', 'NIT': 'ppbv', 'CH3Br': 'pptv', 'C3H7I': 'pptv', 'C3H8': 'ppbC', 'DMS': 'ppbv', 'CH2O': 'ppbv', 'CH3IT': 'pptv', 'NO2': 'ppbv', 'NO3': 'ppbv', 'N2O5': 'ppbv', 'CHBr3': 'pptv', 'DST4': 'ppbv', 'DST3': 'ppbv', 'DST2': 'ppbv', 'DST1': 'ppbv', 'HOCl': 'ppbv', 'NITs': 'ppbv', 'RCHO': 'ppbv', 'C2H6': 'ppbC', 'MPN': 'ppbv', 'INO': 'pptv', 'MP': 'ppbv', 'CH2Br2': 'pptv', 'SALC': 'ppbv', 'NH3': 'ppbv', 'CH2ICl': 'pptv', 'IEPOX': 'ppbv', 'ClO': 'ppbv', 'NO': 'pptv', 'SALA': 'ppbv', 'MOBA': 'ppbv', 'R4N2': 'ppbv', 'BrCl': 'pptv', 'OClO': 'ppbv', 'PMN': 'ppbv', 'CO': 'ppbv', 'CH2IBr': 'pptv', 'ISOP': 'ppbC', 'BCPO': 'ppbv', 'MVK': 'ppbv', 'BrNO2': 'pptv', 'IONO': 'pptv', 'Cl2': 'ppbv', 'HOBr': 'pptv', 'PROPNN': 'ppbv', 'Cl': 'ppbv', 'I2O2': 'pptv', 'I2O3': 'pptv', 'I2O4': 'pptv', 'I2O5': 'pptv', 'MEK': 'ppbC', 'MMN': 'ppbv', 'ISOPN': 'ppbv', 'SO4s': 'ppbv', 'I2O': 'pptv', 'ALK4': 'ppbC', 'MSA': 'ppbv', 'I2': 'pptv', 'Br2': 'pptv', 'IBr': 'pptv', 'MACR': 'ppbv', 'I': 'pptv', 'AERI': 'pptv', 'HOI': 'pptv', 'BrO': 'pptv', 'NH4': 'ppbv', 'SO2': 'ppbv', 'SO4': 'ppbv', 'IO': 'pptv', 'H2O2': 'ppbv', 'BCPI': 'ppbv', 'ICl': 'pptv', 'GLYC': 'ppbv',
    # Extra diagnostics to allow for simplified processing 
    'CH3I':'pptv', 'Iy':'pptv', 'PSURF': 'hPa', 'OH':'pptv', 'HO2':'pptv', 'TSKIN':'K', 'GMAO_TEMP': 'K', 'GMAO_VWND' :'m/s','GMAO_UWND': 'm/s', 'RO2': 'pptv', 'U10M':'m/s','V10M': 'm/s' 
    } 
    units = tra_unit[x]

    if scale: 
        scaleby= get_unit_scaling( units )
    if scale:
        return units, scaleby
    else:
        return units

# --------   
# 4.14 - Get scaling for a given unit
# --------
def get_unit_scaling( units ):
    if any( [ (units ==  i) for i in 'pptv', 'pptC' ]):
        scaleby = 1E12
    if any( [ (units ==  i) for i in 'ppbv', 'ppbC' ]):
        scaleby = 1E9
    if any( [units ==i for i in  'K', 'm/s'  ] ):
        scaleby = 1
    return scaleby

# --------------
# 4.99 - Reference data, (inc. grid data) from GChem - credit: GK (Gerrit Kuhlmann )
# ------------- 

# --------------------
# ---  gchemgrid 
#---- credit: Gerrit Kuhlmann  
#! /usr/bin/env python
# coding: utf-8

# Python Script Collection for GEOS-Chem Chemistry Transport Model (gchem)
# Copyright (C) 2012 Gerrit Kuhlmann
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


#This module contains (some) grid coordinates used with GEOS-Chem as numpy
#arrays

def gchemgrid(input=None, return_dict=False, debug=False):
    d = {
    # 4x5
    'c_lon_4x5' : np.arange(-180, 175+5, 5) ,
    'e_lon_4x5' : np.arange(-182.5, 177.5+5, 5) ,
    'c_lat_4x5' : np.array( [-89]+ range(-86, 90, 4)+ [89] ),
    'e_lat_4x5' : np.array( [-90]+ range(-88, 92, 4)+ [90] ),
    
    # 2x2.5
    'e_lon_2x25' : np.arange(-181.25, 178.75+.5,2.5), 
    'c_lon_2x25' : np.arange(-180, 177.5+2.5, 2.5),
    'e_lat_2x25' :np.array( [-90]+range(-89, 91, 2)+[90] ) ,
    'c_lat_2x25' : np.array( [-89.5]+range(-88, 90, 2)+[89.5] ),

    # China nested grid
    'c_lon_05x0667_CH' : np.arange(70, 150+0.667, 0.667) ,
    'c_lat_05x0667_CH' : np.arange(-11, 55+.5, 0.5) ,
    'e_lon_05x0667_CH' : np.arange(70-(0.667/2), 150+0.667, 0.667) ,
    'e_lat_05x0667_CH' : np.arange(-11-(.5/2), 55+.5, 0.5) ,

    # Europe nested grid
    'c_lon_05x0667_EU'  : np.array([float(i) for i in np.arange(-30., 50+0.667, 0.667 ) ]) ,
    'c_lat_05x0667_EU' : np.arange(30, 70+.5, 0.5),
    'e_lon_05x0667_EU' : np.arange( -30.333, 50.333+0.666, 0.667),
    'e_lat_05x0667_EU' : np.arange(29.75, 70.25+0.25, 0.5 ),

    # generic arrays
    'e_lon_generic' : np.arange(-180.0,181.0),
    'c_lon_generic' : np.arange(-179.5,180.0),
    'e_lat_generic' : np.arange(-90.0,91.0),
    'c_lat_generic' : np.arange(-89.5,90.0),
    
    # Altitude - Grid box level edges (eta coordinate):
    'e_eta_geos5_r' : np.array([
            1.00179600e+00,   9.86769000e-01,   9.71665000e-01,
            9.56562000e-01,   9.41459000e-01,   9.26356000e-01,
            9.11253000e-01,   8.96152000e-01,   8.81051000e-01,
            8.65949000e-01,   8.50848000e-01,   8.35748000e-01,
            8.20648000e-01,   8.00515000e-01,   7.75350000e-01,
            7.50186000e-01,   7.25026000e-01,   6.99867000e-01,
            6.74708000e-01,   6.36974000e-01,   5.99251000e-01,
            5.61527000e-01,   5.23819000e-01,   4.86118000e-01,
            4.48431000e-01,   4.10759000e-01,   3.73114000e-01,
            3.35486000e-01,   2.85974000e-01,   2.42774000e-01,
            2.06167000e-01,   1.75170000e-01,   1.48896000e-01,
            1.26563000e-01,   1.07578000e-01,   9.14420000e-02,
            7.77260000e-02,   5.58200000e-02,   3.97680000e-02,
            2.80770000e-02,   1.95860000e-02,   9.19100000e-03,
            4.02600000e-03,   1.62500000e-03,   6.01000000e-04,
            1.99000000e-04,   5.50000000e-05,   0.00000000e+00]) ,

#Grid box level edges [km]:
    'e_km_geos5_r' : np.array([
            6.00000000e-03,   1.35000000e-01,   2.66000000e-01,
            3.99000000e-01,   5.33000000e-01,   6.69000000e-01,
            8.06000000e-01,   9.45000000e-01,   1.08600000e+00,
            1.22900000e+00,   1.37400000e+00,   1.52000000e+00,
            1.66900000e+00,   1.87100000e+00,   2.12800000e+00,
            2.39200000e+00,   2.66300000e+00,   2.94100000e+00,
            3.22800000e+00,   3.67300000e+00,   4.14000000e+00,
            4.63100000e+00,   5.14900000e+00,   5.69800000e+00,
            6.28300000e+00,   6.91000000e+00,   7.58700000e+00,
            8.32400000e+00,   9.41100000e+00,   1.05050000e+01,
            1.15780000e+01,   1.26330000e+01,   1.36740000e+01,
            1.47060000e+01,   1.57310000e+01,   1.67530000e+01,
            1.77730000e+01,   1.98550000e+01,   2.20040000e+01,
            2.42400000e+01,   2.65960000e+01,   3.17160000e+01,
            3.75740000e+01,   4.42860000e+01,   5.17880000e+01,
            5.99260000e+01,   6.83920000e+01,   8.05810000e+01]) ,

#Grid box level edges [hPa]:
    'e_hPa_geos5_r' : np.array([
            1.01181400e+03,   9.96636000e+02,   9.81382000e+02,
            9.66128000e+02,   9.50874000e+02,   9.35621000e+02,
            9.20367000e+02,   9.05114000e+02,   8.89862000e+02,
            8.74610000e+02,   8.59358000e+02,   8.44107000e+02,
            8.28856000e+02,   8.08522000e+02,   7.83106000e+02,
            7.57690000e+02,   7.32279000e+02,   7.06869000e+02,
            6.81458000e+02,   6.43348000e+02,   6.05247000e+02,
            5.67147000e+02,   5.29062000e+02,   4.90984000e+02,
            4.52921000e+02,   4.14873000e+02,   3.76851000e+02,
            3.38848000e+02,   2.88841000e+02,   2.45210000e+02,
            2.08236000e+02,   1.76930000e+02,   1.50393000e+02,
            1.27837000e+02,   1.08663000e+02,   9.23660000e+01,
            7.85120000e+01,   5.63880000e+01,   4.01750000e+01,
            2.83680000e+01,   1.97920000e+01,   9.29300000e+00,
            4.07700000e+00,   1.65100000e+00,   6.17000000e-01,
            2.11000000e-01,   6.60000000e-02,   1.00000000e-02]) ,

#Grid box level centers (eta-coordinates)
    'c_eta_geos5_r' : np.array([
            9.94283000e-01,   9.79217000e-01,   9.64113000e-01,
            9.49010000e-01,   9.33908000e-01,   9.18805000e-01,
            9.03703000e-01,   8.88601000e-01,   8.73500000e-01,
            8.58399000e-01,   8.43298000e-01,   8.28198000e-01,
            8.10582000e-01,   7.87933000e-01,   7.62768000e-01,
            7.37606000e-01,   7.12447000e-01,   6.87287000e-01,
            6.55841000e-01,   6.18113000e-01,   5.80389000e-01,
            5.42673000e-01,   5.04968000e-01,   4.67274000e-01,
            4.29595000e-01,   3.91937000e-01,   3.54300000e-01,
            3.10730000e-01,   2.64374000e-01,   2.24471000e-01,
            1.90668000e-01,   1.62033000e-01,   1.37729000e-01,
            1.17070000e-01,   9.95100000e-02,   8.45840000e-02,
            6.67730000e-02,   4.77940000e-02,   3.39230000e-02,
            2.38320000e-02,   1.43890000e-02,   6.60900000e-03,
            2.82500000e-03,   1.11300000e-03,   4.00000000e-04,
            1.27000000e-04,   2.80000000e-05]) , 

#Grid box level centers [km]
    'c_km_geos5_r' : np.array([
            7.10000000e-02,   2.01000000e-01,   3.32000000e-01,
            4.66000000e-01,   6.01000000e-01,   7.37000000e-01,
            8.75000000e-01,   1.01600000e+00,   1.15700000e+00,
            1.30100000e+00,   1.44700000e+00,   1.59400000e+00,
            1.76900000e+00,   1.99900000e+00,   2.25900000e+00,
            2.52700000e+00,   2.80100000e+00,   3.08400000e+00,
            3.44800000e+00,   3.90400000e+00,   4.38200000e+00,
            4.88600000e+00,   5.41900000e+00,   5.98500000e+00,
            6.59100000e+00,   7.24100000e+00,   7.94700000e+00,
            8.84800000e+00,   9.93800000e+00,   1.10210000e+01,
            1.20860000e+01,   1.31340000e+01,   1.41700000e+01,
            1.51980000e+01,   1.62220000e+01,   1.72430000e+01,
            1.87270000e+01,   2.08360000e+01,   2.30200000e+01,
            2.53070000e+01,   2.86540000e+01,   3.40240000e+01,
            4.01660000e+01,   4.71350000e+01,   5.48340000e+01,
            6.30540000e+01,   7.21800000e+01]) , 

#Grid box level centers [hPa]
    'c_hPa_geos5_r' : np.array([
            1.00422500e+03,   9.89009000e+02,   9.73755000e+02,
            9.58501000e+02,   9.43247000e+02,   9.27994000e+02,
            9.12741000e+02,   8.97488000e+02,   8.82236000e+02,
            8.66984000e+02,   8.51732000e+02,   8.36481000e+02,
            8.18689000e+02,   7.95814000e+02,   7.70398000e+02,
            7.44984000e+02,   7.19574000e+02,   6.94163000e+02,
            6.62403000e+02,   6.24298000e+02,   5.86197000e+02,
            5.48105000e+02,   5.10023000e+02,   4.71952000e+02,
            4.33897000e+02,   3.95862000e+02,   3.57850000e+02,
            3.13844000e+02,   2.67025000e+02,   2.26723000e+02,
            1.92583000e+02,   1.63661000e+02,   1.39115000e+02,
            1.18250000e+02,   1.00514000e+02,   8.54390000e+01,
            6.74500000e+01,   4.82820000e+01,   3.42720000e+01,
            2.40800000e+01,   1.45420000e+01,   6.68500000e+00,
            2.86400000e+00,   1.13400000e+00,   4.14000000e-01,
            1.39000000e-01,   3.80000000e-02]) , 
    }
    if (debug):
        print 'gchemgrid called'

    if return_dict:
        return d
    else:
        return d[input]

