import plotly.plotly as py
from plotly.graph_objs import *
from plotly import tools
import numpy as np 
from mpl_toolkits.basemap import Basemap
import modules

# Create random data with numpy
import numpy as np

N = 1000
random_x = np.random.randn(N)
random_y = np.random.randn(N)

lat_c,lat_e,lon_c,lon_e = modules.model_grids('GEOSCHEM_4x5')

data = np.random.rand(len(lat_c),len(lon_c))

# Create a trace
trace1 = Heatmap(
    z=data,
    x=lon_e,
    y=lat_e,
    colorscale="RdBu",
    zauto=False,  # custom contour levels
    zmin=0,      # first contour level
    zmax=1        # last contour level  => colorscale is centered about 0
)


m = Basemap() 

# Make trace-generating function (return a Scatter object)
def make_scatter(x,y):
    return Scatter(
        x=x,
        y=y,
        mode='lines',
        line=Line(color="black"),
        name=' '  # no name on hover
    )

# Functions converting coastline/country polygons to lon/lat traces
def polygons_to_traces(poly_paths, N_poly):
    ''' 
    pos arg 1. (poly_paths): paths to polygons
    pos arg 2. (N_poly): number of polygon to convert
    '''
    traces = []  # init. plotting list 

    for i_poly in range(N_poly):
        poly_path = poly_paths[i_poly]
        
        # get the Basemap coordinates of each segment
        coords_cc = np.array(
            [(vertex[0],vertex[1]) 
             for (vertex,code) in poly_path.iter_segments(simplify=False)]
        )
        
        # convert coordinates to lon/lat by 'inverting' the Basemap projection
        lon_cc, lat_cc = m(coords_cc[:,0],coords_cc[:,1], inverse=True)
        
        # add plot.ly plotting options
        traces.append(make_scatter(lon_cc,lat_cc))
     
    return traces

# Function generating coastline lon/lat traces
def get_coastline_traces():
    poly_paths = m.drawcoastlines().get_paths() # coastline polygon paths
    N_poly = 91  # use only the 91st biggest coastlines (i.e. no rivers)
    return polygons_to_traces(poly_paths, N_poly)

# Function generating country lon/lat traces
def get_country_traces():
    poly_paths = m.drawcountries().get_paths() # country polygon paths
    N_poly = len(poly_paths)  # use all countries
    return polygons_to_traces(poly_paths, N_poly)

# Get list of of coastline and country lon/lat traces
traces_cc = get_coastline_traces()+get_country_traces()

#trace1 = Scatter(
#    x=[1, 2, 3],
#    y=[2, 1, 2])

trace2 = Scatter(
    x=[1, 2, 3],
    y=[2, 1, 2])

data = []
#data.append(trace1)
data.append(Data([trace1]+traces_cc))
data.append(trace2)

fig = tools.make_subplots(rows=2)

axis_style = dict(zeroline=False,showline=False,showgrid=False,ticks='',showticklabels=False)

layout = Layout(height=1000, width=1000, title='Stacked subplots',showlegend=False,hovermode="closest",xaxis=XAxis(axis_style),yaxis=YAxis(axis_style),autosize=False)

fig = { 'data':data, 'layout':layout }

url = py.plot(fig, filename='test_map',auto_open=False)
