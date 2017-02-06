#!/usr/bin/python
# a Python script that uses pyKML to create a Hello World example
from lxml import etree
from pykml.factory import KML_ElementMaker as KML
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import numpy as np
from matplotlib.colors import ColorConverter 

data = [5.1,6.3,2,3,4]
time = [0,1,2,3,4]
lons  = [134,136,138,140,142]
lats = [0,1,2,1,0]
alts = [1000,2000,4000,5000,3000]

def rgb_to_hex(rgb):
    return '#%02x%02x%02x' % rgb

c = mcolors.ColorConverter().to_rgb

colors = data
plt.figure(figsize=(9, 1.5))
plt.scatter(time,alts, c=colors, cmap=plt.cm.coolwarm)
plt.gca().set_visible(False)
cax = plt.axes([0.1, 0.2, 0.8, 0.6])
plt.colorbar(orientation='horizontal',cax=cax)
#plt.colorbar(orientation='horizontal')
#plt.show()
plt.savefig("colourbar.png")

jet = cm = plt.get_cmap('coolwarm')
cNorm  = mcolors.Normalize(vmin=np.min(data), vmax=np.max(data))
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
rgba_vals = scalarMap.to_rgba(data)
print rgba_vals
rgb_vals = (rgba_vals[:,0:3])*255
rgb_vals = tuple(map(tuple, rgb_vals))
print rgb_vals
hex_points = []
for point in rgb_vals:
    hex_points = np.append(hex_points,rgb_to_hex(point))

hex_points = [i[1:] for i in hex_points]
print hex_points

#reorder hex for kml read - kml hex format is aabbggrr, a = alpha
r_hex = [i[0:2] for i in hex_points]
g_hex = [i[2:4] for i in hex_points]
b_hex = [i[4:] for i in hex_points]

bg_hex = [a+b for a,b in zip(b_hex,g_hex)]
print bg_hex
bgr_hex = [a+b for a,b in zip(bg_hex,r_hex)]
hex_points = bgr_hex

# create a document element with a single label style
kmlobj = KML.kml(
    KML.Document(
        KML.Style(
            KML.LabelStyle(
                KML.scale(6)
            ),
        )
    )
)

kmlobj.Document.append(
    KML.GroundOverlay(
        KML.name("Large-scale overlay on terrain"),
        KML.visibility("1"),
        KML.Icon(
          KML.href("http://gdata1.sci.gsfc.nasa.gov/daac-bin/G3/giovanni-wmx.cgi?SERVICE=WMS&WMTVER=1.0.0&REQUEST=GetMap&SRS=EPSG:4326&EXCEPTIONS=INIMAGE&LAYERS=MAMO_CHLO_4km.CR::l3m_data:max=30.0:min=0.05&STYLES=latlonplot:top=surface:bottom=surface:ctype=predefined:palette=Rainbow&BBOX=107.400,-24.000,177.600,17.400&TIME=2012-01-01T00:00:00Z/2012-02-01T00:00:00Z&FORMAT=GIF&transparent=true&WIDTH=800&HEIGHT=400"),
        ),
        KML.LatLonBox(
          KML.north("17.400"),
          KML.south("-24.000"),
          KML.east("177.600"),
          KML.west("107.400"),
        ),
    ),
)

kmlobj.Document.append(
    KML.ScreenOverlay(
        KML.name("Background colour bar"),
        KML.visibility("1"),
        KML.Icon(
          KML.href("l3m_data_MAMO_CHLO_4km_CR_AreaMap_2012-01_cb.png"),
        ),
        KML.overlayXY(x="0",y="1",xunits="fraction",yunits="fraction",),
        KML.screenXY(x="0",y="1",xunits="fraction",yunits="fraction",),
        KML.rotationXY(x="0",y="0",xunits="fraction",yunits="fraction",),
        KML.size(x="0.5",y="0.08",xunits="fraction",yunits="fraction",),
      )
)


kmlobj.Document.append(
    KML.ScreenOverlay(
        KML.name("Conc. Colour Bar"),
        KML.visibility("1"),
        KML.Icon(
          KML.href("colourbar.png"),
        ),
        KML.overlayXY(x="0",y="0.92",xunits="fraction",yunits="fraction",),
        KML.screenXY(x="0",y="0.92",xunits="fraction",yunits="fraction",),
        KML.rotationXY(x="0",y="0",xunits="fraction",yunits="fraction",),
        KML.size(x="0.5",y="0.08",xunits="fraction",yunits="fraction",),
      )
)

# add placemarks to the Document element
for i in range(len(data)):
    hex_color = 'ff%s'%(hex_points[i])
    print hex_color
    print data[i]
    kmlobj.Document.append(
        KML.Style(
            KML.IconStyle(
                KML.color(hex_color),
                KML.scale(1.0),
                KML.Icon(
                    KML.href("http://maps.google.com/mapfiles/kml/pal2/icon26.png"),
                ),
                id="mystyle"
            ),
            id="pushpin"
        ),
    )
    kmlobj.Document.append(
        KML.Placemark(
            KML.name(data[i]),
            KML.styleUrl("#pushpin"),
            KML.Point(
                KML.extrude(1),
                KML.altitudeMode('absolute'),
                KML.coordinates('{lon},{lat},{alt}'.format(
                        lon=lons[i],
                        lat=lats[i],
                        alt=alts[i],
                    ),
                ),
            ),
        )
    )
    print hex_color
    print data[i]



s_file = etree.tostring(etree.ElementTree(kmlobj),pretty_print=True)
text_file = open("python_kml.kml", "w")
text_file.write("%s"%s_file)
text_file.close()
