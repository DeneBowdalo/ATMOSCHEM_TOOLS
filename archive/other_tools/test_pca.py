import numpy as np
from netCDF4 import Dataset
from sklearn.decomposition import PCA
from sklearn.datasets import load_iris
from sklearn.preprocessing import StandardScaler
import numpy as np
import matplotlib.pyplot as plt

o3_g = Dataset('GEOS_CHEM_SURFACE_O3_2005_2010_v90103_2x2.5_GEOS5.nc')
co_g = Dataset('GEOS_CHEM_SURFACE_CO_2005_2010_v90103_2x2.5_GEOS5.nc')
no_g = Dataset('GEOS_CHEM_SURFACE_NO_2005_2010_v90103_2x2.5_GEOS5.nc')
no2_g = Dataset('GEOS_CHEM_SURFACE_NO2_2005_2010_v90103_2x2.5_GEOS5.nc')
temp_g = Dataset('GEOS_CHEM_SURFACE_GMAO_TEMP_2005_2010_v90103_2x2.5_GEOS5.nc')
isop_g = Dataset('GEOS_CHEM_SURFACE_ISOP_2005_2010_v90103_2x2.5_GEOS5.nc')
ro2_g = Dataset('GEOS_CHEM_SURFACE_RO2_2005_2010_v90103_2x2.5_GEOS5.nc')

o3 = o3_g.variables['o3'][:]
co = co_g.variables['co'][:]
no = no_g.variables['no'][:]
no2 = no2_g.variables['no2'][:]
temp = temp_g.variables['gmao_temp'][:]
isop = isop_g.variables['isop'][:]
ro2 = ro2_g.variables['ro2'][:]
lat_c = o3_g.variables['lat_centre'][:]
lon_c = o3_g.variables['lon_centre'][:]

all_lats = []
for i in lat_c:
    all_lats = np.append(all_lats,[i]*len(lon_c))

print all_lats

print o3.shape

o3 = np.average(o3,axis=0)
co = np.average(co,axis=0)
no = np.average(no,axis=0)
no2 = np.average(no2,axis=0)
temp = np.average(temp,axis=0)
isop = np.average(isop,axis=0)
ro2 = np.average(ro2,axis=0)
print o3.shape

#o3 = o3[:1000]
#co = co[:1000]
#no = no[:1000]
#no2 = no2[:1000]

all_var = np.vstack((np.ravel(o3),np.ravel(co),np.ravel(no),np.ravel(no2),np.ravel(temp),np.ravel(isop)))

all_var = np.transpose(all_var)

print all_var.shape

def scikit_pca(X):

    # Standardize
    X_std = StandardScaler().fit_transform(X)

    o3_std = StandardScaler().fit_transform(o3)
    
    print 'standardised'

    # PCA
    sklearn_pca = PCA(n_components=2)

    print 'pca 1'    

    X_transf = sklearn_pca.fit_transform(X_std)

    print X_transf
    print X_transf.shape

    print 'pca 2'

    # Plot the data
    plt.scatter(X_transf[:,0],X_transf[:,1],c=ro2)
    plt.colorbar()
    plt.title('PCA via scikit-learn (using SVD)')
    plt.show()


scikit_pca(all_var)
