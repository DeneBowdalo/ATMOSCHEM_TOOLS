import numpy as np

start_year = 1988
end_year = 1993

typ = 'MH_obs_annual_magnitudes.npy'
big_vals = []

for count in range(21):

	filename = '%s_%s/%s'%(start_year,end_year,typ)

	vals = np.load(filename)
	#conv = 24/(2*np.pi)
	#vals = vals*conv	

	#print vals	
	big_vals = np.append(big_vals,vals)
	print big_vals

	start_year+=1
	end_year+=1

	np.save('annual_magnitudes',big_vals)
		
	
 
