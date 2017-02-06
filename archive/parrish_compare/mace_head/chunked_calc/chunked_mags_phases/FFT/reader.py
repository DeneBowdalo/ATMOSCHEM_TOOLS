import numpy as np

start_year = 1988
end_year = 1993

typ = 'MH_obs_annual_phases.npy'
big_vals = []

for count in range(21):

	filename = '%s_%s/%s'%(start_year,end_year,typ)

	values = np.load(filename)
	print values
	#convert phase from -pi:pi to 0:2pi
	for num in range(len(values)):
		if values[num] < 0:
			diff = np.pi- np.abs(values[num])
			values[num] = diff
		elif values[num] > 0:
			values[num] = np.pi + values[num]

	print values
	conv = 365.25/(2*np.pi)
	vals = values*conv	

	#print vals	
	big_vals = np.append(big_vals,vals)
	print big_vals

	start_year+=1
	end_year+=1

	np.save('annual_phases',big_vals)
		
	
 
