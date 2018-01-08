# this assignment uses the Allen Brain Observatory SDK to retrieve experimental 
# recordings of experiment 511510667 cell repsonse to natural movie 3. This Data 
# will be clustered to form an inisghtful analysis

from allensdk.core.brain_observatory_cache import BrainObservatoryCache
from allensdk.brain_observatory.natural_movie import NaturalMovie
from allensdk.brain_observatory.drifting_gratings import DriftingGratings
import allensdk.brain_observatory.dff as dff 
import allensdk.brain_observatory.stimulus_info as stim_info 
import json 
import pandas as pd
import pprint as p
import matplotlib.pyplot as plt 
import numpy as np 
import h5py 


# PART I: FIND CELLS TO INVESTIGATE
# filter the experiments to only retrieve cell activity in VISp
## I need to find out what the filters all are
filter_json = """
[
	{
		"field": "area",
		"op": "in",
		"value": [
			"VISp"
		]
	}
]
"""

filters = json.loads(filter_json) 
# download the list of experiments 
boc = BrainObservatoryCache(manifest_file='boc/manifest.json') 
# retrieve the cell information and convert it to a dataframe 
cells = boc.get_cell_specimens(filters=filters)
cells = pd.DataFrame.from_records(cells)
print 'the number of cell data retrieved is: %d' % len(cells)

############# not needed ##############
# download the experiment containers for VISp experiment containers 
visp_experiment_containers = boc.get_experiment_containers(targeted_structures=['VISp'])
print 'the number of VISp experiment containers is: %d' % len(visp_experiment_containers)
# retrieve the VISp experiment ids from the 'visp_experiment_containers'
visp_experiment_ids = [ ec['id'] for ec in visp_experiment_containers]
visp_cells = cells[cells['experiment_container_id'].isin(visp_experiment_ids)]
# p.pprint(visp_cells)
# visp cells turned out to be the same as 'cells'
########################################

# find cells with cell metric reliability_nm3 > 0 natural movie stimuli 
nm3_cells = cells[cells['reliability_nm3'] > 0]
print 'the number of cells stimulated by natural movie 3 is: %d' % len(nm3_cells) 
# this gives 9152 rows 
# ___________________________________________________________________________________________________
# PART II: ANALYSIS 

# find the experiment containers for those cells 
nm3_Experiment_Container_ids = nm3_cells['experiment_container_id'].unique()
print 'total number of nm3 experiment containers is: %d' % len(nm3_Experiment_Container_ids)

# Download the ophys experiments containing the nm3 stimulus for VISp experiment containers 
nm3_Experiments = boc.get_ophys_experiments(experiment_container_ids=nm3_Experiment_Container_ids, stimuli=[stim_info.NATURAL_MOVIE_THREE])
print 'VISp natural movie 3 ophys experiments: %d' % len(nm3_Experiments)

# print '\nthis is the first of the of the natural movie 3 experiments:'
# p.pprint(nm3_Experiments[0])

# Download experiment data for a cell (retrieve NWB file)
nm3_cell = nm3_cells.iloc[0] #nm3_cell contains a dictionary of field : values for this cell 
print '\nexperiment data cell number:'
p.pprint(nm3_cell['cell_specimen_id'])

# Find which opyhs experiment has the Natural Movie Three files 
cell_experiment = boc.get_ophys_experiments(cell_specimen_ids=[nm3_cell['cell_specimen_id']], stimuli=[stim_info.NATURAL_MOVIE_THREE])[0]
print '\nhere is the ophys experiment that contains Natural Movie Three: '
p.pprint(cell_experiment)

data_set = boc.get_ophys_experiment_data(cell_experiment['id'])
print '\nhere is the metadata set for the cell experiment above:'
p.pprint(data_set.get_metadata())

# print '\nhere is the stimuli available for this file:'
# print(data_set.list_stimuli())
# ___________________________________________________________________________________________________
# PART III: Mess around with the Natural Movie class
print '\n******************'
print '     PART III     '
print '******************\n'

# access the NaturalMovie class
nm = NaturalMovie(data_set, stim_info.NATURAL_MOVIE_THREE)
# print 'Natural Movie class has the following characteristics'
# print nm.movie_name
# print nm._sweeplength
# print nm._sweep_response
# print str(nm.__dict__)
# print '----------------------------'
sweep_response = pd.DataFrame.from_records(nm.get_sweep_response()) #turned the numpy array into a pandas dataframe
# print sweep_response #the sweep response numpy array contains the dF/F response for each cell 
# print 'sweep length is: ' + str(nm.sweeplength)
# print ' lets try this peak thing again \n' + str(pd.DataFrame.from_records(nm.get_peak()))

# with sweep response, we can export that to a csv file to get a good look at it
sweep_response.to_csv('Sweep_Response_Data.csv',index = True, header = True)

# try to grab some of the information from sweep_response
# print 'here is the 0 col:'
# print sweep_response['0']

# print '\nhere is the 1 col'
# print sweep_response['1']

# print 'here is the 0 row:'
# print sweep_response.iloc[0]

# trial 1 information 
row0 = sweep_response.iloc[0]
# print ' here is row 0 values'
# print np.array(row0)

# ___________________________________________________________________________________________________
# PART IV.A: Compute dF/F
print '\n******************'
print '     PART IV.A    '
print '******************\n'
 
print 'uncomment section to see fluorescence trace'
print 'this section is only necessary to use if client is working with their own data' 
# cell_id = nm3_cell['cell_specimen_id']
# time, corrected_traces = data_set.get_corrected_fluorescence_traces(cell_specimen_ids=[cell_id])

# print 'here is what the numpy array of the corrected_traces look like:'
# print np.array(corrected_traces)
# # print len(corrected_traces[0])

# plt.figure(figsize=(14,4))
# plt.title("dF/F Trace")
# avg_fluorescence = dff.compute_dff(np.array(corrected_traces))
# plt.plot(time, avg_fluorescence[0,:])
# # plt.show()

# ___________________________________________________________________________________________________
# PART IV.B: Start Df/F Traces 
print '\n******************'
print '     PART IV.B    '
print '******************\n'

cell_id = nm3_cell['cell_specimen_id']
time, dff_traces = data_set.get_dff_traces(cell_specimen_ids=[cell_id])

plt.figure(figsize=(14,14))
plt.title("dF/F Trace")
plt.plot(time[:len(dff_traces[0])], dff_traces[0])

print 'here are the dff_traces' 