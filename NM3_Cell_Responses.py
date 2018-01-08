from allensdk.core.brain_observatory_cache import BrainObservatoryCache 
from allensdk.core.nwb_data_set import NwbDataSet 
from allensdk.brain_observatory.natural_movie import NaturalMovie 
from allensdk.brain_observatory.stimulus_analysis import StimulusAnalysis 
from allensdk.brain_observatory.dff import compute_dff 
from matplotlib import pyplot as plt 
from operator import add 
from matplotlib import style 
style.use("ggplot")
from sklearn.cluster import KMeans 

from mpl_toolkits.mplot3d import Axes3D

import allensdk.brain_observatory.stimulus_info as stim_info 
import matplotlib.pyplot as plt 
import matplotlib.patches as patches 
import pandas as pd
import pprint as p 
import numpy as np
import h5py



# _____________________________FUNCTIONS_________________________________________________________

def plot_stimulus_table(stim_table, title):
    fstart = stim_table.start.min()
    fend = stim_table.end.max()
    
    fig = plt.figure(figsize=(15,1))
    ax = fig.gca()
    for i, trial in stim_table.iterrows():    
        x1 = float(trial.start - fstart) / (fend - fstart)
        x2 = float(trial.end - fstart) / (fend - fstart)            
        ax.add_patch(patches.Rectangle((x1, 0.0), x2 - x1, 1.0, color='r'))
    ax.set_xticks((0,1))
    ax.set_xticklabels((fstart, fend))
    ax.set_yticks(())
    ax.set_title(title)
    ax.set_xlabel("frames")


# _________________________________________________________________________________________________

# this class uses the 'manifest' to keep track of downloaded metadata 
boc = BrainObservatoryCache(manifest_file='boc/manifest.json')

print '\n##############'
print ' Part I'
print '##############\n'

# 1: Download the experiments for a container 
experiment_containers = boc.get_experiment_containers()
# 2: Find a single experiment container 
exp_container = experiment_containers[0]
# 3. Retrieve the experiment container id 
exp_container_id = exp_container['id']
# 4: retrieve the three experiments associated with this experiment container
#    keep in mind that the experiments are contained in ophy_experiments
# 	 therefore, this line will search for the experiments assoc. with 'exp_container_id.' 
#    experiments = boc.get_ophys_experiments(experiment_container_ids=[exp_container_id])
# 5. one of the experiments should have a field value pair of 'session_type': 'three_session_A'
# 	 Find the experiment with the natural movie 3 stimlus 
#    therefore, rather than use #4, line 21, I can just use the line below: 
session_with_NM3_stim = boc.get_ophys_experiments(experiment_container_ids=[exp_container_id],
												stimuli=[stim_info.NATURAL_MOVIE_THREE]) # This is a list containing a dict 
# 6. For this specific experiment, I need to download the NWB file
#    this is done through retrieving the 'id' field of the 'session_with_NM3_stim'
session_id = session_with_NM3_stim[0]['id']
session_A_data = boc.get_ophys_experiment_data(session_id) # this is a dictionary 
# p.pprint(session_A_data.get_metadata()) # to view the metadata: 
# print session_A_data.list_stimuli() # to check the stimuli for this session:
# __________________________________________________________________________________
print '\n##############'
print ' Part II'
print '##############\n'

# 	 goal: With the session A data retrieved and stored relative to the manifest file, 
# 	 use the natural movie class to receive cell traces

## 	 Remember: Keep in mind that in part 1.6 we have data from all of session A, & not
##   just the information from the natural_movie_3 stimuli alone.
##   therefore, my first task is to get only the experiment data associated with 
##	 natural movie 3.  

# 1: read in the stimulus table, which describes when a given frame is displayed 
#    make a csv file to check out what the data looks like. &  
#    find when natural movie 3 is used in the stimulus set  
stim_table =  session_A_data.get_stimulus_table(stim_info.NATURAL_MOVIE_THREE)
stim_table.to_csv('stim_data', index = True, header = True)
frame_range = [0, 3599]
# plot_stimulus_table(stim_table, "frames %d -> %d " % (frame_range[0], frame_range[1]))
# plt.show()
# __________________________________________________________________________________
print '\n##############'
print ' Part III'
print '##############\n'

#	 goal: Figure out how to get the cell specimen data specifically for nm3 session A stimuli associated
#         with ophy_experiment_id: 570305847, experiment_container_id: 566759225

##	 Remember: you have downloaded the NWB data for an experiment. You need to get the fluorecent traces 
##   of the cells from this data

# 1: begin analysis using the 
# 	 natural_movie class. nm3 will be a NaturalMovie object.  
##   It will be important to ensure whether or not passing in the moview_name arguement
##   is necessary for data analysis of information 'I know' is recorded from the 
##   natural movie three stimulus. I believe the 'movie_name' param is neccessary 
##   because of the whole tuning analysis to a specific movie which is a big part of this class
nm3_object = NaturalMovie(session_A_data, stim_info.NATURAL_MOVIE_THREE )

#  2: grab the sweep response from the nm3_object and convert it to a pandas dataframe
sweep_response = nm3_object.get_sweep_response()
sweep_response = pd.DataFrame.from_records(sweep_response)

# print len(sweep_response.iloc[0][0])
# print sweep_response.iloc[1][0]
# __________________________________________________________________________________
print '\n##############'
print ' Part IV'
print '##############\n'

# 	 goal: using the sweep responce information, create traces of a cell for all 10 trials 
	
# 	 remember: the sweep_response dataframe is a trial(i) X cell matrix(j), 0 = i < 10, and 0 = j < 181.  
#    Also, Response(ij) in the matrix is an array. for each Response(ij), there are 3622 point in the array,
#    which represents 3622 frames, and is approx 2 minutes worth of data. Since there are 10 rows, for 10
#	 trials, we have a total of 19.921 minutes worth of data. 

def avg_sweep_response(cell_index, get_array = None):
#   1: find the average sweep response for a specific cell response over 10 trials 
#   each cell has it's own column that consists of 10 trials 
#   retrieve the data for each trial for cell_index at 0:  
    t1 = sweep_response.iloc[0][cell_index]
    t2 = sweep_response.iloc[1][cell_index]
    t3 = sweep_response.iloc[2][cell_index]
    t4 = sweep_response.iloc[3][cell_index]
    t5 = sweep_response.iloc[4][cell_index]
    t6 = sweep_response.iloc[5][cell_index]
    t7 = sweep_response.iloc[6][cell_index]
    t8 = sweep_response.iloc[7][cell_index]
    t9 = sweep_response.iloc[8][cell_index]
    t10 = sweep_response.iloc[9][cell_index]

    sum_array = [sum(x) for x in zip(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10)]
    avg_sweep_response_array = [x / 10 for x in sum_array]

    if get_array == True:
        return avg_sweep_response_array

    # 2: plot sweep response traces for this cell 
    stop = (3622 * .033) / 60 # conver 3622 frames to minutes  
    x = np.linspace(0, stop, num=3622)
    y = avg_sweep_response_array
    title = "Cell " + str(cell_index + 1) + ": Average Sweep Response"

    plt.figure(figsize=(20, 10))
    plt.title(title)
    plt.xlabel("Time (min)")
    plt.ylabel("avg DF/F")
    plt.plot(x, y)
    plt.show()

# cell_1_sweep_response = avg_sweep_response(0)
# cell_2_sweep_response = avg_sweep_response(1)

#_________________________________________________________________
print '\n##############'
print ' Part V'
print '##############\n'

# goal: We made it, now let's cluster this data Fool!!!!!!!!!!! :*) 
# combined the two arrays
cell_1_array = avg_sweep_response(0, True)
cell_2_array = avg_sweep_response(1, True)
cell_3_array = avg_sweep_response(2, True)


plt.scatter(cell_1_array, cell_2_array)
plt.show()

# for each element in list one and list two 
# create a new list which will contain a list of 
# the elements a specific index
x = []
for i in range(len(cell_1_array)):
    val_one = cell_1_array[i]
    val_two = cell_2_array[i]
    lists_for_x = [val_one, val_two]
    x.append(lists_for_x)

kmeans = KMeans(n_clusters = 2)
kmeans.fit(x)

#visualize the data 
centroids = kmeans.cluster_centers_
labels = kmeans.labels_ 

# print(centroids)
# print(labels[0:20])

colors = ["g.", "r.", "m."]

for i in range(len(x)):
    print("coordinate:", x[i], "label:", labels[i])
    plt.plot(x[i][0], x[i][1], colors[labels[i]], markersize = 10)

plt.scatter(centroids[:,0], centroids[:,1], marker = "*", s=150, linewidths = 5, zorder = 10)
plt.show()

#_____________3d visualization________________________

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
 
# x = []
# for i in range(len(cell_1_array)):
#     val_one = cell_1_array[i]
#     val_two = cell_2_array[i]
#     val_three = cell_3_array[i]
#     lists_for_x = [val_one, val_two, val_three]
#     x.append(lists_for_x)

# kmeans = KMeans(n_clusters = 3)
# kmeans.fit(x)

# #visualize the data 
# centroids = kmeans.cluster_centers_
# labels = kmeans.labels_ 

# # print(centroids)
# # print(labels[0:20])

# colors = ["g", "r.", "m."]


# ax.scatter(centroids[:,0], centroids[:,1], centroids[:,2], marker = "*", linewidths = 5, zorder = 10)
# plt.show()
