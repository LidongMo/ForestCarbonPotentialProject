from os import listdir
from os.path import isfile, join
import numpy as np

fout=open("Data/ForestCoverScaler/20230126_Merged_Table_with_Forest_Cover_Scaler_with_1km_resolution.csv","w+")


sampled_data = listdir("Data/ForestCoverScaler/ExtractedSubsets/")
sampled_data = list(filter(lambda f: f.endswith('.csv'), sampled_data))

# first file:
for line in open("Data/ForestCoverScaler/ExtractedSubsets//"+sampled_data[0]):
    fout.write(line)

# now the rest:
for num in range(1,len(sampled_data)):
	try:
		f = open("Data/ForestCoverScaler/ExtractedSubsets//"+sampled_data[num])
		f.__next__()
		for line in f:
			fout.write(line)
		f.close()
	except IOError:
	    pass

fout.close()
