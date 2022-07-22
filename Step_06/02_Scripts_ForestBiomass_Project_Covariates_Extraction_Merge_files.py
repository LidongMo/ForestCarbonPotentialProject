from os import listdir
from os.path import isfile, join
import numpy as np

fout=open("CovariatesTable/20201007_Merged_Covariates_sampled_dataset.csv","a")


sampled_data = listdir("CovariatesTable/CovariatesExtractedTables")
sampled_data = list(filter(lambda f: f.endswith('.csv'), sampled_data))

# first file:
for line in open("CovariatesTable/CovariatesExtractedTables/"+sampled_data[0]):
    fout.write(line)

# now the rest:
for num in range(1,len(sampled_data)):
	try:
		f = open("CovariatesTable/CovariatesExtractedTables/"+sampled_data[num])
		f.__next__()
		for line in f:
			fout.write(line)
		f.close()
	except IOError:
	    pass

fout.close()