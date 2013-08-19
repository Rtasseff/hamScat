import numpy as np


inName = 'f75736c0e604deaf3b8f491bcf6e86f6.out.txt'
outName = 'results/DF4_2_additional_20130819'

inDir = '/titan/cancerregulome9/workspaces/golems/master'
outDir ='/titan/cancerregulome9/ITMI_PTB/users/rtasseff/DF4/DF4_2/PPC'
sampMetaName = 'sampleMeta.dat'
indexName = 'index.tab'


pathIn = inDir+'/'+inName
pathOut = outDir+'/'+outName+'_hamScat.tsv'
sampMetaPath = outDir+'/'+sampMetaName
sampMeta = np.loadtxt(sampMetaPath,dtype=str)
index = np.loadtxt(outDir+'/'+indexName,dtype=str)

# set up for ordering all runs completed
names = index[:,0]
del index
nNames = len(names)

# start reading the input file need the features
fin = open(pathIn)
line = fin.next()
data = line.strip().split('\t')
# check if the line for reading
if data[0]!='>>': raise ValueError('found an unexpeected line in '+inName+', line number = 0')
# get the expected labels 
labels=np.array(data[2::2],dtype=str)
m = len(labels)
# initialize the output matrix
repData = np.zeros((nNames,m)) - 1
# find location of data
tsID = data[1]
ind = names == tsID
repData[ind] = data[3::2]
# do this for all other data in input
count = 0
for line in fin:
	count += 1
	data = line.strip().split('\t')
	# check if the line for reading
	if data[0]!='>>': raise ValueError('found an unexpeected line in '+inName+', line number = '+str(count))
	# compare labels 
	if np.any(~(labels==np.array(data[2::2],dtype=str))): 
		raise ValueError('label missmatch in '+inName+', line number = '+str(count))
	# find location of data
	tsID = data[1]
	ind = names == tsID
	repData[ind] = data[3::2]
	

fin.close()

# lets write the report
fout = open(pathOut,'w')
# write the labels
fout.write('region_ID')
for i in range(m):
	fout.write('\t'+labels[i])
fout.write('\n')
# record missing data
foutMiss = open(pathOut+'_missing.dat','w')
for i in range(nNames):
	tsID = names[i]
	if np.any(repData[i] < 0):foutMiss.write(tsID+'\n')

	fout.write(tsID)
	for j in range(m):
		if repData[i,j] < 0: repData[i,j] = np.nan
		fout.write('\t'+str(repData[i,j]))
	fout.write('\n')


foutMiss.close()
fout.close()

