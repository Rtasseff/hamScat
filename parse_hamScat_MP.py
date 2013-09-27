import numpy as np


inName = '7d94fc54cf1c4338d907ac67285df451.out.txt'
outName = 'results/tmp_HSSEP'

inDir = '/titan/cancerregulome9/workspaces/golems/master'
outDir ='/titan/cancerregulome9/ITMI_PTB/users/rtasseff/DF4/DF4_2/PPC'
sampMetaName = 'sampleMeta.dat'
indexName = 'index.tab'


pathIn = inDir+'/'+inName
pathOut = outDir+'/'+outName+'_hamScat_p.tsv'
pathOutSep = outDir+'/'+outName+'_hamScat_sep.tsv'
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
# check if the line is for reading
if data[0]!='>>': raise ValueError('found an unexpeected line in '+inName+', line number = 0')
# get the expected labels 
labels=np.array(data[2::3],dtype=str)
m = len(labels)
# initialize the output matrix
repDataP = np.zeros((nNames,m)) - 1
repDataS = np.zeros((nNames,m))
# find location of data
tsID = data[1]
ind = names == tsID
repDataP[ind] = data[3::3]
repDataS[ind] = data[4::3]
# do this for all other data in input
count = 0
for line in fin:
	count += 1
	data = line.strip().split('\t')
	# check if the line for reading
	if data[0]!='>>': raise ValueError('found an unexpeected line in '+inName+', line number = '+str(count))
	# compare labels 
	if np.any(~(labels==np.array(data[2::3],dtype=str))): 
		raise ValueError('label missmatch in '+inName+', line number = '+str(count))
	# find location of data
	tsID = data[1]
	ind = names == tsID
	repDataP[ind] = data[3::3]
	repDataS[ind] = data[4::3]
	

fin.close()

# lets write the report
fout = open(pathOut,'w')
# doing sep in spe file, if we start adding more info we should combine into one file
foutSep = open(pathOutSep,'w')
# write the labels
fout.write('region_ID')
foutSep.write('region_ID')

for i in range(m):
	fout.write('\t'+labels[i])
	foutSep.write('\t'+labels[i])
fout.write('\n')
foutSep.write('\n')
# record missing data
foutMiss = open(pathOut+'_missing.dat','w')
for i in range(nNames):
	tsID = names[i]
	if np.any(repDataP[i] < 0):foutMiss.write(tsID+'\n')

	out = tsID
	outSep = tsID
	for j in range(m):
		if repDataP[i,j] < 0 or np.isnan(repDataP[i,j]) or np.isnan(repDataS[i,j]): 
			out = out+'\tnan'
			outSep = outSep+'\tnan'
		else: 
			out = '%s\t%05.4E' % (out,repDataP[i,j])
			if repDataS[i,j]==np.inf:
				outSep = outSep+'\tinf'
			else:
				outSep = '%s\t%05.4E' % (outSep,repDataS[i,j])
	fout.write(out+'\n')
	foutSep.write(outSep+'\n')


foutMiss.close()
fout.close()
foutSep.close()
