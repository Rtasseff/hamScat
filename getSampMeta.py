# Python 2.7 script for collecting data and puting it in the right format.
# This script is very specific to curent setup given the input paths and 
# formats that are expected, it is not meant to generalize, although the 
# ouput format may be usefull
# 
# RAT 20130816

import numpy as np

dataDir = '/titan/cancerregulome9/ITMI_PTB/users/rtasseff/DF4/clin'

targetDataPath = dataDir+'/targetphenotypes.tsv'
originalDataPath = dataDir+'/2013_07_12_clinical_ptb_race_seq'
race14DataPath = dataDir+'/df4_admix_14subpop.txt'

sampleIDDir = '/titan/cancerregulome9/ITMI_PTB/users/rtasseff/DF4/DF4_2/PPC'
sampleIDPath = sampleIDDir+'/sampleID.dat'

outputDir = sampleIDDir
outputName = 'sampleMeta'
outputPath = outputDir+'/'+outputName+'.dat'
outputSelectPath = outputDir+'/'+outputName+'_selectFeature.dat'


# get some of the family data in the original and 
# put the race values in
data = np.loadtxt(originalDataPath,delimiter='\t',dtype=str)
n,m = data.shape
fam2meta = {}
for i in range(1,m):
	fam = data[0,i][:7]
	dataBit = data[1:,i]
	race = [dataBit[3],dataBit[4],dataBit[5]]
	fam2meta[fam] = [race]

labels = ['C:ADMX:Admix_80_Percent']

# get the family wise target features
# be sure to parse out the member specific features
data = np.loadtxt(targetDataPath,delimiter='\t',dtype=str)
n,m = data.shape
for i in range(1,m):
	fam = data[0,i][:7]
	for j in range(1,14):
		fam2meta[fam].append(data[j,i])
		if i==1:labels.append(data[j,0].split(':')[0]+':'+data[j,0].split(':')[1]+':'+data[j,0].split(':')[2])
	# order, member wise F M NB
	fam2meta[fam].append(data[14:17,i])
	if i==1:labels.append('N:QCTL:ORDER')
	# version, member wise F M NB
	fam2meta[fam].append(data[17:20,i])
	if i==1:labels.append('C:QCTL:PL_VERSION')

	

# lets get the sample wise 14 pop measure
# we are going to do max here
data =  np.loadtxt(race14DataPath,delimiter='\t',dtype=str)
samp2race = {}
n,m = data.shape
races = data[0,1:]
for i in range(1,n):
	dataBit = np.array(data[i,1:],dtype=float)
	ind = np.argmax(dataBit)
	samp2race[data[i,0]]=races[ind]




# get the sample order from the VCF
samp = np.loadtxt(sampleIDPath,dtype = str)
n = len(samp)

# lets do the output
fout = open(outputPath,'w')


# write the labels
fout.write('Family_ID\tFamily_Member\tC:ADMX:Admix_14_max')
for i in range(len(labels)):
	fout.write('\t'+labels[i])
fout.write('\n')

# lets also create a run file to turn features on and off
foutTmp = open(outputSelectPath,'w')
foutTmp.write('Family_ID\tFalse\nFamily_Member\tFalse\nC:ADMX:Admix_14_max\tFalse')
for i in range(len(labels)):
	foutTmp.write('\n'+labels[i]+'\tFalse')
foutTmp.close()


# record features 
for i in range(n):
	# record the family ID
	fout.write(samp[i][:7])
	# get and record the member 
	memb = samp[i][8:]
	if len(memb)>2:
		memb = memb[:2]
	fout.write('\t'+memb)
	# the 14 cat race feature 
	fout.write('\t'+samp2race[samp[i]])
	# get the extra meta data and record
	meta = fam2meta[samp[i][:7]]
	for j in range(len(meta)):
		# some are family wise and need to be chosen based on member
		if (type(meta[j])==np.ndarray or type(meta[j])==list) and len(meta[j])==3:
			if memb=='F':data=meta[j][0]
			elif memb=='M':data=meta[j][1]
			elif memb=='NB':data=meta[j][2]
			else: raise ValueError('error 001 for samp '+samp[i]+' feature '+labels[j]+' -> '+str(meta[j]))
		elif type(meta[j])==np.string_ or type(meta[j])==str: data=meta[j]
		else: raise ValueError('error 002 for samp '+samp[i]+' feature '+labels[j]+' -> '+str(meta[j]))
		if data=='NaN' or data=='NA' or data=='na': data = 'nan'
		fout.write('\t'+data)
	fout.write('\n')
fout.close()




	
	
	
	


	


