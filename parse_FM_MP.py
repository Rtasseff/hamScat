import numpy as np

FM = 'MINA'
inDir = '/titan/cancerregulome9/ITMI_PTB/users/rtasseff/runMP'
outDir ='/titan/cancerregulome9/ITMI_PTB/users/rtasseff/DF4/DF4_2/PPC'
inName = '3193e0a3e9b3ec9645f7b35f15c09a7a.out.txt'
outName = 'results/DF4_2_test_20130628'
sampMetaName = 'sampleMeta.dat'
info = 'PPC'

pathIn = inDir+'/'+inName
pathOut = outDir+'/'+outName
sampMetaPath = outDir+'/'+sampMetaName
sampMeta = np.loadtxt(sampMetaPath,dtype=str)



memberList = ['M','F','NB']

foutFM = [open(pathOut+'_'+FM+'_'+memberList[0]+'_FM.tsv','w'), open(pathOut+'_'+FM+'_'+memberList[1]+'_FM.tsv','w'), open(pathOut+'_'+FM+'_'+memberList[2]+'_FM.tsv','w')]






j = 0
for member in memberList:
	keepInd = sampMeta[:,1]==member
	sampList = sampMeta[keepInd,0]
	foutFM[j].write('.')
	for samp in sampList:
		foutFM[j].write('\t'+samp)
	foutFM[j].write('\n')
	j = j+1






fin = open(pathIn)

for line in fin:
	data = line.rstrip().lstrip().split('\t')
	tsID = data[0]
	if tsID.find(FM) <0:
		raise ValueError('Feature type '+tsID+' found, expected '+FM'.')
	data = np.array(data[1:])
	j = 0
	for member in memberList:
		keepInd = sampMeta[:,1]==member	
		tmpData = data[keepInd]

		foutFM[j].write(tsID+':'+member+'::::'+info)
		for i in range(len(tmpData)):
			foutFM[j].write('\t'+tmpData[i])

		foutFM[j].write('\n')
		j = j+1
	
	
	

foutFM[0].close()
foutFM[1].close()
foutFM[2].close()


