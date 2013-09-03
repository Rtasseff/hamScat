import numpy as np

FM = 'HMST'
#inDir = '/titan/cancerregulome9/ITMI_PTB/users/rtasseff/runMP'
inDir = '/titan/cancerregulome9/workspaces/golems/master'
outDir ='/titan/cancerregulome9/ITMI_PTB/users/rtasseff/DF4/DF4_2/miRNA'
inName = '71c9744ac681ec94233b57358d4dc25a.out.txt'
outName = 'results/DF4_2_miRNA_20130903'
sampMetaName = 'sampleMeta.dat'
info = 'miRNA'

pathIn = inDir+'/'+inName
pathOut = outDir+'/'+outName
sampMetaPath = outDir+'/'+sampMetaName
sampMeta = np.loadtxt(sampMetaPath,dtype=str)[1:]



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
		raise ValueError('Feature type '+tsID+' found, expected '+FM+'.')
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


