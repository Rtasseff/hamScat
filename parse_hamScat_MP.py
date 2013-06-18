import numpy as np


inDir = '/titan/cancerregulome9/ITMI_PTB/users/rtasseff/runMP'
outDir ='/titan/cancerregulome9/ITMI_PTB/users/rtasseff/DF4/test_20130618'
name = 'test_HS_M_20130618.dat'
sampMetaName = 'sampleMeta.dat'

pathIn = inDir+'/'+name
pathOut = outDir+'/'+name
member = 'M'
sampMetaPath = outDir+'/'+sampMetaName
sampMeta = np.loadtxt(sampMetaPath,dtype=str)

keepInd = sampMeta[:,1]==member

sampMeta = sampMeta[keepInd]
fin = open(pathIn)

foutRep = open(pathOut+'_sepReport.tsv','w')
foutMS = open(pathOut+'_meanScat.tsv','w')
foutHHR = open(pathOut+'_HHRef.tsv','w')

foutRep.write('transcriptID\tsep_race\tsep_ptb\tsep_ptb_cat\tmeanScat1S_ptb\tmeanScat2S_ptb\n')

n = len(sampMeta)

foutHHR.write('.')
for i in range(n):
	foutHHR.write('\t'+sampMeta[i,0])
foutHHR.write('\n')


foutMS.write('.')
for i in range(n):
	foutMS.write('\t'+sampMeta[i,0])
foutMS.write('\n')


for line in fin:
	data = line.rstrip().lstrip().split('\t')
	tsID = data[0]
	foutRep.write(tsID)
	for i in range(1,6):
		foutRep.write('\t'+data[i])
	foutRep.write('\n')
	
	foutMS.write(tsID)
	for i in range(6,6+n):
		foutMS.write('\t'+data[i])
	foutMS.write('\n')
	
	foutHHR.write(tsID)
	for i in range(6+n,len(data)):
		foutHHR.write('\t'+data[i])
	foutHHR.write('\n')


foutMS.close()
foutRep.close()
foutHHR.close()


