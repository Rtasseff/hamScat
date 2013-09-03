import numpy as np

label = 'DF4_2_miRNA_hamScat_20130830'
myPath = '/titan/cancerregulome9/ITMI_PTB/users/rtasseff/runMP'
dataPath = '/titan/cancerregulome9/ITMI_PTB/users/rtasseff/DF4/DF4_2/miRNA'
tagPath = dataPath+'/sampleMeta.dat'
#tmpInd = [2,4,7]
index = np.loadtxt(dataPath+'/index.tab',dtype=str)#[tmpInd]
selFeatPath = dataPath+'/sampleMeta_selectFeature.dat'

f=open('golemJobList_'+label+'.txt','w')

for i in range(len(index)):
	f.write('1 /tools/bin/python2.7 '+myPath+'/run_hamScat_MP.py '+ \
		index[i,0]+' '+ \
		dataPath+'/'+index[i,5]+' '+ \
		index[i,4]+' '+ \
		tagPath+' '+ \
		selFeatPath+'\n')

f.close()


