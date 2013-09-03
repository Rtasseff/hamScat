import numpy as np

label = 'DF4_2_hamScat_addFeat_2_missed_20130903_test'
myPath = '/titan/cancerregulome9/ITMI_PTB/users/rtasseff/runMP'
dataPath = '/titan/cancerregulome9/ITMI_PTB/users/rtasseff/DF4/DF4_2/PPC'
tagPath = dataPath+'/sampleMeta.dat'
index = np.loadtxt(dataPath+'/index.tab',dtype=str)
selFeatPath = myPath+'/sampleMeta_selectFeature_test.dat'
missing = np.loadtxt(dataPath+'/results/DF4_2_additional2_20130830_hamScat.tsv_missing.dat',dtype=str)

f=open('golemJobList_'+label+'.txt','w')

for i in range(len(index)):
	if np.any(index[i,0]==missing):
		f.write('1 /tools/bin/python2.7 '+myPath+'/run_hamScat_MP.py '+ \
			index[i,0]+' '+ \
			dataPath+'/'+index[i,5]+' '+ \
			index[i,4]+' '+ \
			tagPath+' '+ \
			selFeatPath+'\n')

f.close()


