import numpy as np

#FM = 'MINA'
#FM = 'HHRF'
FM = 'HMST'
label = 'DF4_2_HMST_20130708'
myPath = '/titan/cancerregulome9/ITMI_PTB/users/rtasseff/runMP'
dataPath = '/titan/cancerregulome9/ITMI_PTB/users/rtasseff/DF4/DF4_2/PPC'
tagPath = 'sampleMeta.dat' # can be commented out if FM != HMST


index = np.loadtxt(dataPath+'/index.tab',dtype=str)#[:10]

f=open('golemJobList_'+label+'.txt','w')

for i in range(len(index)):
	f.write('1 /tools/bin/python2.7 '+myPath+'/run_'+FM+'_FM_MP.py '+ \
		index[i,0]+' '+ \
		dataPath+'/'+index[i,5]+' '+ \
		index[i,4]+' '+ \
		dataPath+'/'+tagPath+'\n')
#		'\n')

f.close()
