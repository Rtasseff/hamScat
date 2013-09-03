import numpy as np

#FM = 'MINA'
#FM = 'HHRF'
FM = 'HMST'
label = 'DF4_2_miRNA_HMST_20130903'
myPath = '/titan/cancerregulome9/ITMI_PTB/users/rtasseff/runMP'
dataPath = '/titan/cancerregulome9/ITMI_PTB/users/rtasseff/DF4/DF4_2/miRNA'
tagPath = 'sampleMeta.dat' 


index = np.loadtxt(dataPath+'/index.tab',dtype=str)#[:10]

f=open('golemJobList_'+label+'.txt','w')

for i in range(len(index)):
	if int(index[i,4])>1:
		f.write('1 /tools/bin/python2.7 '+myPath+'/run_'+FM+'_FM_MP.py '+ \
			index[i,0]+' '+ \
			dataPath+'/'+index[i,5]+' '+ \
			index[i,4]+' '+ \
			dataPath+'/'+tagPath+'\n')


f.close()
