import numpy as np

label = 'DF4_2_hamScat_F_20130712_missed'
myPath = '/titan/cancerregulome9/ITMI_PTB/users/rtasseff/runMP'
dataPath = '/titan/cancerregulome9/ITMI_PTB/users/rtasseff/DF4/DF4_2/PPC'
member = 'F'
tagPath = 'sampleMeta.dat'
index = np.loadtxt(dataPath+'/index.tab',dtype=str)#[:10]
missed = np.loadtxt(dataPath+'/results/DF4_2_20130712_scat_report_F_missing.dat',dtype=int)
index = index[missed]

f=open('golemJobList_'+label+'.txt','w')

for i in range(len(index)):
	f.write('1 /tools/bin/python2.7 '+myPath+'/run_hamScat_MP.py '+ \
		index[i,0]+' '+ \
		dataPath+'/'+index[i,5]+' '+ \
		index[i,4]+' '+ \
		dataPath+'/'+tagPath+' '+ \
		member+'\n')

f.close()


