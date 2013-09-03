# short script to run the hamming scatter analysis 
# RAT 20130617
# runing sevral diffrent tests on a given transcript
# lots of assumptions here
# last edited:
# RAT 20130628

import sys
import numpy as np
import hamScat

# get the unique transcript ID
transID = sys.argv[1]
# the path to the transcript data matrix
matPath = sys.argv[2]
# number of variants in transcript data matrix
nVar = int(sys.argv[3])
# path to the tag infromation, sample meta data
tagPath = sys.argv[4]



# break out meta data into tags for given family member
tag = np.loadtxt(tagPath,dtype=str)

meanScat = np.zeros(len(tag))


for memb in ['M','F','NB']:
	membInd = tag[:,1] == memb


	# get X from bianary matrix
	X = hamScat.readBiMat(matPath,nVar)
	#recall X is twice the size, 2 entries for each sample
	membIndX = np.zeros(2*len(membInd),dtype=bool)
	membIndX[::2]= membInd
	membIndX[1::2]= membInd

	# get only the ones for our current member
	X = X[:,membIndX]

	# get distances
	D = hamScat.calcNHDMat(X)

	meanScat[membInd] = hamScat.calcMeanScat(D)





# print to std out for golem to catch 
out = 'N:HMST:'+transID
for i in range(len(meanScat)):
	out = '%s\t%05.4E' % (out,meanScat[i])
out = out+'\n'

print out



