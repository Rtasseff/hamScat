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
tag = np.loadtxt(tagPath,dtype=str)[1:]
m = len(tag)

out = 'N:MINA:'+transID

if nVar>1:

	# get X from bianary matrix
	X = hamScat.readBiMat(matPath,nVar)


	sumMA = hamScat.calcSumMA(X)


	# print to std out for golem to catch 
	out = 'N:MINA:'+transID
	for i in range(len(sumMA)):
		out = '%s\t%05.4E' % (out,sumMA[i]) 
	out = out+'\n'
else:
	for i in range(m):
		out = out+'\tnan'
	out = out+'\n'		

print out



