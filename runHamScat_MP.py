# short script to run the hamming scatter analysis 
# RAT 20130617
# runing sevral diffrent tests on a given transcript
# lots of assumptions here

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
# family member prefix
memb  = sys.argv[5]

# break out meta data into tags for given family member
tag = np.loadtxt(tagPath,dtype=str)
membInd = tag[:,1] == memb
tagR = np.array(tag[membInd,2],dtype=str)
tagPC = np.array(tag[membInd,3],dtype=str)
tagP = np.array(tag[membInd,4],dtype=str)


# get X from bianary matrix
X = hamScat.readBiMat(matPath,nVar)
#recall X is twice the size, 2 entries for each sample
membIndX = np.zeros(2*len(membInd),dtype=bool)
membIndX[::2]= membInd
membIndX[1::2]= membInd

X = X[:,membIndX]

# get the homo/hetro/ref score
Xmerg = X[:,::2] + X[:,1::2]
hhrScore = np.max(Xmerg,0)
hhrScore[hhrScore<3] = 0
hhrScore[hhrScore==3] = 1
hhrScore[hhrScore==4] = 2

# get distances
D = hamScat.calcNHDMat(X)

# Race test
pR = hamScat.empTestDS(D,tagR)

# PTB cat test 
pPC = hamScat.empTestDS(D,tagPC)

# PTB bi test, DS and MS( 1 and 2 sided)
pPD = hamScat.empTestDS(D,tagP)
p1sPM,p2sPM,meanScat = hamScat.rTestMS_cust(D,tagP)

# print to std out for golem to catch 
out = transID+'_'+memb+'\t'+str(pR)+'\t'+str(pPC)+'\t'+str(pPD)+'\t'+str(p1sPM)+'\t'+str(p2sPM)
for i in range(len(meanScat)):
	out = out+'\t'+str(meanScat[i])

for i in range(len(hhrScore)):
	out = out+'\t'+str(hhrScore[i])


out = out+'\n'

print out



