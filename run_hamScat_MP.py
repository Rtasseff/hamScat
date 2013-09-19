# short script to run the hamming scatter analysis 
# RAT 20130617
# runing sevral diffrent tests on a given transcript
# lots of assumptions here
# last edited:
# RAT 20130816
# added easier feature selection
# more verbose output so not to mix up what results are
# doing each family member for each feature

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
# selected Feature path 
selFeatPath = sys.argv[5]

# load relevent tags
tag = np.loadtxt(tagPath,dtype=str,delimiter='\t')
# get selected tags
selFeat = np.loadtxt(selFeatPath,dtype=str,delimiter='\t')
# check for consistancy
tmp = selFeat[:,0]==tag[0,:]
if np.any(~tmp):
	raise ValueError('feature selections, '+selFeatPath+', do not match feature matrix, '+tagPath)
tag = tag[1:]

# get the correct tests features
featInd = np.arange(len(selFeat))[np.array(selFeat[:,1] == 'True')]

# setup output
out = '>>\t'+transID

# run on each family member
membList = ['F','M','NB']
for memb in membList:

	if nVar>1:

		# break out meta data into tags for given family member
		membInd = tag[:,1] == memb



		# get X from bianary matrix
		X = hamScat.readBiMat(matPath,nVar)
		#recall X is twice the size, 2 entries for each sample
		membIndX = np.zeros(2*len(membInd),dtype=bool)
		membIndX[::2]= membInd
		membIndX[1::2]= membInd

		X = X[:,membIndX]
		# get distances
		D = hamScat.calcNHDMat(X)
		
		# run tests
		for i in range(len(featInd)):
			ind = featInd[i]
			label = 'HSSEP|'+selFeat[ind,0]+'|'+memb
			p,sep = hamScat.empTestDS(D,tag[:,ind][tag[:,1] == memb],permList=[100,1E3,1E4,1E5,1E6])
			out = '%s\t%s\t%05.4E\t%05.4E'%(out,label,p,sep)


		# special one, diffrent algo, similar to tests on HMST
		#p1sPM,p2sPM = hamScat.rTestMS_cust1(D,tag[:,5])
		#out = '%s\t%s\t%05.4E\t%s\t%05.4E'%(out,'RMHMST1|B:CLIN:Preterm|'+memb,p1sPM,'RMHMST2|B:CLIN:Preterm|'+memb,p2sPM)

	else:
		for i in range(len(featInd)):
			ind = featInd[i]
			label = 'HSSEP|'+selFeat[ind,0]+'|'+memb
			out = out+'\t'+label+'\tnan\tnan'

		



# print to std out for golem to catch

# nan and inf get formated oddly 
out = out.replace('00INF','inf')
out = out.replace('00NAN','nan')

print out



