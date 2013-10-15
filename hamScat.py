import numpy as np
import statsUtil
from gpdPerm import gpdPerm
import struct

eps = 1E-26




### hamming distances ###

def _hamDist(gen_i,gen_j):
	dist = 0
	if gen_i[0]!=gen_j[0]:dist = dist+1
	if gen_i[1]!=gen_j[1]:dist = dist+1
	return(dist)

def calcHD(gen_i,gen_j):
	""" find genetic distance between 2 varriants.
	"""
	gen_i1 = gen_i[0]
	gen_i2 = gen_i[1]
	gen_j1 = gen_j[0]
	gen_j2 = gen_j[1]

	dist = _hamDist((gen_i1,gen_i2),(gen_j1,gen_j2))
	if gen_i1!=gen_i2:
		tmp = _hamDist((gen_i2,gen_i1),(gen_j1,gen_j2))
		if tmp < dist: dist = tmp
	
	return(dist)


# calc the normalized hamming dist matrix
def calcNHDMat(X,countMiss=True):
	"""calculat the the normalized hamming distance matrix
	This is specific to the genotype matrices.
	Currently we assume that each row is a variant
	each column is 1/2 a sample ie two calls for each sample
	so every 2 col is a sample.
	Assuming integers
	Missing data is a zero 
	Comparisions with any missing data are ignored.
	coutMiss	if true - count missing in normalization,
			but unable to calc dif so equivlant to
			setting as no observed dif
			else - ignore and do not calc dif and
			do not count in normalization
	"""
	if X.ndim==2:
		nVar,m = X.shape
	else:
		nVar = 1
		m = len(X)
		tmp = np.zeros((2,m),dtype='i')
		tmp[0,:] = X
		X = tmp
	nSamp = m/2
	# initialize matrices
	HD = np.zeros((nSamp,nSamp),dtype= 'i')
	N = np.zeros((nSamp,nSamp),dtype= 'i')
	if countMiss:
		# counting all calls, even missing, in normalization
		N = N+nVar*2
	# if not add only when non missing comparison is made, elsewhere
	
	# loop through all variants
	for k in range(nVar):
		# loop through all samples
		for _i in range(nSamp):
			i = _i*2
			# check missing
			if X[k,i] > 0 and X[k,i+1] > 0:
				gen_i = (X[k,i],X[k,i+1])
				for _j in range(_i):
					j = _j*2
					if X[k,j] > 0 and X[k,j+1] > 0:
						gen_j = (X[k,j],X[k,j+1])
						dist = calcHD(gen_i,gen_j)
						HD[_i,_j] = HD[_i,_j]+dist
						if not countMiss: N[_i,_j] =  N[_i,_j]+2
	N = N+1E-32
	NHD = HD/N
	NHD = NHD + NHD.T
	return(NHD)


### utilties ###

def readBiMat(path,nRow):
	"""Reads a bianary format and spits out a matrix.
	assumes 8 bit unsigned char, will put it 
	an np array by dividing into nRows
	type = unsigned 8 bit integer
	"""
	# open and read the bianary file
	tmp = open(path,'rb').next()
	# one byte for each char 
	nData = len(tmp)
	# get col
	nCol = nData/nRow 

	# need to setup format for conversion
	fmt = str(nData)+'B'
	# save data as np array, uint8
	data = np.array(struct.unpack(fmt,tmp),dtype='uint8')
	# put into matrix 
	X = np.reshape(data,(nRow,nCol))
	
	return(X)



### Scatter calculaitons and tests ###

def empTestDS(D,tagTmp,permList=[100,1E4,1E6]):
	"""perform an empirical test on direct scatter.
	Specifically we find the within and between class distances,
	we calculate the ratio sb/sw and the calculate a one sided permutation test
	that the ratio is larger than random.
	Permutations used to estimate null,
	GPDPerm used when possible to improve estimate.
	D	distance matrix, nxn np array
	tag	class labels as integers 0,..,K, n np int vector 
	permList	list of permutations to be used in order (highest last)
	"""
	minExc = 10
	# pre-process 
	D,tag = pp(D,tagTmp)
	
	# get separation score
	scat_b,scat_w,sep = calcDS(D,tag)

	# calculate p via permutation 
	# apply GPD when possible
	p = np.nan
	k = 0
	for nPerm in permList:
		nPerm = int(nPerm)
		null_sep = permuteSep(D,tag,nPerm)
		# estimate without GDP if possible
		nExc = np.sum(sep<=null_sep)
		if nExc >= minExc: p = float(nExc)/nPerm
		elif not np.isinf(sep): p = gpdPerm.est(sep,null_sep)
		if p == 0: p = np.nan
		if not np.isnan(p): break
	# check to see if still not set
	if np.isnan(p):
		# not enough perms, est best possible
		# note that if p was worse, it would have been set already
		p = 10./permList[-1]


	return(p,sep)

def permuteSep(D,tag,nSamp=1000,v=False):
	"""Calculate the null distribution
	of the separability (ratio of between to within)
	by permuting the class tags
	returns a vector of size nSamp
	for the distribution. 
	D	distance matrix, nxn np array
	tag	class labels as integers 0,..,K, n np int vector
	nSamp 	number of samples in bs, int
	"""
	n = len(tag)
	null_sep = np.zeros(nSamp)
	for i in range(nSamp):
		tmpTag = np.random.permutation(tag)
		scat_b,scat_w,sep = calcDS(D,tmpTag,v)
		null_sep[i] = sep

	return(null_sep)


def calcDS(D,tag,v=False):
	"""calculate the direct scatter (between and within) classes
	and the separability (ratio of between to within).
	D	distance matrix, nxn np array
	tag	class labels as integers 0,..,K, n np int vector 
	"""

	n = len(tag)
	K = np.max(tag)
	nb = 0
	nw = 0
	sb = 0 
	sw = 0
	for i in range(K+1):
		sb = sb + np.sum(D[tag==i][:,tag!=i])
		sw = sw + np.sum(D[tag==i][:,tag==i])
		nk = np.sum(tag==i)
		nb = nb + nk*(n-nk)
		nw = nw + nk**2 - nk

	scat_b = sb/nb
	scat_w = sw/nw
	if scat_b < eps and scat_w < eps:
		sep = 1
		if v: print 'warning no variation in data.'
	elif scat_w < eps:
		sep = np.inf
		if v: print 'warning no variation in within class.'
	else:
		sep = scat_b/scat_w

	return(scat_b,scat_w,sep)



def rTestDS(D,tag):
	"""perform a robust, non-parametric, rank test on direct scatter.
	Specifically we find the within and between class distances,
	we calculate the median and the calculate a one sided Rank Sum test
	that the between class distances are larger then the pair wise distance.
	D	distance matrix, nxn np array
	tag	class labels as integers 0,..,K, n np int vector 
	NOTE: that the distances are pair-wise so they are not 
	independent, ie the samples in SW or SB 
	
	"""
	D,tag = pp(D,tag)

	sw,sb = sepScat(D,tag)

	p,_ = statsUtil.rankSum1S(sb,sw)


	return(p)


def rTestMS(D,tag):
	"""This is a custom method for a specific analysis
	in which we get the mean scatter of D 
	and test if tag==1 > tag==0 in terms of a 
	ranked sum test
	Only use this for the ptb binary label.
	"""
	newTag = np.array(np.zeros(len(tag)),dtype=int)
	newTag[tag=='TRUE'] = 1
	D,tag = pp(D,tag)
	meanScat = calcMeanScat(D)
	# because the preproc (pp) does not 
	# care about the values we have to be careful
	# we are 
	x = meanScat[newTag==1]
	y = meanScat[newTag==0]
	p1s,_ = statsUtil.rankSum1S(x,y)
	p2s,_ = statsUtil.rankSum(x,y)

	return(p1s,p2s)


def rTestMS_cust(D,tag):
	"""This is a custom method for a specific analysis
	in which we get the mean scatter of D 
	and test if tag==1 > tag==0 in terms of a 
	ranked sum test
	Only use this for the ptb binary label.
	"""
	newTag = np.array(np.zeros(len(tag)),dtype=int)
	newTag[tag=='TRUE'] = 1
	D,tag = pp(D,tag)
	meanScat = calcMeanScat(D)
	# because the preproc (pp) does not 
	# care about the values we have to be careful
	# we are 
	x = meanScat[newTag==1]
	y = meanScat[newTag==0]
	p1s,_ = statsUtil.rankSum1S(x,y)
	p2s,_ = statsUtil.rankSum(x,y)

	return(p1s,p2s,meanScat)

def rTestMS_cust1(D,tag):
	"""This is a custom method for a specific analysis
	in which we get the mean scatter of D 
	and test if tag==1 > tag==0 in terms of a 
	ranked sum test
	Only use this for the ptb binary label.
	"""
	newTag = np.array(np.zeros(len(tag)),dtype=int)
	newTag[tag=='TRUE'] = 1
	D,tag = pp(D,tag)
	meanScat = calcMeanScat(D)
	# because the preproc (pp) does not 
	# care about the values we have to be careful
	# we are 
	x = meanScat[newTag==1]
	y = meanScat[newTag==0]
	p1s,_ = statsUtil.rankSum1S(x,y)
	p2s,_ = statsUtil.rankSum(x,y)

	return(p1s,p2s)


def pp(D,tag):
	"""pre process data by checking the distance matrix
	and converting tags to ints
	Tag must be an array of strings
	Also removes missing 'nan' values
	"""
	if not D[0,0]<1E-52 or not np.abs(D[1,0]-D[0,1])<1E-52:
		D = D.T+D
		if not D[0,0]<1E-52 or not np.abs(D[1,0]-D[0,1])<1E-52:
			raise ValueError('Distance Matrix not correct')

	# remove nan values
	# needs to be np array
	tag = np.array(tag,dtype=str) # if already np, this does nothing
	# assuming string so look for 'nan' NOT np.nan
	missInd = tag=='nan'
	tagNew = tag[~missInd]

	
	tagNew,_ = statsUtil.cat2int(tagNew)

	return(D[~missInd][:,~missInd],tagNew)


def sepScat(D,tag):
	"""Separate the scatter measures into 
	between and within classes.  
	Eliminate the redundant pairs and the 
	self comparisons
	returns
	sw - array of all within class scatters
	sb - array of all between class scatters
	"""
	n = len(tag)
	K = np.max(tag)
	sb = np.array([])
	sw = np.array([])
	for i in range(K+1):
		tmp = D[tag==i][:,tag!=i]
		nTmp,mTmp = tmp.shape
		sb = np.append(sb,np.reshape(tmp,nTmp*mTmp))
		tmp = D[tag==i][:,tag==i]
		nTmp,mTmp = tmp.shape
		sw = np.append(sw,np.reshape(tmp,nTmp*mTmp))
	# remove redundancies and self comparisons
	sb.sort()
	sw.sort()
	sw = sw[n:]
	sw = sw[::2]
	sb = sb[::2]
	
	return(sw,sb)


def calcMeanScat(D):
	"""Calculate the mean scatter for each sample
	just average D down the cols, not including self,
	and return
	"""
	sumD = np.sum(D,0)
	n = float(len(D))
	meanScat = sumD/(n-1)
	return(meanScat)


## burden tests and feaure generation
# all code here assumes BINARY features!!!!!
# all code here assumes 2 copies, diploid!!!!!

def getMA(x):
	"""Find the minor allele from the genomic vector x
	assumes binary with 0 as missing, so 1, or 2.
	"""
	ma = 2
	if np.sum(x==2)>np.sum(x==1): ma=1
	return(ma)

def calcHHRef(X):
	"""Calcualte a catigorical variable
	0 1 2, if there is a heterozygous non ref
	allele, 1, or a homozygous non ref allele,2,
	anywhere in the passed region, X 
	(an int geomic vector).
	FOr binary variants this is similar to 
	minor allele (just on occasion 0 and 2 may be switched 
	for a classifier that has no concept of 0 and 2 this is
	a non issue."""
	# merge the columns into one number per sample per varriant 
	# copy merge .
	Xmerg = X[:,::2] + X[:,1::2]
	# assuming binary and 0 = missing we know the following:
	hhrScore = np.max(Xmerg,0)
	hhrScore[hhrScore<3] = 0
	hhrScore[hhrScore==3] = 1
	hhrScore[hhrScore==4] = 2
	return hhrScore

def calcSumMA(X):
	"""returns a vector which contains the 
	sum of minor alleles for each sample over the 
	region defined by X (int genomic matrix).
	"""
	# this does not generally assume binary/ however getMA does 

	# two copies for each sample
	if X.ndim==2:
		n,m=X.shape
		nSamp = m/2
		# init vector 
		sumMA = np.zeros(nSamp,dtype=int)
		# loop through variants
		for i in range(n):
			ma = getMA(X[i])
			sumMA = sumMA + np.array(X[i,::2]==ma,dtype='int8') # check the first copy
			sumMA = sumMA + np.array(X[i,1::2]==ma,dtype='int8') # check the second copy 
	elif X.ndim==1:
		nSamp = len(X)/2
		sumMA = np.zeros(nSamp,dtype=int)
		ma = getMA(X)
		sumMA = sumMA + np.array(X[::2]==ma,dtype='int8') # check the first copy
		sumMA = sumMA + np.array(X[1::2]==ma,dtype='int8') # check the second copy
	
	return(sumMA)
