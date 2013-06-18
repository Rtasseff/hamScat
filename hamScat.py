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
	""" find genetic distance between 2 varriants
	assume single digit identifiers seperated by 
	a single character sep ie 0/1 or 0|1.
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
def calcNHDMat(X):
	"""calculat the the normalized hamming distance matrix
	This is specific to the genotype matrices.
	Currently we assume that each row is a variant
	each column is 1/2 a sample ie two calls for each sample
	so every 2 col is a sample.
	Assuming integers
	Missing data is a zero 
	Comparisions with any missing data are ignored.
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
						N[_i,_j] =  N[_i,_j]+2
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

def empTestDS(D,tag,minPerm=100,maxPerm=1E8,stepPerm=100):
	"""perform an empirical test on direct scatter.
	Specifically we find the within and between class distances,
	we calculate the ratio sb/sw and the calculate a one sided permutation test
	that the ratio is larger than random.
	Permutations used to estimate null,
	GPDPerm used when possible to improve estimate.
	D	distance matrix, nxn np array
	tag	class labels as integers 0,..,K, n np int vector 
	"""
	# pre-process 
	D,tag = pp(D,tag)
	
	# get separation score
	scat_b,scat_w,sep = calcDS(D,tag)

	# calculate p via permutation 
	# apply GPD when possible
	nPerm = minPerm
	p = np.nan
	while np.isnan(p):
		if nPerm > maxPerm: break
		null_sep = permuteSep(D,tag,nPerm)
		p = gpdPerm.est(sep,null_sep)
		if p == 0: p = np.nan
		nPerm = nPerm * stepPerm

	return(p)

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


def pp(D,tag):
	"""pre process data by checking the distance matrix
	and converting tags to ints
	"""
	if not D[0,0]<1E-52 or not np.abs(D[1,0]-D[0,1])<1E-52:
		D = D.T+D
		if not D[0,0]<1E-52 or not np.abs(D[1,0]-D[0,1])<1E-52:
			raise ValueError('Distance Matrix not correct')
	
	tag,_ = statsUtil.cat2int(tag)

	return(D,tag)


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

