#!/usr/bin/env python
#
# 
#     Copyright (C) 2003-2012 Institute for Systems Biology
#                             Seattle, Washington, USA.
# 
#     This library is free software; you can redistribute it and/or
#     modify it under the terms of the GNU Lesser General Public
#     License as published by the Free Software Foundation; either
#     version 2.1 of the License, or (at your option) any later version.
# 
#     This library is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     Lesser General Public License for more details.
# 
#     You should have received a copy of the GNU Lesser General Public
#     License along with this library; if not, write to the Free Software
#     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
# 
# 20130613 RAT
import numpy as np
import scipy.stats as stats
import scipy.linalg
import operator as op
import itertools

def cat2int(y):
	"""change a set of str catigory lables with n
	unique values into an int vector with unique 
	values of numbers from 0 to n-1
	returns the new int vector
		a list of catigories corrisponding to int value
	"""
	unique = list(set(y))
	yNew = np.array(np.zeros(len(y)),dtype=int)
	for i in range(len(unique)):
		yNew[y==unique[i]]=i

	return(yNew,unique)
def rankSum(x,y,forceExact=False):
	"""This is a two-sided Wilcoxon rank sum test,
	robust and nonparametric way to estimate if 
	two samples x and y have diffrent medians.
	Test assumes that two samples x,y are independent but 
	from distributions with the same median.  
	This is tested agains the alternate hyp that 
	the two samples are indipendent but come from 
	distributions with diffrent medians.
	Identical to Mann-Whitney U test.
	By defult the normal approximation is
	made for the test statistic IF the 
	sum of the samples has 15 or more obs;
	ortherwise, we enumerate the combinations
	to find the exact p-value.
	If forceExact==True then the exact p-value
	is calculated regardless of sample size.
	returns 
	p-value - sacalr
	z, approximate test statistic - scalar

	This function is similar to matlab ranksum.
	"""
	# find smallest 
	xLen = len(x)
	yLen = len(y)
	if xLen<=yLen:
		n = xLen
		sampSm = x
		sampLrg = y
	else:
		n = xLen
		sampSm = x
		sampLrg = y
	
	sampComb = np.append(sampSm,sampLrg)
	# check for no varriation
	if np.var(sampComb) == 0: return(np.nan,np.nan)
	# get ranks, possibly tied
	ranks,tieAdj = tiedRank(sampComb)
	ranks = ranks +1
	# get the stat
	xRank = ranks[:n]
	w = np.sum(xRank)
	# calculate the z statistic
	cor = 2*tieAdj/( (xLen+yLen)*(xLen+yLen-1) )
	wVar = xLen*yLen*((xLen+yLen+1) - cor)/12.
	wMean = n*(xLen+yLen+1)/2.
	wC = w-wMean
	z = (wC - .5 * np.sign(wC))/np.sqrt(wVar)

	# more complex methods exist, but I am using this 
	# simple way to deal with small samples
	if forceExact or (xLen+yLen)<15:
		allComb = chooseAllComb(ranks,n)
		allCombSum = np.sum(allComb,1)
		allCombLen = len(allCombSum)
		plo = np.sum(allCombSum<=w)/float(allCombLen)
		phi = np.sum(allCombSum>=w)/float(allCombLen)
		p = np.min([plo,phi]) 
		p = np.min([2*p,1])
	else:
		p = 2 * stats.norm.cdf(-np.abs(z))

	
		
	return(p,z)

def rankSum1S(x,y,forceExact=False):
	"""This is a one-sided Wilcoxon rank sum test,
	robust and nonparametric way to estimate if 
	the median of sample x > y.
	Test assumes that x,y are independent.  
	This is tested aginst the alternate hyp that 
	the two samples are indipendent but come from 
	distributions with the same medians.
	Identical to a one-sided Mann-Whitney U test.
	By defult the normal approximation is
	made for the test statistic IF the 
	sum of the samples has 15 or more obs;
	ortherwise, we enumerate the combinations
	to find the exact p-value.
	If forceExact==True then the exact p-value
	is calculated regardless of sample size.
	returns 
	p-value - sacalr
	z, approximate test statistic - scalar

	This function is similar to matlab ranksum.
	"""
	# find smallest 
	xLen = len(x)
	yLen = len(y)
	n = xLen
	
	sampComb = np.append(x,y)
	# check for no varriation
	if np.var(sampComb) == 0: return(np.nan,np.nan)
	# get ranks, possibly tied
	ranks,tieAdj = tiedRank(sampComb)
	ranks = ranks +1
	# get the stat
	xRank = ranks[:n]
	w = np.sum(xRank)
	# calculate the z statistic
	cor = 2*tieAdj/( (xLen+yLen)*(xLen+yLen-1) )
	wVar = xLen*yLen*((xLen+yLen+1) - cor)/12.
	wMean = n*(xLen+yLen+1)/2.
	wC = w-wMean
	z = (wC - .5 * np.sign(wC))/np.sqrt(wVar)

	# more complex methods exist, but I am using this 
	# simple way to deal with small samples
	if forceExact or (xLen+yLen)<15:
		allComb = chooseAllComb(ranks,n)
		allCombSum = np.sum(allComb,1)
		allCombLen = len(allCombSum)
		p = np.sum(allCombSum>=w)/float(allCombLen)
	else:
		p = stats.norm.cdf(-z)

	
		
	return(p,z)


def tiedRank(x):
	"""Calculates ranks for vector x while considering 
	ties.  If values are the same then the rank assigned
	is the average of the ranks they would span.
	Also returns the tie adjustment used in some tests.
	"""
	n = len(x)
	# sor the values 
	ind = np.argsort(x)
	xSort = x[ind]
	# get the nan values, at the end
	nNan = np.sum(np.isnan(x))
	xLen = n-nNan
	
		
	# find all ties
	tie = xSort[:xLen-1] == xSort[1:xLen]	
	tieInd = np.arange(xLen-1)[tie]
	tieInd = np.append(tieInd,xLen+2)
	nTies = len(tieInd)

	rankTmp = np.arange(float(n)) # the nan does not seem to work with ints 
	rankTmp[xLen:] = np.nan
	tieAdj = 0
	count = 0
	while count<nTies-1:
		tieStart = tieInd[count]
		nTied = 2
		# count number of ties, multiple ties will have consecuitive numbers
		while tieInd[count+1] == tieInd[count]+1:
			count = count + 1
			nTied = nTied + 1
			
		tieAdj = tieAdj + nTied*(nTied-1)*(nTied+1)/2.
		# compute average 
		rankTmp[tieStart:tieStart+nTied] = np.sum(rankTmp[tieStart:tieStart+nTied])/nTied
		count = count+1

	# put in order
	rank = np.zeros(n)
	rank[ind] = rankTmp
	return(rank,tieAdj)

	
	

def nck(n, k):
	"""Calculate n choose k
	N!/K!(N-K)!.  No warning for 
	large numbers that will take 
	forever!"""
	k = min(k, n-k)
	if k == 0: return 1
	numer = reduce(op.mul, xrange(n, n-k, -1))
	denom = reduce(op.mul, xrange(1, k+1))
	return numer//denom

def chooseAllComb(V,k):
	"""choose all possible combinations of 
	k items from the vector V
	blows up fast so be careful
	"""
	return(np.array(list(itertools.combinations(V,k))))
	
