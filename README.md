hamScat
=======

Tools and utilities to run Hamming distance and scatter tests on variant data.


Functionality 
---------------
Python 2.7 scripts.
Parts of this code were written for a very specific purpose, while others are more general.
It is designed to read a special binary version of variant matrices.
Although not required each matrix represents a transcript.
Using this matrix a Hamming distance is calculated between selected samples.
Then various stats and summaries are calculated, 
with determining separability of groups being the primary calculation.
Separability in this case signifies how distant samples of different classes are
compared to the distance of samples in the same class.
The primary output is a p-value which tests this ratio against 
the null distribution (random) as determined by permutation.
It is designed to run multiple processes easily by calling "runHamScat\_MP.py"
for each transcript, with the distribution code GOLEM in mind.
Separate scripts are available to parse and check the results.

Contents Summary
------------------
note the terms scat or hamScat signify the main purpose, calculation of hamming scatter separability, while the term FM signifies extra code written to generate feature matrices without using sample meta data.

CreatJobList\_\*.py	scripts to create the appropriate golem job list    
getSampleMeta.py	project specific file to gather the appropriate meta data   
hamScat.py		main code for running the tests, includes functions to run the tests, calculate hamming distance and open appropriate binary variant region files.    
parse\_\*\_MP.py	parses the resulting output from golem into an easy to read tsv file which is ordered and checked for errors    
run\_\*\_MP.py		code called by golem to run the appropriate functions    
sampleMeta\_selectFeature.dat	listing of features in the sample meta data (specific to project) and can be used to select features for testing   
statsUtil.py		python functions for some statistical utilities   


Dependencies
------------
Python 2.7
Requires gpdPerm (https://github.com/Rtasseff/gpdPerm.git),
Numpy, scipy.
GOLEM is expected for the MP but not required.

Code to generate the binary variant matrices was written by RK and is available only on ISB servers.
Email RTasseff for more info.



License
-------

Copyright (C) 2003-2013 Institute for Systems Biology, Seattle, Washington, USA.
 
The Institute for Systems Biology and the authors make no representation about the suitability or accuracy of this software for any purpose, and makes no warranties, either express or implied, including merchantability and fitness for a particular purpose or that the use of this software will not infringe any third party patents, copyrights, trademarks, or other rights. The software is provided "as is". The Institute for Systems Biology and the authors disclaim any liability stemming from the use of this software. This software is provided to enhance knowledge and encourage progress in the scientific community. 
 
This is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
You should have received a copy of the GNU Lesser General Public License along with this library; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
http://www.gnu.org/licenses/lgpl.html


draft of run protocol
------------

Get current VCF
currently we require a single concatenated file for the VCF

Get sample info 
For the code, we will need information on the samples and we will be relying on the order of the samples in the genomic data to match the clinical data.
Get an ordered list of sample ID’s. 
Produce a tab delimitated file of the sample ids from the current VCF file. Currently this is a bit messy I use: $ zcat df4\_merged.vcf.gz | grep -m 1 CHR | cut -f 10- > sampleID.dat
(note the above assumes 9 prefix fields, so we need 10 and on)
consider leaving this output in the DF dir
Get Sample meta data
This is written into a script in /titan/cancerregulome9/ITMI\_PTB/users/rtasseff/DF4/ called getSampMeta.py
open script to fill in the names and dirs and then run 
consider putting the output in the DF dir
relies on a file received from clinical side, may have to edit if file is different
puts the output into a standard format for later use

Parse out transcript mats
We need to parse the existing VCF into binary matrices, one for each transcript.  This process is done in two parts, that was created by RK.  The latest version will always be in /proj/ilyalab/rkramer/ and I have local copies, new as of this time, on /Users/rtasseff/Projects/INOVA/DF4/, which also contains a readme called readme\_split.txt.
First run the split code. find the number of columns using: $ zcat df4\_merged.vcf.gz | head -n 100 | awk -F $'\t' '{print NF}' | uniq
this is used as an input 
run the script split-transcripts, for example: nohup zcat df4\_merged.vcf.gz | ./split-transcripts -columns <column count> -root=<path>  1> out 2> err &
help can be found using: $ ./split-transcripts -help
consider using a sub dir under DF dir for the matrices
Run the indexing parser,
the script is called mkindex.sh 
it will produce a file called index.tab which gives the info on the transcripts, including where they are in the folder structure. 
example call: ./mkindex.sh out <sample count>
sample count can be found by subtracting 9 from the column count, at least as of this time

Running hamScat
You need to get the sample data.  Currently we have the script getSampMeta.dat
you can find this setup on the server /titan/cancerregulome9/ITMI\_PTB/users/rtasseff/DF4/
it is very specific to this to this setup
produces the sampleMeta.dat and the sampleMeta\_selectFeature.dat
Move all needed files to the run dir, including the hamScat files. we suggest to also move the sampleMeta\_selectFeature.dat file, as it may be changed for different runs
Chose the target features you would like to test against by including the term ‘True’ next to corresponding feature names in sampleMeta\_selectFeature.dat.
create the golem job list using creatJobList\_scat.py
run the job list on golem, this will call the program run\_hamScat\_MP.py, note that other test can be run but are not setup to run using the simple selection of features so they have to be manually entered into run\_hamScat\_MP.py (some commented out examples are present.
note the job id, you will need it to parse the file output
parse results by running parse\_hamScat\_MP.py,
open the file and edit to include input and output names and dirs
will check for completeness and consistency as well as order the results by whats in the original index.tab file


