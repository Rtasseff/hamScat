hamScat
=======

Tools and utilities to run Hamming distance and scatter tests on variant data.


Functionality 
---------------
This code is highly specific.  
It is designed to read a special binary version of variant matrices.
Although not required each matrix represents a transcript.
Using this matrix a Hamming distance is calculated.
Then various stats and summaries are calculated, 
with determining separability of groups being the goal.
It is designed to run multiple processes easily by calling "runHamScat\_MP.py"
for each transcript.  Designed with GOLEM in mind.
Results are then parsed.

Dependencies
------------
Requires gpdPerm (https://github.com/Rtasseff/gpdPerm.git),
Numpy, scipy.
GOLEM is expected for the MP but not required.



License
-------

Copyright (C) 2003-2013 Institute for Systems Biology, Seattle, Washington, USA.
 
The Institute for Systems Biology and the authors make no representation about the suitability or accuracy of this software for any purpose, and makes no warranties, either express or implied, including merchantability and fitness for a particular purpose or that the use of this software will not infringe any third party patents, copyrights, trademarks, or other rights. The software is provided "as is". The Institute for Systems Biology and the authors disclaim any liability stemming from the use of this software. This software is provided to enhance knowledge and encourage progress in the scientific community. 
 
This is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
You should have received a copy of the GNU Lesser General Public License along with this library; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
http://www.gnu.org/licenses/lgpl.html


