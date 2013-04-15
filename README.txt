** LIBGeFP -- Geometrical FLIRT Phrases for Large Scale Place Recognition in 2D Range Data
** G.D. Tipaldi, L.Spinello, W. Burgard -- Int. Conf. Robotics and Automation (ICRA) 2013
** coded: L. Spinello 2013

Why are Geometrical FLIRT Phrases useful?
=========================================

Place recognition, i.e., the problem of recognizing
if the robot is navigating in an already visited place, is a
fundamental problem in mobile robot navigation. Efﬁcient
solutions to this problem are relevant for effectively localizing
robots and for creating maps in real time. Relatively few methods
have been proposed to efﬁciently solve this problem in very large
environments using 2D range data. In this paper, we introduce
geometrical FLIRT phrases (GFPs) as a novel retrieval method
for very efﬁcient and precise place recognition. GFPs perform
approximate 2D range data matching, have low computational
cost, can handle complicated partial matching patterns and are
robust to noise. Experiments carried out with publicly available
datasets demonstrate that GFPs largely outperform state-of-theart approaches in 2D range-based place recognition in terms of
efﬁciency and recall. We obtain retrieval performances with more
than 85% recall at 99% precision in less than a second, even on
data sets obtained from several kilometer long runs.

see the <a href="http://www.informatik.uni-freiburg.de/~spinello/tipaldiICRA13.pdf">paper</a> 

Installation
=========================================

Download the library. The library relies on cmake to generate the Makefiles.

Go to the root directory of your project and run
$ mkdir build
$ cd build  
$ cmake .. 
$ make
Binaries are generated in build/bin and build/lib

Library can be installed in your system by running
$ make install

The software depends on the following external libraries
Boost >= 1.4 (special_functions/binomial)


How to cite LIBGeFP
=========================================

"Geometrical FLIRT Phrases for Large Scale Place Recognition in 2D Range Data", G. D. Tipaldi, L. Spinello, W. Burgard -- Int. Conf. Robotics and Automation (ICRA) 2013
 

License
=========================================

GFP LIBGeFP Copyright (c) 2013 Luciano Spinello, licensed under GPL ver. 2.0

GFP is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version. 
GFP is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with GFP.  If not, see <http://www.gnu.org/licenses/>.
 */
