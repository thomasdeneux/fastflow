--------------------------------------------------------------------------------
CONTENT

The FASTFLOW program is a Matlab graphical user interface (GUI) toolbox to 
estimate blood velocity in optical imaging recording. 

Matlab code and a demo dataset are available at 
https://trac.int.univ-amu.fr/fastflow/.

--------------------------------------------------------------------------------
LICENSE

Copyright (C) 2005-2011 Thomas Deneux.
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

The full GNU General Public License can be found in the license.pdf file.
See also <http://www.gnu.org/licenses/>.

--------------------------------------------------------------------------------
CONTACT

Please contact Thomas Deneux (thomas.deneux@incm.cnrs-mrs.fr) for help and bug
reports.

--------------------------------------------------------------------------------
STARTUP

No installation is needed. 

Once the main program folder has been copy at any desired location, run the file
'StartFastFlow.m'.

Note that this file will automatically add the appropriate folders to the Matlab
path prior to opening the graphic interface.

Note also that it will attempt to compile two MEX-file funtions, 
'fast_regenergy' and 'fast_regenergy_single'. If the compilation should fail
(for exemple, if no compiler is installed), it will attempt compiling each time 
the program is started again, and it will not be possible to perform 
coregistration (other functionalities will be unaffected however).

--------------------------------------------------------------------------------
DOCUMENTATION

No documentation exists yet.

Download the demonstration data from https://trac.incm.cnrs-mrs.fr/fastflow/ and
follow instructions from the README_data file to get used to the software. 



