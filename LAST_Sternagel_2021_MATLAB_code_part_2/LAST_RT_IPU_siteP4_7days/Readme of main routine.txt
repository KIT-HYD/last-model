### 
General:

This code comprises routines considering solute transport under
preferential flow conditions and mixing between macropores and soil
matrix and therewith extending the original particle based
Lagrangian Model of Zehe and Jackisch (2016). This extension was done
within the scope of the master thesis of Alexander Sternagel (2018).
Further information on the theoretical background and the development can
be obtained from these two references. The functionality and handling of
this code/model are explained in the different readme-files within the
folder structure.

Zehe, E. and Jackisch, C.: A Lagrangian model for soil water dynamics during rainfall-driven conditions,
Hydrol. Earth Syst. Sci., 20, 3511–3526, doi:10.5194/hess-20-3511-2016, 2016

Sternagel, A.: Modelling of Solute Transport and Travel Time
Distributions Using a Particle-Based Lagrangian Model, 2018

###
Folder Structure:

STLM: "Solute Transport Lagrangian Model", the main model script/routine
init_STLM: initialisation routine of the main routine in which general and not adjustable parameters are calculated
raw data: contains original data set
source: contains subfolders with all input data used in the main script STLM
	every subfolder contains an own readme file describing the respective input data
functions: contains subfolders with all functions used in the main script STLM
	   every subfolder contains an own readme file describing the respective function

###
Predefining parameters:

global and soil matrix parameters can be predefined at the beginning of the main script STLM
parameters of the preferential flow domain (pfd) can be predefined in the subscript "prefflow_domain" in the folder 'source\initial'

###
List of variables and parameters:
