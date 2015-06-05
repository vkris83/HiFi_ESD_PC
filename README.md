# HiFi_ESD_PC
High Fidelity Energy Storage Dispatch in Production Costing Model

This MATLAB code was co-developed by Dr. Trishna Das and Dr. Venkat Krishnan (http://www.ece.iastate.edu/~vkrish/index.html), during our research work at Iowa State University with professor James D. McCalley. Dr. Venkat Krishnan currently works at National Renewable Energy Laboratory (http://www.nrel.gov/analysis/staff/venkat_krishnan.html).

This is a production costing (PC) program integrated with energy storage dispatch model. The master branch contains the codes and data for performing A 2-DAY (48 HOUR) PRODUCTION COSTING SIMULATION. The optimization problem is formulated as a network flow optimization problem, and the data files (nodesinitial.txt and arcsinitial.txt) model the various nodes and arcs of IEEE 24 bus modified RTS system. The user can feed any other system data in a similar format to the PC model. The master branch also contains many folders with various versions of the basic PC code, which can be used to perform many other analysis. They include:

1. A folder for codes and relevant data which can perform YEARLY SIMULATION (i.e., run 182 2-day simulation sequentially) and
2. Another folder for 2-DAY PC SIMULATION WITH 5-MIN ECONOMIC DISPATCH. 
3. One of the folders also contains codes that create a framework which can be used to perform OPTIMAL ALLOCATION OF ENERGY STORAGE in a grid.

The repository also contains a PPT file (PC_Demonstration_Manual.pptx), which will walk through the various aspects of the code, the data required, the I/O aspects and the way one could use it.

The basic MATLAB code for PC has come down through many students to us, upon which the high-fidelity energy storage dispatch model was developed and used for many analysis. Tomlab optimization environment is used to perform the CPLEX based MILP (Unit Committment) and LP (Economic Dispatch) optimization in Matlab. The model development and the consequent analysis studies are all documented in the various publications below. 


References:


1. Das, Trishna, "Performance and Economic Evaluation of Storage Technologies" (2013).Graduate Theses and Dissertations. Paper 13047


2. T. Das, V. Krishnan, and J. D. McCalley, High-Fidelity Dispatch Model of Storage Technologies for Production Costing Studies, IEEE Transactions on Sustainable Energy, vol.5, no.4, pp.1242–1252, Oct. 2014


3. T. Das, V. Krishnan, and J. McCalley, Incorporating cycling costs in generation dispatch program — an economic value stream for energy storage, International Journal of Energy Research, Wiley Online Library, Volume 38, Issue 12, pages 1551–1561, 10 October 2014


4. T. Das, V. Krishnan, and J. D. McCalley, Assessing the benefits and economics of bulk energy storage technologies in the power grid, Applied Energy, Volume 139, pp. 104–118 , 1 February 2015


5. V. Krishnan and T. Das, Optimal allocation of energy storage in a co-optimized electricity market: Benefits assessment and deriving indicators for economic storage ventures, Energy, Available online 8 January 2015


6. D. Nock, V. Krishnan, and J. McCalley, Dispatching Intermittent Wind Resources for Ancillary services via Wind Control and its Impact on Power System Economics, Renewable Energy, Volume 71, November 2014, Pages 396–400


7. M. Howland, V. Krishnan, N. Brown, and J. McCalley, Assessing the Impact of Power Rate Limitation based Wind Control Strategy, Proceedings of the 2014 IEEE PES Transmission & Distribution Conference & Exposition, Chicago USA, April 2014
