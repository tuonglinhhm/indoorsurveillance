# This is repository for a signal-based surveillance method

Required 
   - Matlab
   - Communications/Propagation and Channel Models toolbox
   
How to build building structure
   1. Build the structure for the complex/building in RayTracing/CCU.csv 
   
   (Each wall represented by two rows of data 1) [Xstart, Ystart, Zstart]; 2) [Xend, Yend, Zend]
   
   2. Set the routes of the target in RayTracing/main.m in the variable Tx.xyz.  
   
How to run the simulation
   1. Run RayTracing/main.m to collect the ``estimated'' location of the targets based on the Ray Tracing method. The parameters for wireless frequency and path loss constants can be found here too. 
   2. Apply the same setting and run Multi-arrayPositioning/main.m to collect the ``estimated'' location of the targets based on the Single-Anchor Radio Positioning method.
   3. Run IndoorLocalization.m in the root folder to find the other estimated location by difference running times.  
   4. Averaging the estimated locations and compare with the pre-defined locations of the targets.

