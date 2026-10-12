# Permafrost-Methane-MultiTracer-and-Monte-Carlo-Simulation

These scripts are for manuscript 'Simulated Detection of Methane Emissions from Ancient Arctic Permafrost with Atmospheric Tracer Measurements'

The scripts to produce Figure 5 and Figure S1 three parameter simulation are in MainSimulation_*.m. 
Sensitivity tests shown in Figure 6 are in SensitivityTest.m, which vary PF vs three parameters as in the main simulation and produce stacked bar uncertainty contribution by removing parameter unc one by one. 
JointConditionEvaluation.m reproduce the joint condition table in Table 3. Figure 1 produce the joint condition for different step of fixed FF. Figure 2 provide joint condition for fixed 100ppb of non-permafrost CH4 and varies PF amount, FF and D14C. Figure 3 shows the joint condition for different step of PF CH4 and their corresponding FF and D14C. The third figure was used to produce the data in Table 3. 

