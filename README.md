# Disneyland Paris Code for Numerical Experiments
-----------------------------------------------------
All numerical experiments were solved in MATLAB using Gurobi optimizer through the Yalmip interface. 

Software required:

  - Matlab (https://www.mathworks.com/products/matlab.html)
  - Gurobi (https://www.gurobi.com/)
  - Yalmip (https://yalmip.github.io/)
  
Numerical Experiments - Section 4
-----------------------------------------------------
The code in file ProjectDisneyland.m solves the optimization problems for the case study in Section 4.1 and the code in file ProjectDisneyland_Family.m is an extension of the last m-file which solves the optimization problems for the case study in Section 4.2. 

In both m-files are defined:
  - The data used for each attraction
  - The decision variables of the optimization problem
  - The constraints of the optimization problem
  - The objective function of the optimization problem
  
The output of ProjectDisneyland.m:
  - Solution of x and u variables
  - Optimal Value
  - Number of visited attractions
  - Total time in minutes spent in the park
  
The output of ProjectDisneyland_Family.m:
  - Solution of x and u variables for all family(in Phase 1), kid with adults(in Phase 2) and for teenagers playing alone(in Phase 2)
  - Optimal Value for all family(in Phase 1), kid with adults(in Phase 2) and for teenagers playing alone(in Phase 2)
  - Number of visited attractions for all family(in Phase 1), kid with adults(in Phase 2) and for teenagers playing alone(in Phase 2)
  - Total time in minutes spent in the park for all family(in Phase 1), kid with adults(in Phase 2) and for teenagers playing alone(in Phase 2)
  - Total Rate for each age group which is included in the family during Phase 1
  - Total Rate for each age group which is included in the family during Phase 1 and Phase 2 together
  
