% Author: Georgia Lazaridou
% Last Date Modified: 03/06/2020
% Description: This m-file is an extension of the ProjectDisneyland.m which introduces the idea of solving the optimization problem as in ProjectDisneyland.m but now for different
% age groups in order to satisfy a whole family. So there are two Phases (1 and 2) where in Phase 1 the family plays as a whole and in Phase 2 Kids are with their Parents and the
% Teenagers are playing alone. There is a for loop (k=1,2,3) in order to solve the problem three different times. In Phase 2 (k=2 for Kid and Parents and k=3 for Teenagers) the
% attractions that have been visited in Phase 1 are excluded. Finally the total rate of each family member as well as the solution of the problems are printed in the console after
% the problem is solved.
%
% -------------- HOW TO RUN THIS PROGRAM - STEPS ------------------------
% 
% STEP 0 : Choose your family, you can choose over the following age
% groups: 
%           1) Preschoolers 
%           2) Kid 
%           3) Tween 
%           4) Teen 
%           5) Adults
%
% STEP 1 :  Choose the right Tmax according to the experiment
%
% STEP 2 :  Selection of Age
%
% STEP 3 :  Selection of the RateGroups
% 
% ----------------------------------------------------------------------

% Initializing the total rating for all family members' age groups
Total_Rate_Age1 = 0;
Total_Rate_Age2 = 0;
Total_Rate_Age3 = 0;
Total_Rate_Age4 = 0;

% ATTENTION!!
% ----------------------------------------------------------------------
% First Group:  --> Age1 + Age3
% Second Group: --> Age2 + Age4
% ----------------------------------------------------------------------

% For loop to repeat the optimization problem 
% 1 = All together, 2 = First group, 3 = Second group
for k = 1 : 3
    
    % -------------- STEP 1: Choose the right Tmax ----------------------
    % Define Tmax for each of the two phases (Phase1: All family - Phase2: Separated)     
    if (k ==1)
        Tmax = 270; % Tmax for the family as a whole --> (Phase 1)
    else
        Tmax = 270; % Tmax for the part that they are separated --> (Phase 2)
    end

    % Define Variables
    N = 46;                         % Number of attractions
    x = binvar (N+2,N+2, 'full');   % Attractions, entrance and exit -> the binary variables of the problem
    u = intvar (N+1,1);             % Useful avoid subtours integer variables u(i)

    % --------- STEP 2: Selection of the Age groups ----------------------
    % Preschoolers = 1, Kid = 2, Tween = 3, Adults/Teens = 4
    % Define some values for the constraints
    % Define the two age groups of the family members
    Age1 = 2;                       
    Age2 = 4;                       
    Age3 = 4;
    Age4 = 4;
    
    
    % --------- STEP 3: Selection of the Rate Groups----------------------
    % Preschoolers = 1, Kid = 2, Tween = 3, Teens = 4, Adults = 5,
    % Kid+Teen+Adult = 6, Kid+Adult = 7, Preschooler+Adult = 8, Teen+Kid = 9, Tween+Adult = 10, 
    % Preschooler+Teen+Adult = 11, Preschooler+Kid+Teen+Adult = 12, Preschooler+Tween+Adult = 13, 
    % Kid+Tween+Adult = 14, Kid+Preschooler+Adult = 15
    % --------------------------------------------------------------------
    RateGroup = 6;          % Whole family rating - used in objective phase 1
    RateGroup1 = 5;         % Rate for Age1        
    RateGroup2 = 5;         % Rate for Age2
    RateGroup3 = 5;         % Rate for Age3
    RateGroup4 = 7;         % Rate for Age4
    
    RateGroup5 = 10;        % Rate for Age1 + Age3 (for first group k = 2)
    RateGroup6 = 7;         % Rate for Age2 + Age4 (for second group k = 3)
    
    % Define Data Used from csv files
    walk = csvread('MinutesWalkingModified.csv');   % Minutes of walking between attractions tij
    Dur = csvread('DurationModified.csv');          % Duration of each attraction
    Wait = csvread('WaitTimeModified.csv');         % Wait time of each attraction 
    Ratings = csvread('RatingsModified3.csv');      % Columns: Preschoolers | Kids | Tweens | Teens | Adults | Kid+Teen+Adult | Kid & Adults | Preschooler & Adult | Teen & Kid | Tween & Adult | Preschooler&Teen&Adult | Preschooler&Kid&Teen&Adult | Preschooler+Tween+Adult | Kid&Tween&Adult | Kid+Preschooler+Adult
    AgeR = csvread('Age RestrictionsModified.csv'); % Columns: Preschoolers | Kids | Tweens | Adults/Teens

    % Calculation of table t which containts the total time (Minutes of walking from i to j + Duration of i + Wait time of i)
    t = zeros(N+2);
    for i = 1 : N+2
        for j = 1 : N+2
            t(i,j) = walk(i,j) + Dur(i) + Wait(i);
        end
    end

    % Define Vectors of Rate and Age Restrictions for the age groups of the family
    Rate = Ratings(:,RateGroup);    % The rate used in objective function (the average) for the Phase 1 k=1
    Rate1 = Ratings(:,RateGroup1);  % Ratings for RateGroup1
    Rate2 = Ratings(:,RateGroup2);  % Ratings for RateGroup2
    Rate3 = Ratings(:,RateGroup3);  % Ratings for RateGroup3
    Rate4 = Ratings(:,RateGroup4);  % Ratings for RateGroup4
    Rate5 = Ratings(:,RateGroup5);  % Ratings for RateGroup5 k=2
    Rate6 = Ratings(:,RateGroup6);  % Ratings for RateGroup6 k=3
    AgeRes1 = AgeR(Age1,1:N);       % Vector consists of restristrictions for age group Age1
    AgeRes2 = AgeR(Age2,1:N);       % Vector consists of restristrictions for age group Age2
    AgeRes3 = AgeR(Age3,1:N);       % Vector consists of restristrictions for age group Age3
    AgeRes4 = AgeR(Age4,1:N);       % Vector consists of restristrictions for age group Age4

    % Define constraints T
    Constraints = [sum(sum(x.*t)) <= Tmax, sum(x(1,2:N+2)) == 1, sum(x(1:N+1,N+2)) == 1, sum(x(:,1)) == 0, sum(x(N+2,:)) == 0]; % Constraints (3.6) and (3.3)
    for i = 1 : N+1
        Constraints = [Constraints, 2 <= u(i) <= N+2];                                                                          % Constraint (3.7)
    end
    for i = 2 : N+1
        Constraints = [Constraints, sum(x(i,2:N+2)) <= 1, sum(x(1:N+1,i)) <= 1, sum(x(i,2:N+2))-sum(x(1:N+1,i)) == 0];          % Constraint (3.5)
    end
    for i = 1 : N+2
        Constraints = [Constraints, x(i,i) == 0];                                                                               % Constraint (3.4)
    end
    for i = 2 : N+2
        for j = 2 : N+2
            Constraints = [Constraints, x(i,j)*(N+1) + u(i-1)-u(j-1) <= N];                                                     % Constraint (3.8)
        end
    end
    for i = 1 : N                                                                                                               % Constraints for age restrictions
        % if statement to take the age restrictions for the specific age group
        % between the N attractions of table AgeRes1 and AgeRes2
        if(k == 1)
            if (AgeRes1(i) == 1 || AgeRes2(i) == 1 || AgeRes3(i) == 1 || AgeRes4(i) == 1)
                for j = 1 : N+1
                    Constraints = [Constraints, x(j,i+1) == 0, x(i+1,j+1) == 0];
                end
            end
        end
        
        if(k == 2)
            if (AgeRes1(i) == 1 || AgeRes3(i) == 1)
                for j = 1 : N+1
                    Constraints = [Constraints, x(j,i+1) == 0, x(i+1,j+1) == 0];
                end
            end
        end
    
        if(k == 3)
            if (AgeRes2(i) == 1 || AgeRes4(i) == 1)
                for j = 1 : N+1
                    Constraints = [Constraints, x(j,i+1) == 0, x(i+1,j+1) == 0];
                end
            end
        end
    end
    
    if(k == 2 || k == 3)                                                                                                       % Constraints for the excluded attractions in Phase 2
        for i = 1 : N
            if (Restrictions(i) == 1)
                for j = 1 : N+1
                    Constraints = [Constraints, x(j,i+1) == 0, x(i+1,j+1) == 0];
                end
            end
        end
    end

    % Define the objective function
    xx = x(1:N+1,2:N+1); % For multiply with rate in the objective function
    % The objective is different each time because for each k (=1,2 or 3) the age groups are different and we take a different average of rating(=Rate, Rate4 or Rate 2 accordingly)
    if (k == 1)
        Objective = -sum(sum(xx*Rate))
    else
        if(k == 2)
            Objective = -sum(sum(xx*Rate5))
        else
            Objective = -sum(sum(xx*Rate6))
        end
    end
        
    % Set some options for YALMIP and solve
    options = sdpsettings('verbose',1,'solver','gurobi','gurobi.timelimit',7200);

    % Solve the problem
    sol = optimize(Constraints,Objective,options);

    % Analyze error flags
    % Total_Rate_AgeGroup: The total rate (in both Phases) achieved by this age group (AgeGroup = Kid, Teen or Adult in this case)
    % Total_Rate_Family_AgeGroup: The total rate (in Phase 1) achieved by this age group (AgeGroup = Kid, Teen or Adult in this case)
    if sol.problem == 0
        % Extract and display value
        display('solution for tmax=');
        disp(Tmax);
        display('solution for k=');
        disp(k);
        solution = value(x);
        Visited = sum(solution)
        Attractions = value(sum(sum(x))-1)
        Time = value(sum(sum(x.*t))) 
        Games = value(xx);
        
        if (k == 1)
            R = sum(solution);
            Restrictions = R(2:N+1);
            solution_family = value(x)
            solution_u_family = value(u);
            Total_Rate_Age1 = Total_Rate_Age1 + sum(sum(Games*Rate1));
            Total_Rate_Age2 = Total_Rate_Age2 + sum(sum(Games*Rate2));
            Total_Rate_Age3 = Total_Rate_Age3 + sum(sum(Games*Rate3));
            Total_Rate_Age4 = Total_Rate_Age4 + sum(sum(Games*Rate4));
            Total_Rate_Family_Age1 = sum(sum(Games*Rate1));
            Total_Rate_Family_Age2 = sum(sum(Games*Rate2));
            Total_Rate_Family_Age3 = sum(sum(Games*Rate3));
            Total_Rate_Family_Age4 = sum(sum(Games*Rate4));
        end
 
        if(k == 2)
            solution_u_age1 = value(u);
            solution_u_age3 = value(u);
            solution_age1_age3 = value(x)
            Total_Rate_Age1 = Total_Rate_Age1 + sum(sum(Games*Rate1));
            Total_Rate_Age3 = Total_Rate_Age3 + sum(sum(Games*Rate3));
        end
 
        if(k == 3)
            solution_u_age2 = value(u);
            solution_u_age4 = value(u);
            solution_age2_age4 = value(x)
            Total_Rate_Age2 = Total_Rate_Age2 + sum(sum(Games*Rate2));
            Total_Rate_Age4 = Total_Rate_Age4 + sum(sum(Games*Rate4));
    
            disp('Total_Rate_Family_Age1 = ')
            disp(Total_Rate_Family_Age1)
            disp('Total_Rate_Family_Age2 = ')
            disp(Total_Rate_Family_Age2)
            disp('Total_Rate_Family_Age3 = ')
            disp(Total_Rate_Family_Age3)
            disp('Total_Rate_Family_Age4 = ')
            disp(Total_Rate_Family_Age4)
            disp('solution_u_family = ')
            disp(solution_u_family)
        
            disp('Total_Rate_Age1 = ')
            disp(Total_Rate_Age1)
            disp('solution_u_age1 = ')
            disp(solution_u_age1)
     
            disp('Total_Rate_Age2 = ')
            disp(Total_Rate_Age2)
            disp('solution_u_age2 = ')
            disp(solution_u_age2)
     
            disp('Total_Rate_Age3 = ')
            disp(Total_Rate_Age3)
            disp('solution_u_age3 = ')
            disp(solution_u_age3)
            
            disp('Total_Rate_Age4 = ')
            disp(Total_Rate_Age4)
            disp('solution_u_age4 = ')
            disp(solution_u_age4)
        end

 
    else
        display('Hmm, something went wrong!');
        sol.info
        yalmiperror(sol.problem)
        display('solution for tmax=');
        disp(Tmax);
        display('solution for k=');
        disp(k);
        solution = value(x);
        Visited = sum(solution)
        Attractions = value(sum(sum(x))-1)
        Time = value(sum(sum(x.*t)))
        Games = value(xx);
        
        if (k == 1)
            R = sum(solution);
            Restrictions = R(2:N+1);
            solution_u_family = value(u);
            solution_family = value(x)
            Total_Rate_Age1 = Total_Rate_Age1 + sum(sum(Games*Rate1));
            Total_Rate_Age2 = Total_Rate_Age2 + sum(sum(Games*Rate2));
            Total_Rate_Age3 = Total_Rate_Age3 + sum(sum(Games*Rate3));
            Total_Rate_Age4 = Total_Rate_Age4 + sum(sum(Games*Rate4));
            Total_Rate_Family_Age1 = sum(sum(Games*Rate1));
            Total_Rate_Family_Age2 = sum(sum(Games*Rate2));
            Total_Rate_Family_Age3 = sum(sum(Games*Rate3));
            Total_Rate_Family_Age4 = sum(sum(Games*Rate4));
        end
 
        if(k == 2)
            solution_u_age1 = value(u);
            solution_u_age3 = value(u);
            solution_age1_age3 = value(x)
            Total_Rate_Age1 = Total_Rate_Age1 + sum(sum(Games*Rate1));
            Total_Rate_Age3 = Total_Rate_Age3 + sum(sum(Games*Rate3)); 
        end
 
        if(k == 3)
            solution_u_age2 = value(u);
            solution_age2 = value(x)
            solution_u_age4 = value(u);
            solution_age4 = value(x)
            Total_Rate_Age2 = Total_Rate_Age2 + sum(sum(Games*Rate2));
            Total_Rate_Age4 = Total_Rate_Age4 + sum(sum(Games*Rate4));
 
            disp('Total_Rate_Family_Age1 = ')
            disp(Total_Rate_Family_Age1)
            disp('Total_Rate_Family_Age2 = ')
            disp(Total_Rate_Family_Age2)
            disp('Total_Rate_Family_Age3 = ')
            disp(Total_Rate_Family_Age3)
            disp('Total_Rate_Family_Age4 = ')
            disp(Total_Rate_Family_Age4)
            disp('solution_u_family = ')
            disp(solution_u_family)
     
            disp('Total_Rate_Age1 = ')
            disp(Total_Rate_Age1)
            disp('solution_u_age1 = ')
            disp(solution_u_age1)
     
            disp('Total_Rate_Age2 = ')
            disp(Total_Rate_Age2)
            disp('solution_u_age2 = ')
            disp(solution_u_age2)
     
            disp('Total_Rate_Age3 = ')
            disp(Total_Rate_Age3)
            disp('solution_u_age3 = ')
            disp(solution_u_age3)
            
            disp('Total_Rate_Age4 = ')
            disp(Total_Rate_Age4)
            disp('solution_u_age4 = ')
            disp(solution_u_age4)
        end
    end
    
end % End of for loop k=1:3

% End of the program
