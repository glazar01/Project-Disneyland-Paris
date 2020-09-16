% Author: Georgia Lazaridou
% Last Date Modified: 16/09/2020
% Description: This m-file is an extension of the ProjectDisneyland.m which introduces the idea of solving the optimization problem as in ProjectDisneyland.m but now for different
% age groups in order to satisfy a whole family. So there are two Phases (1 and 2) where in Phase 1 the family plays as a whole and in Phase 2 Kids are with their Parents and the
% Teenagers are playing alone. There is a for loop (k=1,2,3) in order to solve the problem three different times. In Phase 2 (k=2 for Kid and Parents and k=3 for Teenagers) the
% attractions that have been visited in Phase 1 are excluded. Finally the total rate of each family member as well as the solution of the problems are printed in the console after
% the problem is solved.

% Initializing the total rating for all family members' age groups
Total_Rate_Kid = 0;
Total_Rate_Teen = 0;
Total_Rate_Adult = 0;

% For loop to repeat the optimization problem (1 = All together, 2 = Kid with parents, 3 = Teenagers)
for k = 1 : 3
    
    % Define Tmax for each of the two phases (Phase1: All family - Phase2: Separated, kid with adults and teenagers alone)     
    if (k ==1)
        Tmax = 270; % Tmax for the family as a whole (Phase 1)
    else
        Tmax = 270; % Tmax for the part that they are separated (Phase 2)
    end

    % Define Variables
    N = 46;                         % Number of attractions
    x = binvar (N+2,N+2, 'full');   % Attractions, entrance and exit -> the binary variables of the problem
    u = intvar (N+1,1);             % Useful avoid subtours integer variables u(i)

    % Define some values for the constraints
    % Define the two age groups of the family members, in our case study we have only Kid (=2)and Adults/Teens (=4)
    Age1 = 2;                       % Preschoolers = 1, Kid = 2, Tween = 3, Adults/Teens = 4
    Age2 = 4;                       % Preschoolers = 1, Kid = 2, Tween = 3, Adults/Teens = 4
    
    % Selection of the Rate Groups
    % Preschoolers = 1, Kid = 2, Tween = 3, Teens = 4, Adults = 5, Kid+Teen+Adult = 6, Preschooler+Tween+Adult = 7, Kid+Adult = 8
    RateGroup = 6;                  % Select 6 for the average rate of Kid, Teens and Adults
    RateGroup1 = 2;                 % Select 2 for the Kid
    RateGroup2 = 4;                 % Select 4 for the Teens
    RateGroup3 = 5;                 % Select 5 for the Adults
    RateGroup4 = 8;                 % Select 8 for the Adults with the Kid

    % Define Data Used from csv files
    walk = csvread('MinutesWalkingModified.csv');   % Minutes of walking between attractions tij
    Dur = csvread('DurationModified.csv');          % Duration of each attraction
    Wait = csvread('WaitTimeModified.csv');         % Wait time of each attraction 
    Ratings = csvread('RatingsModified2.csv');      % Columns: Preschoolers | Kids | Tweens | Teens | Adults | Kid+Teen+Adult | Preschooler+Tween+Adult | Kid & Adults
    AgeR = csvread('Age RestrictionsModified.csv'); % Columns: Preschoolers | Kids | Tweens | Adults/Teens

    % Calculation of table t which containts the total time (Minutes of walking from i to j + Duration of i + Wait time of i)
    t = zeros(N+2);
    for i = 1 : N+2
        for j = 1 : N+2
            t(i,j) = walk(i,j) + Dur(i) + Wait(i);
        end
    end

    % Define Vectors of Rate and Age Restrictions for the age groups of the family
    Rate = Ratings(:,RateGroup);    % The rate used in objective function (the average) for the Phase 1
    Rate1 = Ratings(:,RateGroup1);  % Ratings for kid
    Rate2 = Ratings(:,RateGroup2);  % Ratings for teens
    Rate3 = Ratings(:,RateGroup3);  % Ratings for adults
    Rate4 = Ratings(:,RateGroup4);  % Ratings for adults and kid
    AgeRes1 = AgeR(Age1,1:N);       % Vector consists of restristrictions for age group Age1
    AgeRes2 = AgeR(Age2,1:N);       % Vector consists of restristrictions for age group Age2

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
        if(k == 1 || k == 2)
            if (AgeRes1(i) == 1 || AgeRes2(i) == 1)
                for j = 1 : N+1
                    Constraints = [Constraints, x(j,i+1) == 0, x(i+1,j+1) == 0];
                end
            end
        end
    
        if(k == 3)
            if (AgeRes2(i) == 1)
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
            Objective = -sum(sum(xx*Rate4))
        else
            Objective = -sum(sum(xx*Rate2))
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
            Total_Rate_Kid = Total_Rate_Kid + sum(sum(Games*Rate1));
            Total_Rate_Teen = Total_Rate_Teen + sum(sum(Games*Rate2));
            Total_Rate_Adult = Total_Rate_Adult + sum(sum(Games*Rate3));
            Total_Rate_Family_Kid = sum(sum(Games*Rate1));
            Total_Rate_Family_Teen = sum(sum(Games*Rate2));
            Total_Rate_Family_Adult = sum(sum(Games*Rate3));
        end
 
        if(k == 2)
            solution_u_kid = value(u);
            solution_u_adult = value(u);
            solution_kid_adult = value(x)
            Total_Rate_Kid = Total_Rate_Kid + sum(sum(Games*Rate1));
            Total_Rate_Adult = Total_Rate_Adult + sum(sum(Games*Rate3));
        end
 
        if(k == 3)
            solution_u_teen = value(u);
            solution_teen = value(x)
            Total_Rate_Teen = Total_Rate_Teen + sum(sum(Games*Rate2));
    
            disp('Total_Rate_Family_Kid = ')
            disp(Total_Rate_Family_Kid)
            disp('Total_Rate_Family_Teen = ')
            disp(Total_Rate_Family_Teen)
            disp('Total_Rate_Family_Adult = ')
            disp(Total_Rate_Family_Adult)
            disp('solution_u_family = ')
            disp(solution_u_family)
        
            disp('Total_Rate_Kid = ')
            disp(Total_Rate_Kid)
            disp('solution_u_kid = ')
            disp(solution_u_kid)
     
            disp('Total_Rate_Teen = ')
            disp(Total_Rate_Teen)
            disp('solution_u_teen = ')
            disp(solution_u_teen)
     
            disp('Total_Rate_Adult = ')
            disp(Total_Rate_Adult)
            disp('solution_u_adult = ')
            disp(solution_u_adult)
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
            Total_Rate_Kid = Total_Rate_Kid + sum(sum(Games*Rate1));
            Total_Rate_Teen = Total_Rate_Teen + sum(sum(Games*Rate2));
            Total_Rate_Adult = Total_Rate_Adult + sum(sum(Games*Rate3));
            Total_Rate_Family_Kid = sum(sum(Games*Rate1));
            Total_Rate_Family_Teen = sum(sum(Games*Rate2));
            Total_Rate_Family_Adult = sum(sum(Games*Rate3));
        end
 
        if(k == 2)
            solution_u_kid = value(u);
            solution_u_adult = value(u);
            solution_kid_adult = value(x)
            Total_Rate_Kid = Total_Rate_Kid + sum(sum(Games*Rate1));
            Total_Rate_Adult = Total_Rate_Adult + sum(sum(Games*Rate3)); 
        end
 
        if(k == 3)
            solution_u_teen = value(u);
            solution_teen = value(x)
            Total_Rate_Teen = Total_Rate_Teen + sum(sum(Games*Rate2));
 
            disp('Total_Rate_Family_Kid = ')
            disp(Total_Rate_Family_Kid)
            disp('Total_Rate_Family_Teen = ')
            disp(Total_Rate_Family_Teen)
            disp('Total_Rate_Family_Adult = ')
            disp(Total_Rate_Family_Adult)
            disp('solution_u_family = ')
            disp(solution_u_family)
     
            disp('Total_Rate_Kid = ')
            disp(Total_Rate_Kid)
            disp('solution_u_kid = ')
            disp(solution_u_kid)
     
            disp('Total_Rate_Teen = ')
            disp(Total_Rate_Teen)
            disp('solution_u_teen = ')
            disp(solution_u_teen)
     
            disp('Total_Rate_Adult = ')
            disp(Total_Rate_Adult)
            disp('solution_u_adult = ')
            disp(solution_u_adult)
        end
    end
    
end % End of for loop k=1:3

% End of the program
