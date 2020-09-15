
Total_Rate_Kid = 0;
Total_Rate_Teen = 0;
Total_Rate_Adult = 0;

for k = 1 : 3
    
% Define Tmax for each of the two groups (all family - separated)     
if (k ==1)
    Tmax = 270; % Tmax for the whole family
else
    Tmax = 270; % Tmax for the part that they are separated
end

% Define Variables
N = 46; %number of attractions
x = binvar (N+2,N+2, 'full'); %attractions, entrance and exit
u = intvar (N+1,1);  %useful avoid subtours variables

% Define Some Values for the constraints
Age1 = 2; % Preschoolers = 1, Kid = 2, Tween = 3, Adults/Teens = 4
Age2 = 4; % Preschoolers = 1, Kid = 2, Tween = 3, Adults/Teens = 4

RateGroup = 6; % Preschoolers = 1, Kid = 2, Tween = 3, Teens = 4, Adults = 5, Kid+Teen+Adult = 6, Preschooler+Tween+Adult = 7, Kid+Adult = 8
RateGroup1 = 2; % for the kid
RateGroup2 = 4; % for the teens
RateGroup3 = 5; % for the adults
RateGroup4 = 8; % for the adults with the kid

% Define Data Used
walk = csvread('MinutesWalkingModified.csv'); % minutes of walking between attractions tij
Dur = csvread('DurationModified.csv');
Wait = csvread('WaitTimeModified.csv');
Ratings = csvread('RatingsModified2.csv'); % Preschoolers | Kids | Tweens | Teens | Adults | Kid+Teen+Adult | Preschooler+Tween+Adult | Kid & Adults
AgeR = csvread('Age RestrictionsModified.csv'); % Preschoolers, Kids, Tweens, Adults/Teens

t = zeros(N+2);
for i = 1 : N+2
    for j = 1 : N+2
        t(i,j) = walk(i,j) + Dur(i) + Wait(i);
    end
end

% Define Rate and AgeRes
Rate = Ratings(:,RateGroup); % the rate used in objective function
Rate1 = Ratings(:,RateGroup1); % ratings for kid
Rate2 = Ratings(:,RateGroup2); % ratings for teens
Rate3 = Ratings(:,RateGroup3); % ratings for adults
Rate4 = Ratings(:,RateGroup4); % ratings for adults and kid
AgeRes1 = AgeR(Age1,1:N); % vector consists of restristrictions for age group 1,2,3 or 4
AgeRes2 = AgeR(Age2,1:N); % vector consists of restristrictions for age group 1,2,3 or 4

% Define constraints T
Constraints = [sum(sum(x.*t)) <= Tmax, sum(x(1,2:N+2)) == 1, sum(x(1:N+1,N+2)) == 1, sum(x(:,1)) == 0, sum(x(N+2,:)) == 0];
for i = 1 : N+1
    Constraints = [Constraints, 2 <= u(i) <= N+2];
end
for i = 2 : N+1
    Constraints = [Constraints, sum(x(i,2:N+2)) <= 1, sum(x(1:N+1,i)) <= 1, sum(x(i,2:N+2))-sum(x(1:N+1,i)) == 0];
end
for i = 1 : N+2
    Constraints = [Constraints, x(i,i) == 0];
end
for i = 2 : N+2
    for j = 2 : N+2
        Constraints = [Constraints, x(i,j)*(N+1) + u(i-1)-u(j-1) <= N];
    end
end
for i = 1 : N
    % if statement to take the restrictions for the specific age group
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
% constraints for teen and kid after the first part of the problem
if(k == 2 || k == 3)
    for i = 1 : N
        if (Restrictions(i) == 1)
            for j = 1 : N+1
                Constraints = [Constraints, x(j,i+1) == 0, x(i+1,j+1) == 0];
            end
        end
    end
end

% Define an objective
%r = Rate(:,RateGroup);
xx = x(1:N+1,2:N+1); %for multiply with rate in objective
if (k == 1)
    Objective = -sum(sum(xx*Rate))
else
    if(k == 2)
        Objective = -sum(sum(xx*Rate4))
    else
        Objective = -sum(sum(xx*Rate2))
    end
end
        
% Set some options for YALMIP and solver
% options = sdpsettings('verbose',1,'solver','quadprog','quadprog.maxiter',100);
options = sdpsettings('verbose',1,'solver','gurobi','gurobi.timelimit',7200);

saveampl(Constraints,Objective,'mymodel46_1');

% Solve the problem
sol = optimize(Constraints,Objective,options);


% Analyze error flags
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

end