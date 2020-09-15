% Author: Georgia Lazaridou
% Last Date Modified: 16/09/2020
% Description: 



% Define Variables
N = 46; %number of attractions
x = binvar (N+2,N+2, 'full'); % attractions binary variables
u = intvar (N+1,1);  % subtours integer variables
%tt = t(1:N+1,2:N+2);

% Define Some Values for the constraints
%Tmax = 60; % number of available minutes
Age = 3; % Preschoolers = 1, Kid = 2, Tween = 3, Adults/Teens = 4
RateGroup = 3; % Preschoolers = 1, Kid = 2, Tween = 3, Teens = 4, Adults = 5, Kid+Teen+Adult = 6, Preschooler+Tween+Adult = 7
%RateGroup = 1; 

% Define Data Used
walk = csvread('MinutesWalkingModified.csv'); % minutes of walking between attractions tij
Dur = csvread('DurationModified.csv');
Wait = csvread('WaitTimeModified.csv');
Ratings = csvread('RatingsModified.csv'); % Preschoolers | Kids | Tweens | Teens | Adults | Kid+Teen+Adult | Preschooler+Tween+Adult
AgeR = csvread('Age RestrictionsModified.csv'); % Preschoolers, Kids, Tweens, Adults/Teens

t = zeros(N+2);
for i = 1 : N+2
    for j = 1 : N+2
        t(i,j) = walk(i,j) + Dur(i) + Wait(i);
    end
end

Rate = Ratings(:,RateGroup);
AgeRes = AgeR(Age,1:N); % vector consists of restristrictions for age group 1,2,3 or 4


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
    % between the N attractions of table AgeRes
    if (AgeRes(i) == 1)
        for j = 1 : N+1
            Constraints = [Constraints, x(j,i+1) == 0, x(i+1,j+1) == 0];
        end
    end
end

% Define an objective
%r = Rate(:,RateGroup);
xx = x(1:N+1,2:N+1); %for multiply with rate in objective
Objective = -sum(sum(xx*Rate));   % 

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
 solution = value(x)
 solution_u = value(u)
 Attractions = value(sum(sum(x))-1)
 Time = value(sum(sum(x.*t)))
 
else
 display('Hmm, something went wrong!');
 sol.info
 yalmiperror(sol.problem)
 display('solution for tmax=');
 disp(Tmax);
 solution = value(x)
 solution_u = value(u)
 Attractions = value(sum(sum(x))-1)
 Time = value(sum(sum(x.*t)))
end
end
