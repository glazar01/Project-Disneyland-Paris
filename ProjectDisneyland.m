% Author: Georgia Lazaridou
% Last Date Modified: 16/09/2020
% Description: This m-file consists of a code that solves an orienteering optimization problem in Matlab based on Disneyland's Paris data. There are some attractions with different
% attributes(=Data) as defined below and there is also a time limit(=Tmax). This program defines the constraints as well as the objective function of this problem and solves it 
% using Yalmip (with some extra settings) for a specific age group(=Age) which can be 1,2,3 or 4. Each attraction has a different rate for each age group(table Ratings) where
% column 1 = Preschoolers etc. and the selection of column can be made by choosing the right number for RateGroup variable. Finally, the solution of the problem is printed on 
% the console.


% Define Variables
N = 46; %number of attractions
x = binvar (N+2,N+2, 'full'); % attractions binary variables
u = intvar (N+1,1);  % subtours integer variables

% Define Some Values for the constraints
Tmax = 60; % Number of available time in minutes
Age = 3; % Preschoolers = 1, Kid = 2, Tween = 3, Adults/Teens = 4
RateGroup = 3; % Preschoolers = 1, Kid = 2, Tween = 3, Teens = 4, Adults = 5, Kid+Teen+Adult = 6, Preschooler+Tween+Adult = 7

% Define Data Used for each attraction
walk = csvread('MinutesWalkingModified.csv'); % Minutes of walking between attractions
Dur = csvread('DurationModified.csv'); % Duration of each attraction
Wait = csvread('WaitTimeModified.csv'); % Wait time of each attraction
Ratings = csvread('RatingsModified.csv'); % Columns: Preschoolers | Kids | Tweens | Teens | Adults | Kid+Teen+Adult | Preschooler+Tween+Adult
AgeR = csvread('Age RestrictionsModified.csv'); % Columns: Preschoolers, Kids, Tweens, Adults/Teens

% Calculation of total spending time for each attraction
t = zeros(N+2);
for i = 1 : N+2
    for j = 1 : N+2
        t(i,j) = walk(i,j) + Dur(i) + Wait(i);
    end
end

% Choosing the right column of 'Ratings' according to the RateGroup choosen
Rate = Ratings(:,RateGroup);
% Vector consists of restrictions for age group 1,2,3 or 4 (=Column of AgeR)
AgeRes = AgeR(Age,1:N); 


% Define all the constraints of the problem
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
    % if statement to take only the restrictions for the specific age group between the N attractions of table AgeRes
    if (AgeRes(i) == 1) % If the attraction is not permitted (=1)
        for j = 1 : N+1
            Constraints = [Constraints, x(j,i+1) == 0, x(i+1,j+1) == 0];
        end
    end
end

% Define the objective of the problem
xx = x(1:N+1,2:N+1); % xx for multiply with rate in objective because x includes entrance (=first row)
Objective = -sum(sum(xx*Rate));  

% Set some options for YALMIP and solver
options = sdpsettings('verbose',1,'solver','gurobi','gurobi.timelimit',7200);

% Solve the problem
sol = optimize(Constraints,Objective,options);

% Analyze error flags
if sol.problem == 0
 % Extract and display value
 display('solution for tmax=');
 disp(Tmax);
 solution = value(x)
 solution_u = value(u)
 Attractions = value(sum(sum(x))-1) % Calculate the number of visited attractions
 Time = value(sum(sum(x.*t))) % Calculate the total number of minutes wasted
 
else
 display('Hmm, something went wrong!');
 sol.info
 yalmiperror(sol.problem)
 display('solution for tmax=');
 disp(Tmax);
 solution = value(x)
 solution_u = value(u)
 Attractions = value(sum(sum(x))-1) % Calculate the number of visited attractions
 Time = value(sum(sum(x.*t))) % Calculate the total number of minutes wasted
end
end
