% This MATLAB R2014b code is for EVOLUTIONARY MULTITASKING across minimization problems. 
% For maximization problems, multiply objective function by -1. 

% For suggestions please contact: Wei Zhou (Email: jerryzhou@cqu.edu.cn)
clear

%% Calling the solvers
% For large population sizes, consider using the Parallel Computing Toolbox
% of MATLAB.
% Else, program can be slow.
pop_M=50; % population size 100
gen=2000; % generation count 1000
reps = 30; % repetitions 20
%4,5,6,bad result
n = 1;
for index = 1
    switch index
        case 1
            Tasks = pole2benchmark(1);
            Task = Tasks(1);
        case 2
            Tasks = pole2benchmark(1);
            Task = Tasks(2);
        case 3
            Tasks = pole2benchmark(2);
            Task = Tasks(2);
    end
    
    parfor r =1:reps
        data_result(r)=SOEA(Task,pop_M,gen,'elitist',0,1);   
    end
    save(['SDEpole', num2str(index), '.mat']);   
end