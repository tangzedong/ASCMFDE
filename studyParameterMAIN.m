% This MATLAB R2014b code is for EVOLUTIONARY MULTITASKING across minimization problems. 
% For maximization problems, multiply objective function by -1. 

% For suggestions please contact: Wei Zhou (Email: jerryzhou@cqu.edu.cn)

%% Calling the solvers
% For large population sizes, consider using the Parallel Computing Toolbox
% of MATLAB.
% Else, program can be slow.
pop_M=100; % population size 100
pop_S = pop_M;
gen=2000; % generation count 1000
selection_pressure = 'elitist'; % choose either 'elitist' or 'roulette wheel'
p_il = 0; % probability of individual learning (BFGA quasi-Newton Algorithm) --> Indiviudal Learning is an IMPORTANT component of the MFDE.
rmp=0.3; % random mating probability
reps = 20; % repetitions 20
%4,5,6,bad result
nn = 2;

for dimthreshold = [0.7, 0.9, 1]
    nn = nn + 1;
for index = [8, 18, 24]
    if index == 15
        p_il = 0;
    end
    if index <= 9
        Tasks = mybenchmark(index);
    elseif index <= 18
        Tasks = benchmark(index-9);
    else
        Tasks = benchmark19(index - 18);
    end
    data_GFTMFDE=studyParameterMFDE(Tasks,pop_M,gen,selection_pressure,rmp,p_il,reps,dimthreshold,index);   
    save(['GFTMFDE', num2str(index), '-', num2str(nn), '.mat']);
    % "task_for_comparison_with_SOO" compares performance of corresponding task in MFO with SOO.
    % For Instance, In EXAMPLE 1 ...
    % "task_for_comparison_with_SOO" = 1 --> compares 40+-D Rastrin in MFO with 40-D
    % Rastrigin in SOO.
    % "task_for_comparison_with_SOO" = 2 --> compares 30D Ackley in MFO with
    % 30D Ackley in SOO.
%     task_for_comparison_with_SOO = 1;
%     data_SOO_1(index)=SOEA(Tasks(task_for_comparison_with_SOO),pop_S,gen,selection_pressure,p_il,reps);   
% 
%     task_for_comparison_with_SOO = 2;
%     data_SOO_2(index)=SOEA(Tasks(task_for_comparison_with_SOO),pop_S,gen,selection_pressure,p_il,reps);     
end
end