% This MATLAB R2014b code is for EVOLUTIONARY MULTITASKING across minimization problems. 
% For maximization problems, multiply objective function by -1. 

% For suggestions please contact: Wei Zhou (Email: jerryzhou@cqu.edu.cn)
clear
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
reps = 10; % repetitions 20
%4,5,6,bad result
n = 1;
for drate = [0.5]%[0.8, 0.3, 0.5, 0.6, 0.8]
for index = 16
%     if index == 15
%         p_il = 0;
%     end
    if index <= 9
        Tasks = mybenchmark(index);
    elseif index <= 18
        Tasks = benchmark(index-9);
    else
        Tasks = benchmark19(index - 18);
    end
    for t = 1:reps  
        data_result(t)=MFDEsuccessrate(Tasks,pop_M,gen,selection_pressure,rmp,p_il,1,index,1);
    end
    save(['MFDE-Rnonalign', num2str(index), '.mat'], 'data_result');
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

% clear data_result
% for index = 16
% %     if index == 15
% %         p_il = 0;
% %     end
%     if index <= 9
%         Tasks = mybenchmark(index);
%     elseif index <= 18
%         Tasks = benchmark(index-9);
%     else
%         Tasks = benchmark19(index - 18);
%     end
%     for t = 1:reps  
%         data_result(t)=ASCMFDEgai2(Tasks,pop_M,gen,selection_pressure,rmp,p_il,1,index, 2);
%     end
%     save(['MFDE-Rnoncentering', num2str(index), '.mat'], 'data_result');
%     % "task_for_comparison_with_SOO" compares performance of corresponding task in MFO with SOO.
%     % For Instance, In EXAMPLE 1 ...
%     % "task_for_comparison_with_SOO" = 1 --> compares 40+-D Rastrin in MFO with 40-D
%     % Rastrigin in SOO.
%     % "task_for_comparison_with_SOO" = 2 --> compares 30D Ackley in MFO with
%     % 30D Ackley in SOO.
% %     task_for_comparison_with_SOO = 1;
% %     data_SOO_1(index)=SOEA(Tasks(task_for_comparison_with_SOO),pop_S,gen,selection_pressure,p_il,reps);   
% % 
% %     task_for_comparison_with_SOO = 2;
% %     data_SOO_2(index)=SOEA(Tasks(task_for_comparison_with_SOO),pop_S,gen,selection_pressure,p_il,reps);     
% end

% clear data_result
% for index = 1:28
% %     if index == 15
% %         p_il = 0;
% %     end
%     if index <= 9
%         Tasks = mybenchmark(index);
%     elseif index <= 18
%         Tasks = benchmark(index-9);
%     else
%         Tasks = benchmark19(index - 18);
%     end
%     parfor t = 1:reps  
%         data_result(t)=ASCMFDEgai2(Tasks,pop_M,gen,selection_pressure,rmp,p_il,1,index, 3);
%     end
%     save(['MFDE-Rcentering', num2str(index), '.mat'], 'data_result');
%     % "task_for_comparison_with_SOO" compares performance of corresponding task in MFO with SOO.
%     % For Instance, In EXAMPLE 1 ...
%     % "task_for_comparison_with_SOO" = 1 --> compares 40+-D Rastrin in MFO with 40-D
%     % Rastrigin in SOO.
%     % "task_for_comparison_with_SOO" = 2 --> compares 30D Ackley in MFO with
%     % 30D Ackley in SOO.
% %     task_for_comparison_with_SOO = 1;
% %     data_SOO_1(index)=SOEA(Tasks(task_for_comparison_with_SOO),pop_S,gen,selection_pressure,p_il,reps);   
% % 
% %     task_for_comparison_with_SOO = 2;
% %     data_SOO_2(index)=SOEA(Tasks(task_for_comparison_with_SOO),pop_S,gen,selection_pressure,p_il,reps);     
% end
n = n + 1;
end