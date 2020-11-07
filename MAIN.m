% This MATLAB R2019b code is for ASCMFDE across minimization problems. 
% For maximization problems, multiply objective function by -1. 

% For suggestions please contact: Zedong Tang (Email: omegatangzd@gmail.com)
clear
pop_M=100; % population size 100
gen=2000; % generation count 1000
selection_pressure = 'elitist'; % choose 'elitist', other selection strategies have not yet been implemented
p_il = 0; % probability of individual learning (BFGA quasi-Newton Algorithm) --> Indiviudal Learning is an IMPORTANT component of the MFDE.
rmp=0.3; % random mating probability
reps = 30; % repetitions 30

for index = 1:28
    if index <= 9
        Tasks = mybenchmark(index);
    elseif index <= 18
        Tasks = benchmark(index-9);
    else
        Tasks = benchmark19(index - 18);
    end
    for t = 1:reps  
        data_result(t)=ASCMFDE(Tasks,pop_M,gen,selection_pressure,rmp,p_il,1,index);
    end
end
