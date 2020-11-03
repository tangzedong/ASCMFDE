function [Tasks, g1, g2] = mybenchmark(index)
%BENCHMARK function
%   Input
%   - index: the index number of problem set
%
%   Output:
%   - Tasks: benchmark problem set
%   - g1: global optima of Task 1
%   - g2: global optima of Task 2
    switch(index)
        case 10 % complete intersection with high similarity, Griewank and Rastrigin
            load('CI_H.mat')  % loading data from folder .\Tasks
            dim = 50;
            Tasks(1).dims = dim;    % dimensionality of Task 1
            Tasks(1).fnc = @(x)Griewank(x,Rotation_Task1,GO_Task1);
            Tasks(1).Lb=-100*ones(1,dim);   % Upper bound of Task 1
            Tasks(1).Ub=100*ones(1,dim);    % Lower bound of Task 1
            
            Tasks(2).dims = dim;    % dimensionality of Task 2
            Tasks(2).fnc = @(x)Rastrigin(x,Rotation_Task2,GO_Task2);
            Tasks(2).Lb=-50*ones(1,dim);    % Upper bound of Task 2
            Tasks(2).Ub=50*ones(1,dim);     % Lower bound of Task 2
            
            g1 = GO_Task1;
            g2 = GO_Task2;
        case 11 % complete intersection with medium similarity, Ackley and Rastrigin
            load('CI_M.mat')
            dim = 50;
            Tasks(1).dims = dim;
            Tasks(1).fnc = @(x)Ackley(x,Rotation_Task1,GO_Task1);
            Tasks(1).Lb=-50*ones(1,dim);
            Tasks(1).Ub=50*ones(1,dim);
            
            Tasks(2).dims = dim;
            Tasks(2).fnc = @(x)Rastrigin(x,Rotation_Task2,GO_Task2);
            Tasks(2).Lb=-50*ones(1,dim);
            Tasks(2).Ub=50*ones(1,dim);
            
            g1 = GO_Task1;
            g2 = GO_Task2;
        case 12 % complete intersection with low similarity, Ackley and Schwefel
            load('CI_L.mat')
            dim = 50;
            GO_Task1 = 0*ones(1,dim);
            Tasks(1).dims = dim;
            Tasks(1).fnc = @(x)Ackley(x,Rotation_Task1,GO_Task1);
            Tasks(1).Lb=-50*ones(1,dim);
            Tasks(1).Ub=50*ones(1,dim);
            
            
            Tasks(2).dims = dim;
            Tasks(2).fnc = @(x)Schwefel(x);
            Tasks(2).Lb=-500*ones(1,dim);
            Tasks(2).Ub=500*ones(1,dim);
            
            g1 = GO_Task1;
            g2 = 420.9687*ones(1,dim);
        case 13 % partially intersection with high similarity, Rastrigin and Sphere
            load('PI_H.mat')
            dim = 50;
            Tasks(1).dims = dim;
            Tasks(1).fnc = @(x)Rastrigin(x,Rotation_Task1,GO_Task1);
            Tasks(1).Lb=-50*ones(1,dim);
            Tasks(1).Ub=50*ones(1,dim);
            
            Tasks(2).dims = dim;
            Tasks(2).fnc = @(x)Sphere(x,GO_Task2);
            Tasks(2).Lb=-100*ones(1,dim);
            Tasks(2).Ub=100*ones(1,dim);
            
            g1 = GO_Task1;
            g2 = GO_Task2;
        case 14 % partially intersection with medium similarity, Ackley and Rosenbrock
            load('PI_M.mat')
            dim = 50;
            Tasks(1).dims = dim;
            Tasks(1).fnc = @(x)Ackley(x,Rotation_Task1,GO_Task1);
            Tasks(1).Lb=-50*ones(1,dim);
            Tasks(1).Ub=50*ones(1,dim);
            
            Tasks(2).dims = dim;
            Tasks(2).fnc = @(x)Rosenbrock(x);
            Tasks(2).Lb=-50*ones(1,dim);
            Tasks(2).Ub=50*ones(1,dim);
            
            g1 = GO_Task1;
            g2 = ones(1,dim);
        case 15 % partially intersection with low similarity, Ackley and Weierstrass
            load('PI_L.mat')
            dim = 50;
            Tasks(1).dims = dim;
            Tasks(1).fnc = @(x)Ackley(x,Rotation_Task1,GO_Task1);
            Tasks(1).Lb=-50*ones(1,dim);
            Tasks(1).Ub=50*ones(1,dim);
            
            dim = 25;
            Tasks(2).dims = dim;
            Tasks(2).fnc = @(x)Weierstrass(x,Rotation_Task2,GO_Task2);
            Tasks(2).Lb=-0.5*ones(1,dim);
            Tasks(2).Ub=0.5*ones(1,dim);
            
            g1 = GO_Task1;
            g2 = GO_Task2;
        case 16 % no intersection with high similarity, Rosenbrock and Rastrigin
            load('NI_H.mat')
            dim = 50;
            Tasks(1).dims = dim;
            Tasks(1).fnc = @(x)Rosenbrock(x);
            Tasks(1).Lb=-50*ones(1,dim);
            Tasks(1).Ub=50*ones(1,dim);
            
            Tasks(2).dims = dim;
            Tasks(2).fnc = @(x)Rastrigin(x,Rotation_Task2,GO_Task2);
            Tasks(2).Lb=-50*ones(1,dim);
            Tasks(2).Ub=50*ones(1,dim);
            
            g1 = ones(1,dim);
            g2 = GO_Task2;
        case 17 % no intersection with medium similarity, Griewank and Weierstrass
            load('NI_M.mat')
            dim = 50;
            Tasks(1).dims = dim;
            Tasks(1).fnc = @(x)Griewank(x,Rotation_Task1,GO_Task1);
            Tasks(1).Lb=-100*ones(1,dim);
            Tasks(1).Ub=100*ones(1,dim);
            
            Tasks(2).dims = dim;
            Tasks(2).fnc = @(x)Weierstrass(x,Rotation_Task2,GO_Task2);
            Tasks(2).Lb=-0.5*ones(1,dim);
            Tasks(2).Ub=0.5*ones(1,dim);
            
            g1 = GO_Task1;
            g2 = GO_Task2;
            
        case 18 % no overlap with low similarity, Rastrigin and Schwefel
            load('NI_L.mat')
            dim = 50;
            Tasks(1).dims = dim;
            Tasks(1).fnc = @(x)Rastrigin(x,Rotation_Task1,GO_Task1);
            Tasks(1).Lb=-50*ones(1,dim);
            Tasks(1).Ub=50*ones(1,dim);
            
            Tasks(2).dims = dim;
            Tasks(2).fnc = @(x)Schwefel(x);
            Tasks(2).Lb=-500*ones(1,dim);
            Tasks(2).Ub=500*ones(1,dim);
            
            g1 = GO_Task1;
            g2 = 420.9687*ones(1,dim);
        case 1 % complete intersection with high similarity, Griewank and Rastrigin
            load('CI_H.mat')  % loading data from folder .\Tasks
            dim = 30;
            Rotation_Task1 = eye(dim);
            GO_Task1 = zeros(1,dim);
            GO_Task2 = zeros(1,dim);
            Tasks(1).dims = dim;    % dimensionality of Task 1
            Tasks(1).fnc = @(x)Sphere1(x,Rotation_Task1,GO_Task1);
            Tasks(1).Lb=-50*ones(1,dim);   % Upper bound of Task 1
            Tasks(1).Ub=50*ones(1,dim);    % Lower bound of Task 1
            
            Tasks(2).dims = dim;    % dimensionality of Task 2
            Rotation_Task2 = eye(dim);
            Tasks(2).fnc = @(x)Rosenbrock1(x,Rotation_Task2,GO_Task2);
            Tasks(2).Lb=-50*ones(1,dim);    % Upper bound of Task 2
            Tasks(2).Ub=50*ones(1,dim);     % Lower bound of Task 2
            
            g1 = GO_Task1;
            g2 = GO_Task2;
        case 2 % complete intersection with medium similarity, Ackley and Rastrigin
            load('CI_M.mat')
            dim = 30;
            GO_Task1 = zeros(1,dim);
            GO_Task2 = zeros(1,dim);
            Rotation_Task1 = eye(dim);
            Tasks(1).dims = dim;
            Tasks(1).fnc = @(x)Rosenbrock1(x,Rotation_Task1,GO_Task1);
            Tasks(1).Lb=-50*ones(1,dim);
            Tasks(1).Ub=50*ones(1,dim);
            
            Tasks(2).dims = dim;
            Rotation_Task2 = orth(rand(dim));
            Tasks(2).fnc = @(x)Rastrigin1(x,Rotation_Task2,GO_Task2);
            Tasks(2).Lb=-50*ones(1,dim);
            Tasks(2).Ub=50*ones(1,dim);
            
            g1 = GO_Task1;
            g2 = GO_Task2;
        case 3 % complete intersection with low similarity, Ackley and Schwefel
            load('CI_L.mat')
            dim = 30;
            GO_Task1 = zeros(1,dim);
            GO_Task2 = zeros(1,dim);
            Tasks(1).dims = dim;
            Rotation_Task1 = eye(dim);
            Tasks(1).fnc = @(x)Ackley1(x,Rotation_Task1,GO_Task1);
            Tasks(1).Lb=-50*ones(1,dim);
            Tasks(1).Ub=50*ones(1,dim);
            
            Tasks(2).dims = dim;
            Rotation_Task2 = eye(dim);
            Tasks(2).fnc = @(x)Rastrigin1(x, Rotation_Task2, GO_Task2);
            Tasks(2).Lb=-50*ones(1,dim);
            Tasks(2).Ub=50*ones(1,dim);
            
            g1 = GO_Task1;
            g2 = GO_Task2;
        case 4 % partially intersection with high similarity, Rastrigin and Sphere
            load('PI_H.mat')
            dim = 20;
            GO_Task1 = zeros(1,dim);

            Tasks(1).dims = dim;
            Rotation_Task1 = orth(rand(dim));
            Tasks(1).fnc = @(x)Sphere1(x,Rotation_Task1,GO_Task1);
            Tasks(1).Lb=-50*ones(1,dim);
            Tasks(1).Ub=50*ones(1,dim);
            
            dim = 30;
            GO_Task2 = zeros(1,dim);
            Tasks(2).dims = dim;
            Rotation_Task2 = orth(rand(dim));
            Tasks(2).fnc = @(x)Schwefel12(x);
            Tasks(2).Lb=-500*ones(1,dim);
            Tasks(2).Ub=500*ones(1,dim);
            
            g1 = GO_Task1;
            g2 = GO_Task2;
        case 5 % partially intersection with medium similarity, Ackley and Rosenbrock
            load('PI_M.mat')
            dim = 20;
            GO_Task1 = zeros(1,dim);
            Tasks(1).dims = dim;
            Rotation_Task1 = orth(rand(dim));
            Tasks(1).fnc = @(x)Sphere1(x,Rotation_Task1,GO_Task1);
            Tasks(1).Lb=-50*ones(1,dim);
            Tasks(1).Ub=50*ones(1,dim);
            
            dim = 30;
            GO_Task2 = zeros(1,dim);
            Tasks(2).dims = dim;
            Rotation_Task2 = orth(rand(dim));
            Tasks(2).fnc = @(x)Ackley1(x, Rotation_Task2, GO_Task2);
            Tasks(2).Lb=-50*ones(1,dim);
            Tasks(2).Ub=50*ones(1,dim);
            
            g1 = GO_Task1;
            g2 = GO_Task2;
        case 6 % partially intersection with low similarity, Ackley and Weierstrass
            load('PI_L.mat')
            dim = 30;
            GO_Task1 = zeros(1,dim);
            Tasks(1).dims = dim;
            Rotation_Task1 = orth(rand(dim));
            Tasks(1).fnc = @(x)Griewank1(x,Rotation_Task1,GO_Task1);
            Tasks(1).Lb=-50*ones(1,dim);
            Tasks(1).Ub=50*ones(1,dim);
            
            dim = 20;
            GO_Task2 = zeros(1,dim);
            Tasks(2).dims = dim;
            Rotation_Task2 = orth(rand(dim));
            Tasks(2).fnc = @(x)Rastrigin1(x,Rotation_Task2,GO_Task2);
            Tasks(2).Lb=-50*ones(1,dim);
            Tasks(2).Ub=50*ones(1,dim);
            
            g1 = GO_Task1;
            g2 = GO_Task2;
        case 7 % no intersection with high similarity, Rosenbrock and Rastrigin
            load('NI_H.mat')
            dim = 30;
            GO_Task1 = zeros(1,dim);
            Tasks(1).dims = dim;
            Rotation_Task1 = eye(dim);
            Tasks(1).fnc = @(x)Rosenbrock1(x, Rotation_Task1, GO_Task1);
            Tasks(1).Lb=-50*ones(1,dim);
            Tasks(1).Ub=50*ones(1,dim);
            
            GO_Task2 = 28* ones(1,dim);
            Tasks(2).dims = dim;
            Rotation_Task2 = orth(rand(dim));
            Tasks(2).fnc = @(x)Sphere1(x,Rotation_Task2, GO_Task2);
            Tasks(2).Lb=-50*ones(1,dim);
            Tasks(2).Ub=50*ones(1,dim);
            
            g1 = GO_Task1;
            g2 = GO_Task2;
        case 8 % no intersection with medium similarity, Griewank and Weierstrass
            load('NI_M.mat')
            dim = 30;
            GO_Task1 = zeros(1,dim);
            Tasks(1).dims = dim;
            Rotation_Task1 = eye(dim);
            Tasks(1).fnc = @(x)Rosenbrock1(x, Rotation_Task1, GO_Task1);
            Tasks(1).Lb=-50*ones(1,dim);
            Tasks(1).Ub=50*ones(1,dim);
            
            GO_Task2 = 28* ones(1,dim);
            Tasks(2).dims = dim;
            Rotation_Task2 = orth(rand(dim));
            Tasks(2).fnc = @(x)Rastrigin1(x,Rotation_Task2,GO_Task2);
            Tasks(2).Lb=-50*ones(1,dim);
            Tasks(2).Ub=50*ones(1,dim);
            
            g1 = GO_Task1;
            g2 = GO_Task2;
        case 9 % no overlap with low similarity, Rastrigin and Schwefel
            load('NI_L.mat')
            dim = 30;
            GO_Task1 = 28* ones(1,dim);
            Tasks(1).dims = dim;
            Rotation_Task1 = orth(rand(dim));
            Tasks(1).fnc = @(x)Ackley1(x,Rotation_Task1,GO_Task1);
            Tasks(1).Lb=-50*ones(1,dim);
            Tasks(1).Ub=50*ones(1,dim);
            
            GO_Task2 = zeros(1,dim);
            Tasks(2).dims = dim;
            Rotation_Task2 = orth(rand(dim));
            Tasks(2).fnc = @(x)Weierstrass1(x, Rotation_Task2, GO_Task2);
            Tasks(2).Lb=-50*ones(1,dim);
            Tasks(2).Ub=50*ones(1,dim);
            
            g1 = GO_Task1;
            g2 = GO_Task2;
        case 31
            load('NI_L.mat')
            dim = 30;
            GO_Task1 = 0* ones(1,dim);
            Tasks(1).dims = dim;
            Rotation_Task1 = orth(rand(dim));
            Tasks(1).fnc = @(x)Ackley1(x,Rotation_Task1,GO_Task1);
            Tasks(1).Lb=-50*ones(1,dim);
            Tasks(1).Ub=50*ones(1,dim);
            
            GO_Task2 = zeros(1,dim);
            Tasks(2).dims = dim;
            Rotation_Task2 = orth(rand(dim));
            Tasks(2).fnc = @(x)Weierstrass1(x, Rotation_Task2, GO_Task2);
            Tasks(2).Lb=-50*ones(1,dim);
            Tasks(2).Ub=50*ones(1,dim);
            
            GO_Task3 = 0* ones(1,dim);
            Tasks(3).dims = dim;
            Rotation_Task2 = orth(rand(dim));
            Tasks(3).fnc = @(x)Rastrigin1(x,Rotation_Task2,GO_Task3);
            Tasks(3).Lb=-50*ones(1,dim);
            Tasks(3).Ub=50*ones(1,dim);
            
            g1 = GO_Task1;
            g2 = GO_Task2; 
            g3 = GO_Task3;
        case 32
            load('NI_L.mat')
            dim = 30;
            GO_Task1 = [10* ones(1,dim/2), zeros(1,dim/2)];
            Tasks(1).dims = dim;
            Rotation_Task1 = orth(rand(dim));
            Tasks(1).fnc = @(x)Rastrigin1(x,Rotation_Task1,GO_Task1);
            Tasks(1).Lb=-50*ones(1,dim);
            Tasks(1).Ub=50*ones(1,dim);
            
            GO_Task2 = zeros(1,dim);
            Tasks(2).dims = dim;
            Rotation_Task2 = orth(rand(dim));
            Tasks(2).fnc = @(x)Ackley1(x, Rotation_Task2, GO_Task2);%Weierstrass1
            Tasks(2).Lb=-50*ones(1,dim);
            Tasks(2).Ub=50*ones(1,dim);
            
            GO_Task3 = [zeros(1,dim/2),10* ones(1,dim/2)];
            Tasks(3).dims = dim;
            Rotation_Task3 = eye(dim);
            Tasks(3).fnc = @(x)Rosenbrock1(x, Rotation_Task3, GO_Task3);
            Tasks(3).Lb=-50*ones(1,dim);
            Tasks(3).Ub=50*ones(1,dim);
            
            g1 = GO_Task1;
            g2 = GO_Task2; 
            g3 = GO_Task3;
            
        case 101
            dim = 2;
            GO_Task1 = 0*ones(1,dim);
            Tasks(1).dims = dim;
            Tasks(1).fnc = @(x)Rastrigin1(x,eye(dim),GO_Task1);
            Tasks(1).Lb=-50*ones(1,dim);
            Tasks(1).Ub=50*ones(1,dim);
            
            
            GO_Task2 = 0* ones(1,dim);
            Tasks(2).dims = dim;
            Rotation_Task2 = orth(rand(dim));
            Tasks(2).fnc = @(x)Ackley1(x,Rotation_Task2,GO_Task2);
            Tasks(2).Lb=-50*ones(1,dim);
            Tasks(2).Ub=50*ones(1,dim);
            
%             Tasks(2).dims = dim;
%             Tasks(2).fnc = @(x)Schwefel(x);
%             Tasks(2).Lb=-500*ones(1,dim);
%             Tasks(2).Ub=500*ones(1,dim);
%             
%             g1 = GO_Task1;
%             g2 = 420.9687*ones(1,dim);
    end
end