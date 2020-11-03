function Tasks = benchmark19(index)
%BENCHMARK function
%   Input
%   - index: the index number of problem set
%
%   Output:
%   - Tasks: benchmark problem set

    switch(index)
        case 1
            dim = 50;
            Tasks(1).dims = dim;   
            Tasks(1).fnc = 7;    %f7 in CEC2014 test suite
            Tasks(1).Lb=-100*ones(1,dim);   
            Tasks(1).Ub=100*ones(1,dim);    

            Tasks(2).dims = dim;   
            Tasks(2).fnc = 8;   %f8 in CEC2014 test suite
            Tasks(2).Lb=-100*ones(1,dim);   
            Tasks(2).Ub=100*ones(1,dim);    
            
        case 2             
            dim = 50;
            Tasks(1).dims = dim;   
            Tasks(1).fnc = 7;
            Tasks(1).Lb=-100*ones(1,dim);   
            Tasks(1).Ub=100*ones(1,dim);    

            Tasks(2).dims = dim;   
            Tasks(2).fnc = 12;
            Tasks(2).Lb=-100*ones(1,dim);   
            Tasks(2).Ub=100*ones(1,dim);                
     
        case 3 
            dim = 50;
            Tasks(1).dims = dim;   
            Tasks(1).fnc = 7;
            Tasks(1).Lb=-100*ones(1,dim);   
            Tasks(1).Ub=100*ones(1,dim);    

            Tasks(2).dims = dim;   
            Tasks(2).fnc = 15;
            Tasks(2).Lb=-100*ones(1,dim);   
            Tasks(2).Ub=100*ones(1,dim);                

        case 4 
            dim = 50;
            Tasks(1).dims = dim;   
            Tasks(1).fnc = 8;
            Tasks(1).Lb=-100*ones(1,dim);   
            Tasks(1).Ub=100*ones(1,dim);    

            Tasks(2).dims = dim;   
            Tasks(2).fnc = 15;
            Tasks(2).Lb=-100*ones(1,dim);   
            Tasks(2).Ub=100*ones(1,dim);                       

        case 5
            dim = 50;
            Tasks(1).dims = dim;   
            Tasks(1).fnc = 9;
            Tasks(1).Lb=-100*ones(1,dim);   
            Tasks(1).Ub=100*ones(1,dim);    

            Tasks(2).dims = dim;   
            Tasks(2).fnc = 10;
            Tasks(2).Lb=-100*ones(1,dim);   
            Tasks(2).Ub=100*ones(1,dim);                       

        case 6 
            dim = 50;
            Tasks(1).dims = dim;   
            Tasks(1).fnc = 9;
            Tasks(1).Lb=-100*ones(1,dim);   
            Tasks(1).Ub=100*ones(1,dim);    

            Tasks(2).dims = dim;   
            Tasks(2).fnc = 11;
            Tasks(2).Lb=-1000*ones(1,dim);   
            Tasks(2).Ub=1000*ones(1,dim);           
        case 7 
            dim = 50;
            Tasks(1).dims = dim;   
            Tasks(1).fnc = 10;
            Tasks(1).Lb=-100*ones(1,dim);   
            Tasks(1).Ub=100*ones(1,dim);    

            Tasks(2).dims = dim;   
            Tasks(2).fnc = 16;
            Tasks(2).Lb=-100*ones(1,dim);   
            Tasks(2).Ub=100*ones(1,dim);        
        case 8 
            dim = 50;
            Tasks(1).dims = dim;   
            Tasks(1).fnc = 11;
            Tasks(1).Lb=-100*ones(1,dim);   
            Tasks(1).Ub=100*ones(1,dim);    

            Tasks(2).dims = dim;   
            Tasks(2).fnc = 13;
            Tasks(2).Lb=-100*ones(1,dim);   
            Tasks(2).Ub=100*ones(1,dim);        
        case 9 
            dim = 50;
            Tasks(1).dims = dim;   
            Tasks(1).fnc = 13;
            Tasks(1).Lb=-100*ones(1,dim);   
            Tasks(1).Ub=100*ones(1,dim);    

            Tasks(2).dims = dim;   
            Tasks(2).fnc = 14;
            Tasks(2).Lb=-100*ones(1,dim);   
            Tasks(2).Ub=100*ones(1,dim);                    
        case 10 
            dim = 50;
            Tasks(1).dims = dim;   
            Tasks(1).fnc = 14;
            Tasks(1).Lb=-100*ones(1,dim);   
            Tasks(1).Ub=100*ones(1,dim);    

            Tasks(2).dims = dim;   
            Tasks(2).fnc = 16;
            Tasks(2).Lb=-100*ones(1,dim);   
            Tasks(2).Ub=100*ones(1,dim);        
    end
    fnc1 = Tasks(1).fnc;
    fnc2 = Tasks(2).fnc;
    M1 = load(sprintf('./input_data/M_%d_D%d.txt', Tasks(1).fnc, Tasks(1).dims));
    M2 = load(sprintf('./input_data/M_%d_D%d.txt', Tasks(2).fnc, Tasks(2).dims));
    Shift1 =  load(sprintf('./input_data/shift_data_%d.txt', Tasks(1).fnc));
    Shift2 =  load(sprintf('./input_data/shift_data_%d.txt', Tasks(2).fnc));
    Shuffle1 =  load(sprintf('./input_data/shuffle_data_%d_D%d.txt', Tasks(1).fnc, Tasks(1).dims));
    Shuffle2 =  load(sprintf('./input_data/shuffle_data_%d_D%d.txt', Tasks(2).fnc, Tasks(2).dims));
    if Tasks(1).fnc == 11
        Tasks(1).fnc = @(x)SchwefelM(x', M2(1:Tasks(2).dims,1:Tasks(2).dims), Shift2(1:Tasks(2).dims)');
    else
        Tasks(1).fnc = @(x)mycec14_func(x', fnc1, M1', Shift1', Shuffle1);
    end
    if Tasks(2).fnc == 11
        Tasks(2).fnc = @(x)SchwefelM(x', M2(1:Tasks(2).dims,1:Tasks(2).dims), Shift2(1:Tasks(2).dims)');
    else
        Tasks(2).fnc = @(x)mycec14_func(x', fnc2, M2', Shift2', Shuffle2);
    end
end