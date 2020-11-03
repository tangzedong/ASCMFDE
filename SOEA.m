function data_SOO = SOEA(Task,pop,gen,selection_process,p_il,reps)
%SOEA function: implementation of SOEA algorithm
    clc    
    tic         
    if mod(pop,2) ~= 0
        pop = pop + 1;
    end   
    D = Task.dims;
    options = optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton','MaxIter',2);  % settings for individual learning
    
    fnceval_calls = zeros(1,reps); 
    calls_per_individual=zeros(1,pop);
    EvBestFitness = zeros(reps,gen);    % best fitness found
    TotalEvaluations=zeros(reps,gen);   % total number of task evaluations so fer
    successcounter = 0;
    for rep = 1:reps   
        disp(rep)
        for i = 1 : pop
            population(i) = Chromosome();
            population(i) = initialize(population(i),D);
        end
        for i = 1 : pop
            [population(i),calls_per_individual(i)] = evaluate_SOO(population(i),Task,p_il,options);
        end

        fnceval_calls(rep)=fnceval_calls(rep) + sum(calls_per_individual);
        TotalEvaluations(rep,1)=fnceval_calls(rep);    
        bestobj=min([population.factorial_costs]);
        EvBestFitness(rep,1) = bestobj;    

%         VarSize=[1 D];   % Decision Variables Matrix Size
%         beta_min=0.2;   % Lower Bound of Scaling Factor
%         beta_max=0.8;   % Upper Bound of Scaling Factor
        lb=zeros(1,D); % ����ȡֵ�½�
        ub=ones(1,D); % ����ȡֵ�Ͻ�
        pCR=0.7;
        F=0.5;
        generation=1;
        successflag = false;
        while generation < gen
            generation = generation + 1;
            count=1;
            for i = 1 : pop
                x=population(i).rnvec; % ��ȡ����λ��
                A = randperm(pop);
              
                A(A==i)=[]; % ��ǰ��������λ���ڿգ����������м���ʱ��ǰ���岻���룩
                p1=A(1);
                p2=A(2);
                p3=A(3);
                % ������� Mutation
%                 beta=unifrnd(beta_min,beta_max,VarSize); % ���������������
                y=population(p1).rnvec+F*(population(p2).rnvec-population(p3).rnvec); % �����м���
                % ��ֹ�м���Խ��
                y=max(y,lb);
                y=min(y,ub);
            
                z=zeros(size(x)); % ��ʼ��һ���¸���
                j0=randi([1,numel(x)]); % ����һ��α���������ѡȡ������ά�ȱ��
                for j=1:numel(x) % ����ÿ��ά��
                    if j==j0 || rand<=pCR % �����ǰά���Ǵ�����ά�Ȼ����������С�ڽ������
                    z(j)=y(j); % �¸��嵱ǰά��ֵ�����м����Ӧά��ֵ
                    else
                    z(j)=x(j); % �¸��嵱ǰά��ֵ���ڵ�ǰ�����Ӧά��ֵ
                    end
                end
                
                child(count)=Chromosome();
                child(count).rnvec=z;
                
                count=count+1;
            end        
            for i = 1 : pop            
                [child(i),calls_per_individual(i)] = evaluate_SOO(child(i),Task,p_il,options);           
            end      

            fnceval_calls(rep)=fnceval_calls(rep) + sum(calls_per_individual);
            TotalEvaluations(rep,generation)=fnceval_calls(rep);
            
            intpopulation(1:pop)=population;
            intpopulation(pop+1:2*pop)=child;
            [xxx,y]=sort([intpopulation.factorial_costs]);
            intpopulation=intpopulation(y);
            for i = 1:2*pop
                intpopulation(i).scalar_fitness=1/i;
            end
            if intpopulation(1).factorial_costs<=bestobj
                bestobj=intpopulation(1).factorial_costs;
                bestInd_data(rep)=intpopulation(1);
            end
            EvBestFitness(rep,generation)=bestobj;

            if strcmp(selection_process,'elitist')
                [xxx,y]=sort(-[intpopulation.scalar_fitness]);
                intpopulation=intpopulation(y);
                population=intpopulation(1:pop);            
            elseif strcmp(selection_process,'roulette wheel')
                for i=1:pop
                    population(i)=intpopulation(RouletteWheelSelection([intpopulation.scalar_fitness]));
                end    
            end
            if mod(generation, 100) == 0
            disp(['SOO Generation ', num2str(generation), ' best objective = ', num2str(bestobj)])
            end
            if ~successflag && bestobj<10e-5
                 successflag = true;
            end
             if all(successflag)
                 break
             end
        end    
        if successflag
            successcounter = successcounter + 1;
        end
    end
    data_SOO.wall_clock_time=toc;
    data_SOO.EvBestFitness=EvBestFitness;
    data_SOO.bestInd_data=bestInd_data;
    data_SOO.TotalEvaluations=TotalEvaluations;
    data_SOO.Success = successcounter;
end