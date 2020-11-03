function data_MFDE = ASCMFDE(Tasks,pop,gen,selection_process,rmp,p_il,reps,index)
%Directional Distribution Based MFDE function: implementation of MFDE algorithm
%Author: Zedong Tang
%Initial Data: May 21, 2019
clc
tic
if mod(pop,2) ~= 0
    pop = pop + 1;
end
no_of_tasks=length(Tasks);
if no_of_tasks <= 1
    error('At least 2 tasks required for MFDE');
end
D=zeros(1,no_of_tasks);
for i=1:no_of_tasks
    D(i)=Tasks(i).dims;
end
D_multitask=max(D);
options = optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton','MaxIter',2);  % settings for individual learning

fnceval_calls = zeros(1,reps);
calls_per_individual=zeros(1,pop);
EvBestFitness = zeros(no_of_tasks*reps,gen);    % best fitness found
TotalEvaluations=zeros(reps,gen);               % total number of task evaluations so fer
bestobj=Inf(1,no_of_tasks);
bestFncErrorValue = zeros(100,60);

for rep = 1:reps
    disp(rep)
    for i = 1 : pop
        population(i) = Chromosome();
        population(i) = initialize(population(i),D_multitask);
        population(i).skill_factor=0;
    end
    for i = 1 : pop
        [population(i),calls_per_individual(i)] = evaluate(population(i),Tasks,p_il,no_of_tasks,options);
    end
    
    fnceval_calls(rep)=fnceval_calls(rep) + sum(calls_per_individual);
    TotalEvaluations(rep,1)=fnceval_calls(rep);
    
    factorial_cost=zeros(1,pop);
    for i = 1:no_of_tasks
        for j = 1:pop
            factorial_cost(j)=population(j).factorial_costs(i);
        end
        [xxx,y]=sort(factorial_cost);
        population=population(y);
        for j=1:pop
            population(j).factorial_ranks(i)=j;
        end
        bestobj(i)=population(1).factorial_costs(i);
        EvBestFitness(i+2*(rep-1),1)=bestobj(i);
        bestInd_data(rep,i)=population(1);
    end
    for i=1:pop
        [xxx,yyy]=min(population(i).factorial_ranks);
        x=find(population(i).factorial_ranks == xxx);
        equivalent_skills=length(x);
        if equivalent_skills>1
            population(i).skill_factor=x(1+round((equivalent_skills-1)*rand(1)));
            tmp=population(i).factorial_costs(population(i).skill_factor);
            population(i).factorial_costs(1:no_of_tasks)=inf;
            population(i).factorial_costs(population(i).skill_factor)=tmp;
        else
            population(i).skill_factor=yyy;
            tmp=population(i).factorial_costs(population(i).skill_factor);
            population(i).factorial_costs(1:no_of_tasks)=inf;
            population(i).factorial_costs(population(i).skill_factor)=tmp;
        end
    end
    
    
    lb=zeros(1,D_multitask); % 参数取值下界
    ub=ones(1,D_multitask); % 参数取值上界
    F=0.5;
    CR=0.7;
    
    cm = 0.1;
    generation=1;
    tmppop1 = population;
    tmppop2 = population;
    trapcount = zeros(1, no_of_tasks);
    while generation < gen
        generation = generation + 1;
        count=1;
        tmppop1 = tmppop2;
        tmppop2 = population;
        
        asf=[];bsf=[];
        for j = 1:pop
            if population(j).skill_factor==1
                asf=[asf,j];
            else
                bsf=[bsf,j];
            end
        end
        group={asf,bsf};
        
        P = cell(1, no_of_tasks);
        for i = 1:no_of_tasks
            P{i} = [];
        end
        for j = 1:pop
            P{population(j).skill_factor} = [P{population(j).skill_factor};population(j).rnvec];
        end
        for j = 1:pop
            P{population(j).skill_factor} = [P{population(j).skill_factor};tmppop1(j).rnvec];
        end
        for i = 1:no_of_tasks
            Pm(i,:) = mean(P{i}, 1);
            P{i} = P{i} - repmat(mean(P{i}, 1),[length(P{i}),1]);
        end
        
        dim = Inf;
        for i = 1:no_of_tasks
            [popNew{i}, pcaModel{i}] = ftProc_pca_tr(P{i},[], struct('pcaCoef',0));
            if size(pcaModel{i}.W_prj, 2) < dim
                dim = size(pcaModel{i}.W_prj, 2);
            end
        end
        dim = floor(0.8*dim);
        for i = 1:no_of_tasks
            QQ{i} = pcaModel{i}.W_prj(:, 1:dim);
        end
        %
        M = cell(no_of_tasks, no_of_tasks);
        for i = 1:no_of_tasks
            for j = i:no_of_tasks
                if i == j
                    M{i,j} = eye(dim);
                else
                    M{i,j} = QQ{i}'*QQ{j};
                    M{j,i} = QQ{j}'*QQ{i};
                end
            end
        end
        
        for i = 1:no_of_tasks
            goodF{i} = [];
            goodCR{i} = [];
        end
        
        asf1 = [asf(end),asf(2:end)];
        asf2 = [asf1(end),asf1(2:end)];
        asf3 = [asf2(end),asf2(2:end)];
        
        bsf1 = [bsf(end),bsf(2:end)];
        bsf2 = [bsf1(end),bsf1(2:end)];
        bsf3 = [bsf2(end),bsf2(2:end)];
        
        for i = 1 : pop
            x=population(i).rnvec; % 提取个体位置
            
            asf=cell2mat(group(population(i).skill_factor));
            bsf=cell2mat(group(2/population(i).skill_factor));
            
            B = randperm(length(asf));
            C = randperm(length(bsf));
            asf=asf(B);
            bsf=bsf(C);
            
            childsf=0;
            
            for j=1:length(asf)
                if asf(j)==i
                    asf(j)=[];
                    break;
                end
            end
            
            p1=asf(1);
            
            urmp = rand(1);
            if urmp<=rmp
                p2=ceil(rand()*pop);while(i == p2) p2 = ceil(rand()*pop);end
                p3=ceil(rand()*pop);while(i == p3 && p2 == p3) p3 = ceil(rand()*pop);end
                
                targetskill = population(i).skill_factor;
                sourceskill2 = population(p2).skill_factor;
                sourceskill3 = population(p3).skill_factor;
                
                %把解映射到target subspace
                childsf=1;
                p2pos = (population(p2).rnvec - Pm(sourceskill2,:)) * QQ{sourceskill2};
                p3pos = (population(p3).rnvec - Pm(sourceskill3,:)) * QQ{sourceskill3};
                nn = ceil(3*rand());
                F = 0.5;
                CR = 0.7;

                switch nn
                    case 2
                        newpos = F*(rand(1, dim).*(p2pos - p3pos)) * M{sourceskill2, targetskill}*QQ{targetskill}';
                    case 3
                        newpos = F*(rand(1, dim).*(p2pos - p3pos)) * M{sourceskill2, targetskill}*QQ{targetskill}' + 0.5*F*(Pm(sourceskill2,:)- Pm(sourceskill3,:));
                    case 1
                        newpos = F*((population(p2).rnvec) - (population(p3).rnvec));
                    otherwise
                        warning('Not existing strategy is called.');
                end
                y=population(p1).rnvec+newpos;
            else
                p2=asf(ceil(rand()*length(asf)));while(p2 == i)p2=asf(ceil(rand()*length(asf)));end
                p3=asf(ceil(rand()*length(asf)));while(p2 == p3 && p3 == i)p3=asf(ceil(rand()*length(asf)));end
                p4=asf(ceil(rand()*length(asf)));while(p4 == p2 && p4 == p3 && p4 == i)p4=asf(ceil(rand()*length(asf)));end
                
                y=population(p1).rnvec+F*(population(p2).rnvec-population(p3).rnvec);% + F*(Pm(population(p1).skill_factor,:)-population(p1).rnvec); % 产生中间体
            end
            % 变异操作 Mutation
            y = mutate(y, length(y), 0.01);
            % 防止中间体越界
            y=max(y,lb);
            y=min(y,ub);
            
            z=zeros(size(x)); % 初始化一个新个体
            j0=randi([1,numel(x)]); % 产生一个伪随机数，即选取待交换维度编号
            for j=1:numel(x) % 遍历每个维度
                if j==j0 || rand<=CR % 如果当前维度是待交换维度或者随机概率小于交叉概率
                    z(j)=y(j); % 新个体当前维度值等于中间体对应维度值
                else
                    z(j)=x(j); % 新个体当前维度值等于当前个体对应维度值
                end
            end
            
            child(count)=Chromosome();
            child(count).rnvec=z;
            
            
            if childsf==0
                child(count).skill_factor=population(i).skill_factor;
            else
                u = rand(1);
                child(count).skill_factor(u<=0.5)=population(i).skill_factor;
                child(count).skill_factor(u>0.5)=2/population(i).skill_factor;
            end
            count=count+1;
            
        end
        for i = 1 : pop
            [child(i),calls_per_individual(i)] = evaluate(child(i),Tasks,p_il,no_of_tasks,options);
        end
        fnceval_calls(rep)=fnceval_calls(rep) + sum(calls_per_individual);
        TotalEvaluations(rep,generation)=fnceval_calls(rep);
        
        intpopulation(1:pop)=population;
        intpopulation(pop+1:2*pop)=child;

        factorial_cost=zeros(1,2*pop);
        flag = zeros(no_of_tasks, 1);
        for i = 1:no_of_tasks
            for j = 1:2*pop
                factorial_cost(j)=intpopulation(j).factorial_costs(i);
            end
            [xxx,y]=sort(factorial_cost);
            intpopulation=intpopulation(y);
            for j=1:2*pop
                intpopulation(j).factorial_ranks(i)=j;
            end
            if intpopulation(1).factorial_costs(i)<bestobj(i)
                bestobj(i)=intpopulation(1).factorial_costs(i);
                bestInd_data(rep,i)=intpopulation(1);
                flag(i) = 1;
            end
            EvBestFitness(i+2*(rep-1),generation)=bestobj(i);
            
            if mod(fnceval_calls(rep),3000)==0
                bestFncErrorValue(fnceval_calls(rep)/3000,1)=fnceval_calls(rep);
                bestFncErrorValue(fnceval_calls(rep)/3000,i+2*(rep-1)+1)=bestobj(i);
            end
        end
        for i = 1:no_of_tasks
            if flag(i) == 1
                trapcount(i) = 0;
            else
                trapcount(i) = trapcount(i) + 1;
            end
            if trapcount(i) > 20
                trapcount(i) = 0;
                nn = 0;
                for j = pop:-1:1
                    if population(j).skill_factor == i
                        population(j).rnvec = rand(1, D_multitask);
                        [population(j),~] = evaluate(population(j),Tasks,0,no_of_tasks,options);
                        nn = nn + 1;
                    end
                    if nn > 10
                        break;
                    end
                end
            end
        end
        for i=1:2*pop
            [xxx,yyy]=min(intpopulation(i).factorial_ranks);
            intpopulation(i).skill_factor=yyy;
            intpopulation(i).scalar_fitness=1/xxx;
        end
        
        if strcmp(selection_process,'elitist')
            [xxx,y]=sort(-[intpopulation.scalar_fitness]);
            intpopulation=intpopulation(y);
            population=intpopulation(1:pop);
        elseif strcmp(selection_process,'roulette wheel')
            for i=1:no_of_tasks
                skill_group(i).individuals=intpopulation([intpopulation.skill_factor]==i);
            end
            count=0;
            while count<pop
                count=count+1;
                skill=mod(count,no_of_tasks)+1;
                population(count)=skill_group(skill).individuals(RouletteWheelSelection([skill_group(skill).individuals.scalar_fitness]));
            end
        end
        disp(['MFDE Generation = ', num2str(generation), ' best factorial costs = ', num2str(bestobj), ' taskid = ', num2str(index)]);
    end
end

data_MFDE.wall_clock_time=toc;
data_MFDE.EvBestFitness=EvBestFitness;
data_MFDE.bestInd_data=bestInd_data;
data_MFDE.TotalEvaluations=TotalEvaluations;
end