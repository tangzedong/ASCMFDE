function data_MFDE = ASCMFDE(Tasks,pop,gen,selection_process,rmp,p_il,reps,index)
%Directional Distribution Based MFDE function: implementation of MFDE algorithm
%Author: Zedong Tang
%Initial Data: May 21, 2019
%to do:
%JADE
%Position Distribution �� Directional Distribution����Ǩ�Ʋ��Ե�����Ӧѡ��
%improvements:
%1. ʹ��Gedoicisc Flow Transfer��������ֲ�.
%2. History-Based Remix Parameter Coefficients.
%3. Directional Remix Strategy.
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
    for i = 1:pop
        population(i).survive = false;
        population(i).mt = false;
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
    
    F=0.5;
    lb=zeros(1,D_multitask); % ����ȡֵ�½�
    ub=ones(1,D_multitask); % ����ȡֵ�Ͻ�
    pCR=0.5;
    
    Fm = 0.5*ones(1, 2);
    CRm = 0.5*ones(1, 2);
    cm = 0.1;
    goodF = cell(2,1); goodCR = cell(2,1);
    generation=1;
    tmppop1 = population;
    tmppop2 = population;
    trapcount = zeros(1, no_of_tasks);
    while generation < gen
        generation = generation + 1;
        for i = 1:pop
            population(i).mt = false;
        end
        count=1;
        tmppop1 = tmppop2;
        tmppop2 = population;
        if generation > 0.8*gen %% cec19 cancle local search
            p_il = 0.001;
        end
        for i = 1:no_of_tasks
            if ~isempty(goodF{i}) && ~isempty(goodCR{i})
                Fm(i) = (1-cm)*Fm(i) + cm*mean(goodF{i}.^2)/mean(goodF{i});
                CRm(i) = (1-cm)*CRm(i) + cm*mean(goodCR{i});
            end
        end
        
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
            P{i} = P{i} - mean(P{i}, 1);
        end
        %P1 = P1'; P2 = P2';
        %���µ�1�������ӿռ��inter-task trial vector generation����
        %���ӿռ䣬
        %�ŵ�1��ֻ��ǰ�������ɷ�����任������һ�������򻯵����ã�
        %������Ǩ��ʱ���Լ��ٱ����������������Ⱥ��Ӱ��
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
        %�ŵ�2���б�ʽ�⣬����Ҫ�������
        %�ŵ�3����Ϊ�����������ģ�������Ǩ�ƾ���ʱ����Ҫ���棬���ټ�����
        %�ŵ�4��ֻ��Ҫ����һ�α任���󣬾��ܿ��ٵ������һ������ľ���
        %��Ϊʱ�����������ģ�ֻ��ת�þͿ����ˣ���ʡ������Դ��
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
        Flist = zeros(pop,1);
        CRlist = zeros(pop,1);
        
        asf1 = [asf(end),asf(2:end)];
        asf2 = [asf1(end),asf1(2:end)];
        asf3 = [asf2(end),asf2(2:end)];
        
        bsf1 = [bsf(end),bsf(2:end)];
        bsf2 = [bsf1(end),bsf1(2:end)];
        bsf3 = [bsf2(end),bsf2(2:end)];

        sampledStrategy = [];
        mtcount = 0;
        for i = 1 : pop
            [F, CR] = randFCR(1, CRm(population(i).skill_factor), 0.1, Fm(population(i).skill_factor), 0.1);
            x=population(i).rnvec; % ��ȡ����λ��
            
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
            
            mt = false;
            if urmp<=rmp
                p2=ceil(rand()*pop);while(i == p2) p2 = ceil(rand()*pop);end
                p3=ceil(rand()*pop);while(i == p3 && p2 == p3) p3 = ceil(rand()*pop);end
                
                targetskill = population(i).skill_factor;
                sourceskill2 = population(p2).skill_factor;
                sourceskill3 = population(p3).skill_factor;
                
                %�ѽ�ӳ�䵽target subspace
                childsf=1;
                p2pos = (population(p2).rnvec - Pm(sourceskill2,:)) * QQ{sourceskill2};
                p3pos = (population(p3).rnvec - Pm(sourceskill3,:)) * QQ{sourceskill3};
                nn = ceil(3*rand());
                Fr = 0.5;
                %���µ�2�����������ӿռ�Ǩ�ƺ�ԭ�ռ�Ǩ�ƽ������Ч�����ã��������������Ӧ��Ǩ�Ʒ���
                switch nn
                    case 2
                        newpos = Fr*(rand(1, dim).*(p2pos - p3pos)) * M{sourceskill2, targetskill}*QQ{targetskill}';
                    case 3
                        newpos = Fr*(rand(1, dim).*(p2pos - p3pos)) * M{sourceskill2, targetskill}*QQ{targetskill}' + 0.5*Fr*(Pm(sourceskill2,:)- Pm(sourceskill3,:));
                    case 1
                        newpos = Fr*((population(p2).rnvec) - (population(p3).rnvec));
                    otherwise
                        warning('Not existing strategy is called.');
                end
                y=population(p1).rnvec+newpos;%+F*(Pm(population(p1).skill_factor,:)-population(p1).rnvec); % �����м���
                mt = true;
                mtcount = mtcount + 1;
                % y=population(p1).rnvec+F*(population(p2).rnvec - population(p3).rnvec); % �����м���
            else
                p2=asf(ceil(rand()*length(asf)));while(p2 == i)p2=asf(ceil(rand()*length(asf)));end
                p3=asf(ceil(rand()*length(asf)));while(p2 == p3 && p3 == i)p3=asf(ceil(rand()*length(asf)));end
                p4=asf(ceil(rand()*length(asf)));while(p4 == p2 && p4 == p3 && p4 == i)p4=asf(ceil(rand()*length(asf)));end
                
                y=population(p1).rnvec+F*(population(p2).rnvec-population(p3).rnvec);% + F*(Pm(population(p1).skill_factor,:)-population(p1).rnvec); % �����м���
            end
            % ������� Mutation
            %y=population(p1).rnvec+F*(population(p2).rnvec-population(p3).rnvec); % �����м���
            if true%generation > 0.8*gen
                y = mutate(y, length(y), 0.01);%cec19 0.03
            end
            % ��ֹ�м���Խ��
            y=max(y,lb);
            y=min(y,ub);
            
            z=zeros(size(x)); % ��ʼ��һ���¸���
            j0=randi([1,numel(x)]); % ����һ��α���������ѡȡ������ά�ȱ��
            for j=1:numel(x) % ����ÿ��ά��
                if j==j0 || rand<=CR % �����ǰά���Ǵ�����ά�Ȼ����������С�ڽ������
                    z(j)=y(j); % �¸��嵱ǰά��ֵ�����м����Ӧά��ֵ
                else
                    z(j)=x(j); % �¸��嵱ǰά��ֵ���ڵ�ǰ�����Ӧά��ֵ
                end
            end
            
            child(count)=Chromosome();
            child(count)=initialize(child(count),D_multitask);
            child(count).rnvec=z;
            if mt
                child(count).mt = true;
            end
            
            if childsf==0
                child(count).skill_factor=population(i).skill_factor;
            else
                u = rand(1);
                child(count).skill_factor(u<=0.5)=population(i).skill_factor;
                child(count).skill_factor(u>0.5)=2/population(i).skill_factor;
            end
            
            Flist(i) = F;
            CRlist(i) = CR;
            count=count+1;
            
        end
        for i = 1 : pop
            [child(i),calls_per_individual(i)] = evaluate(child(i),Tasks,p_il,no_of_tasks,options);
        end
        fnceval_calls(rep)=fnceval_calls(rep) + sum(calls_per_individual);
        TotalEvaluations(rep,generation)=fnceval_calls(rep);
        
        intpopulation(1:pop)=population;
        intpopulation(pop+1:2*pop)=child;
        tmpFlist(1:pop) = 0;
        tmpFlist(pop+1:2*pop) = Flist;
        tmpCRlist(1:pop) = 0;
        tmpCRlist(pop+1:2*pop) = CRlist;
        factorial_cost=zeros(1,2*pop);
        flag = zeros(no_of_tasks, 1);
        for i = 1:no_of_tasks
            for j = 1:2*pop
                factorial_cost(j)=intpopulation(j).factorial_costs(i);
            end
            [xxx,y]=sort(factorial_cost);
            intpopulation=intpopulation(y);
            tmpFlist = tmpFlist(y);
            tmpCRlist = tmpCRlist(y);
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
            tmpFlist = tmpFlist(y);
            tmpCRlist = tmpCRlist(y);
            population=intpopulation(1:pop);
            Flist = tmpFlist(1:pop);
            CRlist = tmpCRlist(1:pop);
            for i = 1:pop
                if Flist(i) > 0
                    goodF{population(i).skill_factor} = [goodF{population(i).skill_factor} Flist(i)];
                end
                if CRlist(i) > 0
                    goodCR{population(i).skill_factor} = [goodCR{population(i).skill_factor} CRlist(i)];
                end
            end
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
        survivecounter = 0;
        for i = 1:pop
            if population(i).mt
                survivecounter = survivecounter + 1;
            end
        end
        successrate(generation) = survivecounter/mtcount;
        %disp(['MFDE Generation = ', num2str(generation), ' best factorial costs = ', num2str(bestobj), ' taskid = ', num2str(survivecounter/mtcount)]);
    end
end

dlmwrite(['MTSOO_P',num2str(index),'.txt'],bestFncErrorValue,'precision',6);
data_MFDE.wall_clock_time=toc;
data_MFDE.EvBestFitness=EvBestFitness;
data_MFDE.bestInd_data=bestInd_data;
data_MFDE.TotalEvaluations=TotalEvaluations;
data_MFDE.crsuccessrate = successrate;
end