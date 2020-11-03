function object=mutate(p,dim,mum)
    rnvec_temp=p;
%     for i=1:dim
%         if rand(1)<1/dim
%             u=rand(1);
%             if u <= 0.5
%                 del=(2*u)^(1/(1+mum)) - 1;
%                 rnvec_temp(i)=p(i) + del*(p(i));
%             else
%                 del= 1 - (2*(1-u))^(1/(1+mum));
%                 rnvec_temp(i)=p(i) + del*(1-p(i));
%             end
%         end
%     end  
    for i = 1:dim 
        if rand() < mum
            rnvec_temp(i) = rand();
        end
    end
    rnvec_temp(rnvec_temp>1) = rand();
    rnvec_temp(rnvec_temp<0) = rand();
    object = rnvec_temp;          
end  