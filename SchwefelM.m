function obj = SchwefelM(var,M,b)
%SCHWEFEL function
%   - var: design variable vector
    dim = length(var);
    var = var - b;
    var = var;
    var = M*var;
    var = var + 4.209687462275036e+002;
    sum = 0;
    for i = 1: dim
        if var(i)>500
            sum = sum + (500 - mod(var(i),500)) * sin(sqrt(abs(500 - mod(var(i),500))));
            tmp = (var(i) - 500)/100;
            sum = sum - tmp*tmp/dim;
        elseif var(i)<-500
            sum = sum + (-500 +mod(var(i),500)) * sin(sqrt(abs(-500 + mod(var(i),500))));
            tmp = (var(i) + 500)/100;
            sum = sum - tmp*tmp/dim;
        else
            sum = sum + var(i)*sin(sqrt(abs(var(i))));
        end
    end
    
    obj = 418.9829*dim-sum;
    obj = obj + 1100.0;
end

