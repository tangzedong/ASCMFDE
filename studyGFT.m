mu = [2,3];
sigma = [1,1.5;1.5,3];
% rng default  % For reproducibility
r1 = mvnrnd(mu,sigma,800);
plot(r1(:,1), r1(:, 2), '+');

mu = [-10,-10];
sigma = [3,0;0,0.1];
%rng default  % For reproducibility
r2 = mvnrnd(mu,sigma,800);
hold on
plot(r2(:,1), r2(:, 2), 'o');

Ps = pca(r1 - mean(r1, 1));
Pt = pca(r2 - mean(r2, 1));

[V1, V2, theta, dim] = GFT(Ps, Pt(:, 1));
rr = rand();
A = generateGFT(Ps, V1, V2, theta, dim, rr);

for i = 1:800
    p1 = r1(ceil(rand*size(r1, 1)), :);
    p2 = r2(ceil(rand*size(r2, 1)), :);
    p3 = r1(ceil(rand*size(r1, 1)), :);
    p11 = (p1 - mean(r1, 1))*A;
    p22 = (p2 - mean(r2, 1))*A;
    p33 = (p3 - mean(r1, 1))*A;
    newpos = p3 + (0.5*(p22 - p11))*A';
    offspring(i, :) = newpos;
end
hold on;
plot(offspring(:, 1), offspring(:,2), '.');