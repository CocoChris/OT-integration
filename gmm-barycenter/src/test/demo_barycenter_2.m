clear;clc;
addpath('../barycenter/');
%%
N=30; % number of instances
%m=5; % number of states
d=2; % dimension
Tau = [1 .5; .5 2];
df =10;

mlist = zeros(1,N);
for i=1:N
    mlist(i) = floor(unifrnd(4,7));
end

num = sum(mlist);
fprintf('\t %d', num);

mu = [2,3];
sigma = [1,1.5;1.5,3];

supp=zeros(d+d*d, num);
w = zeros(1,num);

k = 1;
for i=1:N
    m = mlist(i);
    t=rand(m,1);
    w(k:k+m-1) = t/sum(t);
    for j=1:m
        r = mvnrnd(mu,sigma,1)'; %生成高斯分布的数据（1个点）
        S = iwishrnd(Tau,df)*(df-2-1);
        %supp(:,(i-1)*m+j)=[r;S(:)]; 
        supp(:,k)=[r;S(:)];
        k = k+1;
    end
end
%stride=m*ones(1,N);
stride = mlist;
instanceW=rand(1,N); % not necessarily sum2one
%% compute GMM barycenter with B-ADMM method allowing the change of component weights
% initiate centroid from an instance
m0 = 3;
c0.supp=supp(:, 1:m0);
c0.w=w(1:m0);
% set the number of iterations and badmm_rho (no need to change)
options.badmm_max_iters=2000;
options.badmm_rho=10;
% compute the centroid
[c, X]=centroid_sphBregman_GMM(stride, instanceW, supp, w, c0, options);
