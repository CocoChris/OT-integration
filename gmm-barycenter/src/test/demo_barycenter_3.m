clear;clc;
addpath('../barycenter/');
%%
N=30; % number of instances
m=5; % number of states
d=3; % dimension
Tau = [1,0,0;0,1.5,0;0,0,2];
df =10;

mu = [2,3,4];
sigma = ones(3,3);

supp=zeros(d+d*d, m*N);
w = zeros(1,m*N);
for i=1:N
    for j=1:m
        r = mvnrnd(mu,sigma,1)'; %生成高斯分布的数据（1个点）
        S = iwishrnd(Tau,df)*(df-2-1);
        supp(:,(i-1)*m+j)=[r;S(:)];       
    end
    t=rand(m,1);
    w((i-1)*m+1:i*m) = t/sum(t);    
end
stride=m*ones(1,N);
instanceW=rand(1,N); % not necessarily sum2one
%% compute GMM barycenter with B-ADMM method allowing the change of component weights
% initiate centroid from an instance
c0.supp=supp(:, 1:m);
c0.w=w(1:m);
% set the number of iterations and badmm_rho (no need to change)
options.badmm_max_iters=2000;
options.badmm_rho=10;
% compute the centroid
% [c, X]=centroid_sphBregman_GMM(stride, instanceW, supp, w, c0, options);
[c2, X2]=centroid_sphBregman_GMM_test(stride, instanceW, supp, w, c0, options);
