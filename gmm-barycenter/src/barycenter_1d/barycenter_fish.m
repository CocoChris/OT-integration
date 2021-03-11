clear;clc;
%addpath('barycenter/');

load fishmcmc.mat

%%
N=10000; % number of instances
m=5; % number of states
d=1; % dimension

supp=zeros(d+d*d, m*N);
ww = zeros(1,m*N);
for i=1:N
    for j=1:m
        r = mean(i,j);
        S = var(i,j);
        supp(:,(i-1)*m+j)=[r;S(:)];       
    end

    ww((i-1)*m+1:i*m) = w(i,:);    
end

stride=m*ones(1,N);
instanceW=rand(1,N); % not necessarily sum2one
%% compute GMM barycenter with B-ADMM method allowing the change of component weights
% initiate centroid from an instance
c0.supp=supp(:, 1:m);
c0.w=ww(1:m);
% set the number of iterations and badmm_rho (no need to change)
options.badmm_max_iters=2000;
options.badmm_rho=10;
options.tau = 50;
% compute the centroid
[c, X]=centroid_sphBregman_GMM(stride, instanceW, supp, ww, c0, options);
