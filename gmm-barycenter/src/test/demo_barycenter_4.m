clear;clc;
addpath('../barycenter/');

N = 30;
m = 5;

stride=m*ones(1,N);
prefix = '~/Workspace/R/barycenter/txt/';
supp = importdata(strcat(prefix, 'supp.txt'));
instanceW = importdata(strcat(prefix, 'instanceW.txt'));
w = importdata(strcat(prefix, 'w.txt'));
supp = supp.data;
instanceW = instanceW.data;
w = w.data;
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

