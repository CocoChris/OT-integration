clear;clc;
addpath('./barycenter/');

dataset_num=2;
cluster_num=[6,9];
d=20;

no_dims = 2;
initial_dims = 10;
perplexity = 30;

supp=zeros(d+d*d,sum(cluster_num));

% replace with your actual pbmc path
prefix='~/Workspace/code/OT-integration/data/pbmc/';

k=1;
cluster_weights=zeros(1,sum(cluster_num));
data_before_ot = [];
index_list=zeros(1,dataset_num);
techs = char('seqwell', '10x');

for i=1:dataset_num
    tech = techs(i,:);
    cluster_label_filename=strcat(prefix, 'clusters_', tech, '.txt');
    cca_filename=strcat(prefix, 'cca/', tech, '.txt');
    
    cluster_labels=importdata(cluster_label_filename);
%     cluster_labels=raw_data(:,1);
    cca_data=importdata(cca_filename);
    cca_data=cca_data.data(:,1:d);
    
    all_data=[cluster_labels cca_data];
    sample_num=size(cca_data,1);
    
    cca_data_group=cell(1,cluster_num(i));

    for j=1:cluster_num(i)
      cca_data_group{j}=all_data(all_data(:,1)==j,2:end);
    end
    
    for j=1:cluster_num(i)
       cca_data=cca_data_group{j};
       m=size(cca_data,1);
       cluster_weights(k)=m/sample_num;
       
       % compute mean and variance
       data_mean=mean(cca_data, 1);
       if m == 1
           cca_data=[cca_data;cca_data]
       end
       data_cov=cov(cca_data);
       supp(:,k)=[data_mean.';data_cov(:)];
       
       k=k+1;
    end
    
    data_before_ot = [data_before_ot;all_data];
    index_list(i) = sample_num;
end

stride = zeros(1,dataset_num);
for i=1:dataset_num
    stride(i)=cluster_num(i);
end
% instanceW=rand(1,dataset_num); % not necessarily sum2one
instanceW = [1,1];
% compute GMM barycenter with B-ADMM method allowing the change of component weights
% initiate centroid from an instance
m0=6;
c0.supp=supp(:,1:m0);
c0.w=cluster_weights(1:m0);
% set the number of iterations and badmm_rho (no need to change)
options.badmm_max_iters=2000;
options.badmm_rho=10;
% compute the centroid
[c, X]=centroid_sphBregman_GMM_test(stride, instanceW, supp, cluster_weights, c0, options);

index_list = cumsum(index_list);
[data_after_ot, dist_mean] = re_cluster(c.supp, data_before_ot(:,2:d+1));

% save gamma matrix
fid=fopen(strcat(prefix, 'gamma.txt'),'w');
data_output=real(X);
[r,c]=size(data_output);            
 for i=1:r
  for j=1:c
    fprintf(fid,'%f\t',data_output(i,j));
  end
    fprintf(fid,'\n');
 end
fclose(fid);


data_before_tsne = [data_before_ot;data_after_ot];

dataset_label = [];
dataset_label(1:index_list(1)) = 1;
dataset_label((index_list(1)+1):index_list(2)) = 2;

raw_data_seqwell = importdata(strcat(prefix, 'pca/seqwell.txt'));
raw_data_seqwell = raw_data_seqwell.data;
raw_data_10x = importdata(strcat(prefix, 'pca/10x.txt'));
raw_data_10x = raw_data_10x.data;

raw_data_union = [raw_data_seqwell; raw_data_10x];

seurat_clusters = importdata(strcat(prefix, '/clusters/seurat3_clusters.txt'));
seurat_clusters = seurat_clusters.data(:,2);

liger_clusters = importdata(strcat(prefix, '/clusters/liger_clusters.txt'));
liger_clusters = liger_clusters.data(:,2);

seurat_tsne = importdata('/Users/coco/Workspace/code/OT-integration/data/pbmc/tsne/pbmc_seurat3.txt');
seurat_tsne = seurat_tsne.data;

liger_tsne = importdata('/Users/coco/Workspace/code/OT-integration/data/pbmc/tsne/pbmc_liger.txt');
liger_tsne = liger_tsne.data;


raw_data_tsne = tsne(raw_data_union);

% Run t-SNE
mappedX = tsne(data_after_ot(:,2:d+1));
% mappedX = importdata('/Users/coco/Workspace/code/OT-integration/data/pbmc/tsne/pbmc_seurat3.txt');
% mappedX = mappedX.data;

load carsmall;
clr1 = lines(dataset_num);
clr2 = lines(m0);

figure;
% subplot(3,3,1)
% h = gscatter(raw_data_tsne(1:index_list(2),1), raw_data_tsne(1:index_list(2),2), dataset_label);title(strcat('Unintegrated datasets'),'FontSize',20);
% set(h,'MarkerSize',12);
subplot(3,3,2)
h = gscatter(mappedX(:,1), mappedX(:,2), dataset_label, clr1);title(strcat('Batch-corrected datasets'),'FontSize',20);
set(h,'MarkerSize',12);
subplot(3,3,3)
h = gscatter(mappedX(:,1), mappedX(:,2), data_after_ot(:,1), clr2);title(strcat('Integrated datasets'),'FontSize',20);
set(h,'MarkerSize',12);

subplot(3,3,4)
h = gscatter(raw_data_tsne(1:index_list(2),1), raw_data_tsne(1:index_list(2),2), dataset_label, clr1);title(strcat('Unintegrated datasets'),'FontSize',20);
set(h,'MarkerSize',12);
subplot(3,3,5)
h = gscatter(mappedX(:,1), mappedX(:,2), dataset_label, clr1);title(strcat('Batch-corrected datasets'),'FontSize',20);
set(h,'MarkerSize',12);
subplot(3,3,6)
h = gscatter(mappedX(:,1), mappedX(:,2), seurat_clusters, clr2);title(strcat('Integrated datasets'),'FontSize',20);
set(h,'MarkerSize',12);

% subplot(3,3,7)
% h = gscatter(raw_data_tsne(1:index_list(2),1), raw_data_tsne(1:index_list(2),2), dataset_label);title(strcat('Unintegrated datasets'),'FontSize',20);
% set(h,'MarkerSize',12);
subplot(3,3,8)
h = gscatter(liger_tsne(:,1), liger_tsne(:,2), dataset_label, clr1);title(strcat('Batch-corrected datasets'),'FontSize',20);
set(h,'MarkerSize',12);
subplot(3,3,9)
h = gscatter(liger_tsne(:,1), liger_tsne(:,2), liger_clusters, clr2);title(strcat('Integrated datasets'),'FontSize',20);
set(h,'MarkerSize',12);