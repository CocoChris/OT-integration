clear;clc;
addpath('../barycenter/');

dataset_num=4;
cluster_num_list=[14,14,14,14];
d=13;

no_dims = 2;
initial_dims = 10;
perplexity = 30;

supp=zeros(d+d*d,sum(cluster_num_list));

prefix='~/Workspace/code/OT-integration/data/human_pancreas/';

k=1;
cluster_weights=zeros(1,sum(cluster_num_list));
data_before_ot = [];
index_list=zeros(1,dataset_num);

samples = char('human1', 'human2', 'human3', 'human4');

for i=1:dataset_num
    sample_name = samples(i,:);
    cluster_label_filename=strcat(prefix, sample_name, '_labels.txt');
    cca_filename=strcat(prefix, 'cca/', sample_name, '.txt');
    
    cluster_labels=importdata(cluster_label_filename);
%     cluster_labels=raw_data(:,1);
    cca_data=importdata(cca_filename);
    cca_data=cca_data.data(:,1:d);
    
    all_data=[cluster_labels cca_data];
    sample_num=size(cca_data,1);
    
    cca_data_group=cell(1,cluster_num_list(i));

    for j=1:cluster_num_list(i)
      cca_data_group{j}=all_data(all_data(:,1)==j,2:end);
    end
    
    for j=1:cluster_num_list(i)
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
    stride(i)=cluster_num_list(i);
end
% instanceW=rand(1,dataset_num); % not necessarily sum2one
instanceW = [1,1,1,1];
% compute GMM barycenter with B-ADMM method allowing the change of component weights
% initiate centroid from an instance
m0=14;
c0.supp=supp(:,1:m0);
c0.w=cluster_weights(1:m0);
% set the number of iterations and badmm_rho (no need to change)
options.badmm_max_iters=2000;
options.badmm_rho=10;
% compute the centroid
[c, X]=centroid_sphBregman_GMM_test(stride, instanceW, supp, cluster_weights, c0, options);

index_list = cumsum(index_list);
data_after_ot = re_cluster(c.supp, data_before_ot(:,2:d+1));

dataset_label = [];
dataset_label(1:index_list(1)) = 1;
dataset_label((index_list(1)+1):index_list(2)) = 2;
dataset_label((index_list(2)+1):index_list(3)) = 3;
dataset_label((index_list(3)+1):index_list(4)) = 4;

ot_clusters = data_after_ot(:,1);

seurat_clusters = importdata('~/Workspace/code/OT-integration/data/human_pancreas/clusters/seurat_human_clusters.txt');
seurat_clusters = seurat_clusters.data(:,2);

liger_clusters = importdata('~/Workspace/code/OT-integration/data/human_pancreas/clusters/liger_human_clusters.txt');
liger_clusters = liger_clusters.data(:,2);

liger_tsne = importdata('~/Workspace/code/OT-integration/data/human_pancreas/tsne/liger_human.txt');
liger_tsne = liger_tsne.data;

raw_data_human1 = importdata('~/Workspace/code/OT-integration/data/human_pancreas/pca/human1.txt');
raw_data_human1 = raw_data_human1.data;
raw_data_human2 = importdata('~/Workspace/code/OT-integration/data/human_pancreas/pca/human2.txt');
raw_data_human2 = raw_data_human2.data;
raw_data_human3 = importdata('~/Workspace/code/OT-integration/data/human_pancreas/pca/human3.txt');
raw_data_human3 = raw_data_human3.data;
raw_data_human4 = importdata('~/Workspace/code/OT-integration/data/human_pancreas/pca/human4.txt');
raw_data_human4 = raw_data_human4.data;

raw_data_union = [raw_data_human1; raw_data_human2; raw_data_human3; raw_data_human4];

raw_data_tsne = tsne(raw_data_union);

load carsmall;
clr1 = lines(4);

color_candy_hex = char('#EF0000', '#336699', '#FEC211', '#3BC371', '#666699', '#999999', '#FF6666', '#6699CC', '#CC6600', '#009999', '#6B67BC', '#99867A', '#CC3333', '#669999', '#CC9900');
% color_candy_hex = char('#EF0000', '#336699', '#FEC211', '#3BC371', '#CC3333', '#999999', '#669999', '#6699CC', '#CC6600', '#009999', '#6B67BC', '#99867A', '#666699', '#FF6666', '#CC9900');
color_candy = zeros(15, 3);
for i = 1:15
    hex = color_candy_hex(i,:);
    color_candy(i,:) = hex2rgb(hex);
end
    
clr2 = color_candy(1:m0,:);


data_before_tsne = [data_before_ot;data_after_ot];
% Run t-SNE
mappedX = tsne(data_after_ot(:,2:d+1));
% mappedX = importdata('~/Workspace/code/OT-integration/data/human_pancreas/tsne/seurat_human.txt');
% mappedX = mappedX.data;

figure;
% subplot(3,3,1)
% h = gscatter(raw_data_tsne(:,1), raw_data_tsne(:,2), dataset_label, clr1);title(strcat('Unintegrated datasets'), 'FontSize',20);
% set(h,'MarkerSize',12);
subplot(3,3,2)
h = gscatter(mappedX(:,1), mappedX(:,2), dataset_label, clr1);title(strcat('Batch-corrected datasets'), 'FontSize',20);
set(h,'MarkerSize',5);
subplot(3,3,3)
h = gscatter(mappedX(:,1), mappedX(:,2), ot_clusters, clr2);title(strcat('Integrated datasets'), 'FontSize',20);
set(h,'MarkerSize',12);

subplot(3,3,4)
h = gscatter(raw_data_tsne(:,1), raw_data_tsne(:,2), dataset_label, clr1);title(strcat('Unintegrated datasets'), 'FontSize',20);
set(h,'MarkerSize',12);
subplot(3,3,5)
h = gscatter(mappedX(:,1), mappedX(:,2), dataset_label, clr1);title(strcat('Batch-corrected datasets'), 'FontSize',20);
set(h,'MarkerSize',5);
subplot(3,3,6)
h = gscatter(mappedX(:,1), mappedX(:,2), seurat_clusters, clr2);title(strcat('Integrated datasets'), 'FontSize',20);
set(h,'MarkerSize',12);

% subplot(3,3,7)
% h = gscatter(raw_data_tsne(:,1), raw_data_tsne(:,2), dataset_label, clr1);title(strcat('Unintegrated datasets'), 'FontSize',20);
% set(h,'MarkerSize',12);
subplot(3,3,8)
h = gscatter(liger_tsne(:,1), liger_tsne(:,2), dataset_label, clr1);title(strcat('Batch-corrected datasets'), 'FontSize',20);
set(h,'MarkerSize',5);
subplot(3,3,9)
h = gscatter(liger_tsne(:,1), liger_tsne(:,2), liger_clusters, clr2);title(strcat('Integrated datasets'), 'FontSize',20);
set(h,'MarkerSize',12);

fid=fopen(strcat(prefix, 'gamma_human.txt'),'w');
data_output=real(X);
[r,c]=size(data_output);            
 for i=1:r
  for j=1:c
    fprintf(fid,'%f\t',data_output(i,j));
  end
    fprintf(fid,'\n');
 end
fclose(fid);