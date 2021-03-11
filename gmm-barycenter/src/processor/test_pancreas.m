clear;clc;
addpath('../barycenter/');

dataset_num=2;
cluster_num_list=[14,13];
d=15;

no_dims = 2;
initial_dims = 10;
perplexity = 30;

supp=zeros(d+d*d,sum(cluster_num_list));

prefix='/Users/coco/Workspace/Git/OT-integration/data/species/';

k=1;
cluster_weights=zeros(1,sum(cluster_num_list));
data_before_ot = [];
index_list=zeros(1,dataset_num);

samples = char('human', 'mouse');

for i=1:dataset_num
    sample_name = samples(i,:);
    cluster_label_filename=strcat(prefix, 'clusters/', sample_name, '_labels.txt');
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
    
    % Run t-SNE
%     mappedX = tsne(all_count_data, 'Algorithm','exact', 'NumDimensions', no_dims, 'NumPCAComponents', initial_dims, 'Perplexity', perplexity);
    % Plot results
%     subplot(2,2,i)
%     gscatter(mappedX(:,1), mappedX(:,2), cluster_labels);title(strcat('human', num2str(i), ' data'));
end

stride = zeros(1,dataset_num);
for i=1:dataset_num
    stride(i)=cluster_num_list(i);
end
% instanceW=rand(1,dataset_num); % not necessarily sum2one
instanceW = [1,1];
% compute GMM barycenter with B-ADMM method allowing the change of component weights
% initiate centroid from an instance
m0=13;
c0.supp=supp(:,1:m0);
c0.w=cluster_weights(1:m0);
% set the number of iterations and badmm_rho (no need to change)
options.badmm_max_iters=2000;
options.badmm_rho=10;
% compute the centroid
[c, X]=centroid_sphBregman_GMM_test(stride, instanceW, supp, cluster_weights, c0, options);

index_list = cumsum(index_list);
data_after_ot = re_cluster(c.supp, data_before_ot(:,2:d+1));

% cluster_distribution = zeros(14, 14);
% for i=1:index_list(4)
%     old_label = data_before_ot(i,1);
%     new_label = data_after_ot(i,1);
%     cluster_distribution(old_label, new_label) = cluster_distribution(old_label, new_label) + 1;
% end

% disp('total samples: ')
% disp(index_list(end))
% disp('correct samples: ')
% disp(sum(diag(cluster_distribution)))
% disp('correct rate: ')
% disp(sum(diag(cluster_distribution)) / index_list(end))

dataset_label = [];
dataset_label(1:1937) = 1;
dataset_label(1938:3661) = 2;
dataset_label(3662:7266) = 3;
dataset_label(7267:8569) = 4;
dataset_label(8570:9391) = 5;
dataset_label(9392:10455) = 6;

ot_clusters = data_after_ot(:,1);

seurat_clusters = importdata('/Users/coco/Workspace/Git/OT-integration/data/species/clusters/seurat_clusters.txt');
seurat_clusters = seurat_clusters.data(:,2);

ex2 = zeros(9,2);
ex2(:,1) = 1:9;
ex2(:,2) = [5,8,7,3,1,4,2,9,6];

seurat_clusters_final = zeros(1,index_list(2));
for i = 1:index_list(2)
    seurat = seurat_clusters(i);
    if (ismember(seurat, ex2(:,1)) == 1)
        k = find(ex2(:,1) == seurat);
        seurat = ex2(k,2);
    end
    seurat_clusters_final(i) = seurat;
end 

liger_clusters = importdata('/Users/coco/Workspace/Git/OT-integration/data/species/clusters/liger_clusters.txt');
liger_clusters = liger_clusters.data(:,2);

liger_tsne = importdata('/Users/coco/Workspace/Git/OT-integration/data/species/liger/pancreas.txt');
liger_tsne = liger_tsne.data;


raw_data_human = importdata('/Users/coco/Workspace/Git/OT-integration/data/species/pca/human.txt');
raw_data_human = raw_data_human.data;
raw_data_mouse = importdata('/Users/coco/Workspace/Git/OT-integration/data/species/pca/mouse.txt');
raw_data_mouse = raw_data_mouse.data;

raw_data_union = [raw_data_human; raw_data_mouse];

raw_data_tsne = tsne(raw_data_union);


data_before_tsne = [data_before_ot;data_after_ot];
% Run t-SNE
% mappedX = tsne(data_after_ot(:,2:d+1));
mappedX = importdata('/Users/coco/Workspace/Git/OT-integration/data/species/tsne/pancreas.txt');
mappedX = mappedX.data;
% Plot results
% figure;
% subplot(2,3,1)
% gscatter(mappedX(1:index_list(4),1), mappedX(1:index_list(4),2), sample_index);title(strcat('grouped by dataset'));
% subplot(2,3,2)
% gscatter(mappedX(1:index_list(4),1), mappedX(1:index_list(4),2), seurat_clusters);title(strcat('Integrated by Seurat'));
% subplot(2,3,3)
% gscatter(mappedX(1:index_list(4),1), mappedX(1:index_list(4),2), data_after_ot(:,1));title(strcat('Integrated by OT'));

load carsmall;
clr1 = lines(6);

color_candy_hex = char('#EF0000', '#336699', '#FEC211', '#3BC371', '#666699', '#999999', '#FF6666', '#6699CC', '#CC6600', '#009999', '#6B67BC', '#99867A', '#CC3333', '#669999', '#CC9900');
color_candy = zeros(15, 3);
for i = 1:15
    hex = color_candy_hex(i,:);
    color_candy(i,:) = hex2rgb(hex);
end
    
clr2 = color_candy(1:m0,:);

figure;
% subplot(3,3,1)
% h = gscatter(raw_data_tsne(:,1), raw_data_tsne(:,2), dataset_label);title(strcat('Unintegrated datasets'), 'FontSize',20);
% set(h,'MarkerSize',12);
subplot(3,3,2)
h = gscatter(mappedX(:,1), mappedX(:,2), dataset_label, clr1);title(strcat('Batch-corrected datasets'), 'FontSize',20);
set(h,'MarkerSize',12);
subplot(3,3,3)
h = gscatter(mappedX(:,1), mappedX(:,2), ot_clusters, clr2);title(strcat('Integrated datasets'), 'FontSize',20);
set(h,'MarkerSize',12);

subplot(3,3,4)
h = gscatter(raw_data_tsne(:,1), raw_data_tsne(:,2), dataset_label, clr1);title(strcat('Unintegrated datasets'), 'FontSize',20);
set(h,'MarkerSize',12);
subplot(3,3,5)
h = gscatter(mappedX(:,1), mappedX(:,2), dataset_label, clr1);title(strcat('Batch-corrected datasets'), 'FontSize',20);
set(h,'MarkerSize',12);
subplot(3,3,6)
h = gscatter(mappedX(:,1), mappedX(:,2), seurat_clusters, clr2);title(strcat('Integrated datasets'), 'FontSize',20);
set(h,'MarkerSize',12);

% subplot(3,3,7)
% gscatter(raw_data_tsne(:,1), raw_data_tsne(:,2), dataset_label);title(strcat('Unintegrated datasets'));
% set(h,'MarkerSize',12);
subplot(3,3,8)
h = gscatter(liger_tsne(:,1), liger_tsne(:,2), dataset_label, clr1);title(strcat('Batch-corrected datasets'), 'FontSize',20);
set(h,'MarkerSize',12);
subplot(3,3,9)
h = gscatter(liger_tsne(:,1), liger_tsne(:,2), liger_clusters, clr2);title(strcat('Integrated datasets'), 'FontSize',20);
set(h,'MarkerSize',12);

fid=fopen('/Users/coco/Workspace/Git/OT-integration/data/species/clusters/ot_clusters.txt','w');%写入文件路径
for i=1:length(data_after_ot(:,1))
fprintf(fid,'%d\n',data_after_ot(i,1));   
end
fclose(fid);

fid=fopen(strcat(prefix, 'ot/X.txt'),'w');
data_output=real(X);
[r,c]=size(data_output);            
 for i=1:r
  for j=1:c
    fprintf(fid,'%f\t',data_output(i,j));
  end
    fprintf(fid,'\n');
 end
fclose(fid);

origin_pancreas = data_before_ot(1:index_list(2),1);
ot_pancreas = data_after_ot(1:index_list(2),1);
seurat_pancreas = seurat_clusters;
liger_pancreas = liger_clusters;
[AR1,RI1,MI1,HI1]=RandIndex(origin_pancreas,ot_pancreas);
[AR2,RI2,MI2,HI2]=RandIndex(origin_pancreas,seurat_pancreas);
[AR3,RI3,MI3,HI3]=RandIndex(origin_pancreas,liger_pancreas);

