clear;clc;
addpath('../barycenter/');

dataset_num=3;
cluster_num=7;
d=15;

no_dims = 2;
initial_dims = 10;
perplexity = 30;

supp=zeros(d+d*d,dataset_num*cluster_num);

prefix='/Users/coco/Files/PhD/Intern/program/OT/data/pancreatic_islet_scRNA_seq/GSE84133_RAW/';

k=1;
cluster_weights=zeros(1,dataset_num*cluster_num);
data_before_ot = [];
index_list=zeros(1,dataset_num);

for i=1:dataset_num
    cluster_label_filename=strcat(prefix, 'rawdata_cluster_filter/cluster_label/', 'human', num2str(i), '.txt');
    cca_filename=strcat(prefix, 'rawdata_cluster_filter/cca/', 'human', num2str(i), '.txt');
    
    cluster_labels=importdata(cluster_label_filename);
%     cluster_labels=raw_data(:,1);
    cca_data=importdata(cca_filename);
    cca_data=cca_data.data(:,1:d);
    
    all_data=[cluster_labels cca_data];
    sample_num=size(cca_data,1);
    
    cca_data_group=cell(1,cluster_num);

    for j=1:cluster_num
      cca_data_group{j}=all_data(all_data(:,1)==j,2:end);
    end
    
    for j=1:cluster_num
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
    stride(i)=cluster_num;
end
% instanceW=rand(1,dataset_num); % not necessarily sum2one
instanceW = [1,1,1];
% compute GMM barycenter with B-ADMM method allowing the change of component weights
% initiate centroid from an instance
m0=cluster_num;
c0.supp=supp(:,1:m0);
c0.w=cluster_weights(1:m0);
% set the number of iterations and badmm_rho (no need to change)
options.badmm_max_iters=2000;
options.badmm_rho=10;
% compute the centroid
[c, X]=centroid_sphBregman_GMM_test(stride, instanceW, supp, cluster_weights, c0, options);

index_list = cumsum(index_list);
data_after_ot = re_cluster(c.supp, data_before_ot(:,2:d+1));

cluster_distribution = zeros(cluster_num, cluster_num);
total_sample_num = size(data_before_ot,1);
for i=1:total_sample_num
    old_label = data_before_ot(i,1);
    new_label = data_after_ot(i,1);
    cluster_distribution(old_label, new_label) = cluster_distribution(old_label, new_label) + 1;
end

disp('total samples: ')
disp(index_list(end))
disp('correct samples: ')
disp(sum(diag(cluster_distribution)))
disp('correct rate: ')
disp(sum(diag(cluster_distribution)) / index_list(end))

data_before_tsne = [data_before_ot;data_after_ot];
% Run t-SNE
mappedX = tsne(data_before_tsne(:,2:d+1));
% Plot results
figure;
subplot(2,3,1)
gscatter(mappedX(1:index_list(1),1), mappedX(1:index_list(1),2), data_before_ot(1:index_list(1),1));title(strcat('human', num2str(1), ' data'));
subplot(2,3,2)
gscatter(mappedX((index_list(1)+1):index_list(2),1), mappedX((index_list(1)+1):index_list(2),2), data_before_ot((index_list(1)+1):index_list(2),1));title(strcat('human', num2str(2), ' data'));
subplot(2,3,3)
gscatter(mappedX((index_list(2)+1):index_list(3),1), mappedX((index_list(2)+1):index_list(3),2), data_before_ot((index_list(2)+1):index_list(3),1));title(strcat('human', num2str(3), ' data'));
subplot(2,3,4)
gscatter(mappedX(1:index_list(3),1), mappedX(1:index_list(3),2), data_before_ot(:,1));title(strcat('overlap'));
subplot(2,3,5)
gscatter(mappedX((index_list(3)+1):end,1), mappedX((index_list(3)+1):end,2), data_after_ot(:,1));title(strcat('after integration'));

