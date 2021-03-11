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

species = ['human', 'mouse'];

for i=1:dataset_num
    if i == 1
        t = 'human';
    else
        t = 'mouse';
    end
    cluster_label_filename=strcat(prefix, 'clusters/', t, '_labels.txt');
    cca_filename=strcat(prefix, 'cca/', t, '.txt');
    
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
m0=10;
c0.supp=supp(:,1:m0);
c0.w=cluster_weights(1:m0);
% set the number of iterations and badmm_rho (no need to change)
options.badmm_max_iters=2000;
options.badmm_rho=10;
% compute the centroid
[c, X]=centroid_sphBregman_GMM_test(stride, instanceW, supp, cluster_weights, c0, options);

index_list = cumsum(index_list);
data_after_ot = re_cluster(c.supp, data_before_ot(:,2:d+1));

cluster_distribution1 = zeros(cluster_num_list(1), 13);
for i=1:index_list(1)
    old_label = data_before_ot(i,1);
    new_label = data_after_ot(i,1);
    cluster_distribution1(old_label, new_label) = cluster_distribution1(old_label, new_label) + 1;
end

cluster_distribution2 = zeros(cluster_num_list(2), 13);
for i=(index_list(1)+1):index_list(2)
    old_label = data_before_ot(i,1);
    new_label = data_after_ot(i,1);
    cluster_distribution2(old_label, new_label) = cluster_distribution2(old_label, new_label) + 1;
end

% disp('total samples: ')
% disp(index_list(end))
% disp('correct samples: ')
% disp(sum(diag(cluster_distribution)))
% disp('correct rate: ')
% disp(sum(diag(cluster_distribution)) / index_list(end))

sample_index = [];
sample_index(1:1937) = 1;
sample_index(1938:3661) = 2;
sample_index(3662:7266) = 3;
sample_index(7267:8569) = 4;
sample_index(8570:9391) = 5;
sample_index(9392:10455) = 6;

raw_data_human = importdata('/Users/coco/Workspace/Git/OT-integration/data/species/human.txt');
raw_data_human = raw_data_human.data;
raw_data_human = raw_data_human';
raw_data_human = raw_data_human(:,1:1000);
raw_data_mouse = importdata('/Users/coco/Workspace/Git/OT-integration/data/species/mouse.txt');
raw_data_mouse = raw_data_mouse.data;
raw_data_mouse = raw_data_mouse';
raw_data_mouse = raw_data_mouse(:,1:1000);

raw_data_union = [raw_data_human; raw_data_mouse];

raw_data_tsne = tsne(raw_data_union);


data_before_tsne = [data_before_ot;data_after_ot];
% Run t-SNE
mappedX = tsne(data_after_ot(:,2:d+1));
% Plot results
figure;
subplot(2,3,1)
gscatter(raw_data_tsne(1:index_list(2),1), raw_data_tsne(1:index_list(2),2), sample_index);title(strcat('Unintegrated datasets'));
subplot(2,3,2)
gscatter(mappedX(1:index_list(2),1), mappedX(1:index_list(2),2), sample_index);title(strcat('Integrated dataset'));
subplot(2,3,3)
gscatter(mappedX(1:index_list(2),1), mappedX(1:index_list(2),2), data_after_ot(:,1));title(strcat('Integrated dataset'));


% c1_human = data_before_ot(1:index_list(1),1);
% c2_human = data_after_ot(1:index_list(1),1);
% [AR,RI,MI,HI]=RandIndex(c1_human,c2_human);
% 
% c1_mouse = data_before_ot((index_list(1)+1):end,1);
% c2_mouse = data_after_ot((index_list(1)+1):end,1);
% [AR,RI,MI,HI]=RandIndex(c1_mouse,c2_mouse);
