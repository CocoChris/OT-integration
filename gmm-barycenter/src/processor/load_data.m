clear;clc;
profix = '/Users/coco/Files/PhD/Intern/program/OT/data/pancreatic_islet_scRNA_seq/GSE84133_RAW/';
filename = strcat(profix, 'GSM2230757_human1_umifm_counts.txt');

raw_data = importdata(filename);
cluster_labels = raw_data(:,1);
count_data = raw_data(:,2:end);

no_dims = 2;
initial_dims = 50;
perplexity = 30;

% Run t?SNE
% mappedX = tsne(count_data, 'Algorithm','exact', 'NumDimensions', no_dims, 'NumPCAComponents', initial_dims, 'Perplexity', perplexity);
% Plot results
% gscatter(mappedX(:,1), mappedX(:,2), cluster_labels);

% compute mean and variance
% sample_mean = mean(count_data, 1);
% sample_cov = cov(count_data);

cluster_label_uni = unique(cluster_labels);
n = length(cluster_label_uni);
count_data_group = cell(1,n);

for i = 1:n
  count_data_group{i} = raw_data(raw_data(:,1)==i-1,2:end);
end



