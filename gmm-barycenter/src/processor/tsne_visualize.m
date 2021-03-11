% Load data
load 'mnist_train.mat'
ind = randperm(size(train_X, 1));
train_X = train_X(ind(1:500),:);
train_labels = train_labels(ind(1:500));
% Set parameters
no_dims = 2;
initial_dims = 50;
perplexity = 30;
% Run t?SNE
mappedX = tsne(train_X, 'Algorithm','exact', 'NumDimensions', no_dims, 'NumPCAComponents', initial_dims, 'Perplexity', perplexity);
% Plot results
gscatter(mappedX(:,1), mappedX(:,2), train_labels);