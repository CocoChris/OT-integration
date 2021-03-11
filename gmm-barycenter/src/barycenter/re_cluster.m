function [D, dist_mean] = re_cluster(supp, data)
% size(supp) = [d*d+d, m]
% size(data) = [n, d]

d = size(data, 2);
m = size(supp, 2);
mean = supp(1:d,:);
n = size(data, 1);
D = zeros(n, d+1);
D(:,2:d+1) = data;

dist_matrix = zeros(n, m);
dist_mean = 0;

for i=1:n
    dist_list = zeros(1, m);
    for j=1:m
        dist_list(j) = sum((mean(:,j) - data(i,:)').^2);
    end
    dist_matrix(i,:) = dist_list;
%     if i < 4
%         disp(dist_list)
%     end
    [min_dist, cluster_label] = min(dist_list);
    dist_sorted = sort(dist_list);
    d1 = dist_sorted(1);
    d2 = dist_sorted(2);
    dist_mean = dist_mean + sqrt(d1/d2);
    D(i,1) = cluster_label;
end

dist_mean = dist_mean / n;