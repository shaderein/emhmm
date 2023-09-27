function [cluster, cluster_center,energy_in_iteration] = my_weighted_kmeans(cluster_num, it_max, ...
  point, weight, cluster_center)
% (internal function)
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2021-07-29
% Antoni B. Chan, Janet H. Hsiao, Lan Hui
% City University of Hong Kong, University of Hong Kong

% 2021-07-29: v0.80 - initial version

it_num = 0;

[dim_num, point_num] = size(point);

dist2cluster = zeros(cluster_num,point_num);
for j = 1:cluster_num
  dist2cluster(j,:) = sum(bsxfun(@minus, point, cluster_center(:,j)).^2,1); 
end
[dummy,cluster] = min(dist2cluster,[],1);

%  Determine the cluster populations and weights.
%  compute the weighted cluster center
[cluster_center,cluster_population,cluster_weight] = gcentroids(point, cluster, weight, cluster_num);

%  Set the point energies.
%   - 1. compute the distance of each point to its assigned cluster center
%   - 2. Adjust the point energies by a weight factor.
[f,cluster_energy] = genergy(point,weight,cluster_center,cluster_weight,cluster);

it_num = 0;
old_energy = sum(cluster_energy);
energy_in_iteration(1) = old_energy;

while ( it_num < it_max )
    % for each cluster
    %   find non members ~m
    %   compute f(~m) to current cluster -> f'(~m)
    %   save these energy f in a k x n matrix, with each element storing the energy of each point to corresponding cluster
    %   find minimum energy for each point
    %   assign all point to each minimum energy cluster
    %   update cluster center, weight, population and energy

    fmat = zeros(cluster_num,point_num);
    for j = 1:cluster_num
        members = (cluster==j);
        fmat(j,members) = f(members);
        nonmembers = ~members;
        adjust_weight = cluster_weight(j)./(cluster_weight(j)+weight(nonmembers));
        f_nonmembers = sum(bsxfun(@minus,point(:,nonmembers),cluster_center(:,j)).^2,1).*adjust_weight;
        fmat(j,nonmembers) = f_nonmembers;
    end

    [dummy,cluster] = min(fmat,[],1);

    [cluster_center,cluster_population,cluster_weight] = gcentroids(point, cluster, weight, cluster_num);
    [f,cluster_energy] = genergy(point,weight,cluster_center,cluster_weight,cluster);
    new_energy = sum(cluster_energy);
    if abs(new_energy-old_energy) < 1e-6
        break;
    else
        old_energy = new_energy;
        it_num = it_num + 1;
        energy_in_iteration(end+1) = new_energy;
    end 

end





function [centroids,cluster_population,cluster_weight] = gcentroids(point, cluster, weight, cluster_num)

[ndim,~] = size(point);
cluster_population = zeros(cluster_num,1);
cluster_weight = zeros(cluster_num,1);
centroids = zeros(ndim,cluster_num);

for j = 1:cluster_num
  members = (cluster == j);
  cluster_population(j) = sum(members);
  cluster_weight(j)     = sum(weight(members));

  centroids(:,j)        = sum(bsxfun(@times,point(:,members),weight(members)),2);
  if cluster_weight(j) > 0
    centroids(:,j) = centroids(:,j)/cluster_weight(j); % new cluster center
  end
end





function [f,cluster_energy] = genergy(point,weight,cluster_center,cluster_weight,cluster)

[ndim,ndata] = size(point);
cluster_num  = size(cluster_center,2);

f(1:ndata)   = 0.0;
cluster_energy(1:cluster_num) = 0.0;

for j = 1:cluster_num
  members = (cluster == j);
  f(members) = sum(bsxfun(@minus, point(:,members), cluster_center(:,j)).^2,1);
  cluster_energy(j) = sum(weight(members).*f(members));
  f(members) = f(members)*cluster_weight(j)./(bsxfun(@minus,cluster_weight(j),weight(members)));
end
