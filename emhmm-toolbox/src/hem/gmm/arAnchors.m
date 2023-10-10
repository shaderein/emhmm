
function [centers, mincenter, mindist, lower, computed, usedIdx, notUsedIdx] = arAnchors(firstcenter,k,data)
% choose k centers by the furthest-first method
%firstcenter and centers are indices into data, not points in a space.

[n,dim] = size(data);
centers = zeros(k,dim);
lower = zeros(n,k);
mindist = Inf*ones(n,1);
mincenter = ones(n,1);
computed = 0;
centdist = zeros(k,k);
usedIdx(1) = firstcenter;
for j = 1:k
    if j == 1
        newcenter = data(firstcenter,:);
    else
        [maxradius,i] = max(mindist);
        newcenter = data(i,:);
        usedIdx(j) = i;
    end

    centers(j,:) = newcenter;
    centdist(1:j-1,j) = calcdist(centers(1:j-1,:),newcenter);
    centdist(j,1:j-1) = centdist(1:j-1,j)';
    computed = computed + j-1;
    
    inplay = find(mindist > centdist(mincenter,j)/2);
    newdist = calcdist(data(inplay,:),newcenter);
    computed = computed + size(inplay,1);
    lower(inplay,j) = newdist;
        
%    other = find(mindist <= centdist(mincenter,j)/2);
%    if ~isempty(other)
%        lower(other,j) = centdist(mincenter(other),j) - mindist(other);
%    end    
        
    move = find(newdist < mindist(inplay));
    shift = inplay(move);
    mincenter(shift) = j;
    mindist(shift) = newdist(move);
end


notUsedIdx = find(mindist);

%udist = calcdist(data,centers(mincenter,:));
%quality = mean(udist);
%q2 = mean(udist.^2);
%[k toc quality q2 computed]
%toc