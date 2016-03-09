function [Z0] = load_truth_Immuno(dataset, k)
% One of those two delimiters should work
delimiterIn = '\t';
% delimiterIn = ' ';

% TODO delete some stuff to match the other loader.

## headerlinesIn = 0;

raw = load(dataset);
gmints = raw(:,1);
pos = raw(:,2);
Z0 = raw(:,3:5);

if nargin > 1
  % don't know where the gene numbers are starting, so shift them to 1.
  mask = (gmints - min(gmints) + 1) < no_genes;
  Z0 = Z0(mask);
  pos = pos(mask);
end

## raw = importdata(dataset,delimiterIn, headerlinesIn);
## Z0 = raw.data(:, 2:4);
## pos = raw.data(:, 1) + 5000;
## [gmnames, gmints]= unique(raw.textdata,'stable');
## if nargin > 1
## 	% only take first k genes
## 	Z0 = Z0(1:gmints(k+1)-1,:);
## 	pos = pos(1:gmints(k+1)-1);
## 	gmnames = gmnames(1:k);
## 	gmints = gmints(1:k);
## end

## % Jerry messed up, reverse the genes now...
## % We will not need to do this for future datasets
## gmints(length(gmints)+1) = length(pos)+1;
## for i = 1:length(gmints)-1
##     Z0(gmints(i):gmints(i+1)-1,:) = flipud(Z0(gmints(i):gmints(i+1)-1,:));
##     pos(gmints(i):gmints(i+1)-1) = flipud(pos(gmints(i):gmints(i+1)-1));
## end
## gmints(end)=[];
## end
