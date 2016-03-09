function [Z0] = load_truth_Immuno(dataset, k)
% One of those two delimiters should work
delimiterIn = '\t';
% delimiterIn = ' ';

% TODO delete some stuff to match the other loader.

headerlinesIn = 0;
datapath = './data/';

if exist(datapath) ~= 7
  error('cannot find ./data symlinked directory');
end

raw = importdata([datapath dataset],delimiterIn, headerlinesIn);
Z0 = raw.data(:, 2:4);
pos = raw.data(:, 1) + 5000;
[gmnames, gmints]= unique(raw.textdata,'stable');
if nargin > 1
	% only take first k genes
	Z0 = Z0(1:gmints(k+1)-1,:);
	pos = pos(1:gmints(k+1)-1);
	gmnames = gmnames(1:k);
	gmints = gmints(1:k);
end

% Jerry messed up, reverse the genes now...
% We will not need to do this for future datasets
gmints(length(gmints)+1) = length(pos)+1;
for i = 1:length(gmints)-1
    Z0(gmints(i):gmints(i+1)-1,:) = flipud(Z0(gmints(i):gmints(i+1)-1,:));
    pos(gmints(i):gmints(i+1)-1) = flipud(pos(gmints(i):gmints(i+1)-1));
end
gmints(end)=[];
end
