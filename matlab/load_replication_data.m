function load_replication_data(mat_file)
%LOAD_REPLICATION_DATA Load data and normalize date variables for replication.
%
% The frozen paper input contains Matlab datetime objects. The R-generated
% input stores dates as Matlab serial date numbers. This helper supports both.

if nargin < 1
    mat_file = '../data/data_JPR.mat';
end

loaded = load(mat_file);
names = fieldnames(loaded);
for i = 1:numel(names)
    assignin('caller', names{i}, loaded.(names{i}));
end

if isfield(loaded, 'tT') && isnumeric(loaded.tT)
    assignin('caller', 'tT', datetime(loaded.tT, 'ConvertFrom', 'datenum'));
end
end
