function [downSampledMatrix] = downSampleMean(M,binLength)
%DOWNSAMPLEMEAN Summary of this function goes here
%   Detailed explanation goes here

nCols = size(M, 2);
M = reshape(M, binLength, nCols, []);
nBins = size(M, 3);
temp_output = mean(M, 1);
downSampledMatrix = reshape(temp_output, nBins, nCols);

end

