function [indices] = getElectrodeIdx(IDs, options)
% Gets indices of electrodes (used to index spikeTimes, etc.), given
% electrode spatial co-ordinates

arguments
    IDs
    options.method = ""; % "Python" for zero indexing
    options.ref_node = true; % "false" if reference electrode has already been removed
end


try
    channels = load("channels.mat").channels;
catch
    channels = [47,	48,	46,	45,	38,	37,	28,	36,	27,	17,	26,	16,	35,	25,	15,...
        14,	24,	34,	13,	23,	12,	22,	33,	21,	32,	31,	44,	43,	41,	42,	52,	51,...
        53,	54,	61,	62,	71,	63,	72,	82,	73,	83,	64,	74,	84,	85,	75,	65,	86,...
        76,	87,	77,	66,	78,	67,	68,	55,	56,	58,	57];
end

indices = zeros(length(IDs),1);

for e = 1:length(IDs)
    indices(e) = find(channels == IDs(e));
end

if strcmp(options.method,"Python")
    indices = indices - 1;
end

if ~options.ref_node
    indices(indices >= 15) = indices(indices >= 15) - 1;
end

