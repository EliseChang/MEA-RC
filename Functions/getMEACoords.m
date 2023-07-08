function [channel_IDs, channel_idx, coords] = getMEACoords

channel_IDs = [11, 12, 13, 14, 15, 16, 17, 18, ... 
            21, 22, 23, 24, 25, 26, 27, 28, ...
            31, 32, 33, 34, 35, 36, 37, 38, ...
            41, 42, 43, 44, 45, 46, 47, 48, ...
            51, 52, 53, 54, 55, 56, 57, 58, ...
            61, 62, 63, 64, 65, 66, 67, 68, ...,
            71, 72, 73, 74, 75, 76, 77, 78, ...,
            81, 82, 83, 84, 85, 86, 87, 88];

channelsOrdering = ...
[21 31 41 51 61 71 12 22 32 42 52 62 72 82 13 23 33 43 53 63 ... 
73 83 14 24 34 44 54 64 74 84 15 25 35 45 55 65 75 85 16 26 ...
36 46 56 66 76 86 17 27 37 47 57 67 77 87 28 38 48 58 68 78];

channel_idx = ...
[nan 24 26 29 32 35 37 nan 21 22 25 30 31 36 39 40 19 20 23 28 33 38 41 42 ...
16 17 18 27 34 43 44 45 15 14 13 4 57 48 47 46 12 11 8 3 58 53 50 49 ...
10 9 6 1 60 55 52 51 nan 7 5 2 59 56 54 nan];

coords = zeros(length(channel_IDs), 2);
coords(:, 2) = repmat(linspace(1, 0, 8), 1, 8);
coords(:, 1) = repelem(linspace(0, 1, 8), 1, 8);

subset_idx = find(~ismember(channel_IDs, [11, 81, 18, 88]));
channel_IDs = channel_IDs(subset_idx);

reorderingIdx = zeros(length(channel_IDs), 1);
for n = 1:length(channel_IDs)
    reorderingIdx(n) = find(channelsOrdering(n) == channel_IDs);
end 

coords = coords(subset_idx, :);

% Re-order the channel IDs and coordinates to match the original
% ordering
channel_IDs = channel_IDs(reorderingIdx);
coords = coords(reorderingIdx, :);
