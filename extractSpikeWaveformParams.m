% Directories
HomeDir = 'D:\MATLAB\MPhil_scripts\MEA-RC'; % analysis folder to home directory
spikeDetectedData = 'D:\MATLAB\MPhil_scripts\MEA-RC\organoids'; % path to output folder containing 1_SpikeDetection folder
cd(HomeDir)
addpath(spikeDetectedData)
addpath Functions\

% Input and output filetype
spreadsheet_filename = 'readable.csv';
spreadsheet_file_type = 'csv';
csvRange = [2, Inf]; % read the data in the range [StartRow EndRow], e.g. [2 Inf] means start reading data from row 2
setUpSpreadSheet

method = 'bior1p5';
spikeWaveformParams = zeros(length(ExpName),3);
% store minAmp in first col., minSlope in second col. and maxSlope in third
% col.
for Exp = 1:length(ExpName)
    cd(spikeDetectedData)
    load(strcat(ExpName{Exp}, '_spikes.mat'), "spikeWaveforms")
    cd(HomeDir)

    [minAmp, minSlope, maxSlope] = getSpikeWaveform(spikeWaveforms, method);
    spikeWaveformParams(Exp,1) = minAmp;
    spikeWaveformParams(Exp,2) = minSlope;
    spikeWaveformParams(Exp,3) = maxSlope;
    
    clear minAmp minSlope maxSlope
end

mean(spikeWaveformParams)