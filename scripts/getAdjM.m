homeDir = ("D:\MATLAB\MEA-RC");
cd(homeDir)
addpath Functions
spreadsheet_filename = "mecp2timepoints.csv";
spreadsheet_file_type = 'csv';
csvRange = [2, Inf];
Params.channelLayout = 'MCS60';
setUpSpreadSheet
srcDir = "D:\MATLAB\MEA-NAP\outputs\OutputData27Oct2023\ExperimentMatFiles";
srcDate = "27Oct2023";
outputDir = "C:\Users\elise\Python\ReservoirComputing\data\connectivity";

lags = [10, 25, 50];

cd(srcDir)

for l = lags

    for ExN = 1:length(ExpName)
        load(strcat(ExpName{ExN}, "_", srcDate), "-mat", "adjMs")
        adjM = adjMs.(strcat("adjM", num2str(l), "mslag"));
        mat2np(adjM, fullfile(outputDir, strcat(ExpName{ExN}, "_", num2str(l), "mslag")),...
            'float64')
    end

end

cd(homeDir)