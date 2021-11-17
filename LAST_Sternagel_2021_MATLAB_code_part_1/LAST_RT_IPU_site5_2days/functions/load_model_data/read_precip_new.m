function [prec_int, prec_time] = read_precip_new(precip_file)

%% Initialize variables.
delimiter = ' ';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Format for each line of text:

formatSpec = '%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(precip_file,'r');

%% Read columns of data according to the format.

dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Create output variable
prec = [dataArray{1:end-1}];
prec_int = prec(:,2);
prec_time = prec(:,1);

%% Convert timeseries and precipitation intensity into SI units

conv_prec_int = (1/1000)/3600; % conversion factor for mm/h -> m³/m²*s
prec_int = prec_int * conv_prec_int;

conv_prec_time = 60; % conversion factor for min -> s
prec_time = prec_time * conv_prec_time;
