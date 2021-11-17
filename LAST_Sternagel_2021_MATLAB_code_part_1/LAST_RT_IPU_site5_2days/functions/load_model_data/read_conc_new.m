function conc = read_conc_new(conc_file)

%% Initialize variables.
delimiter = ' ';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Format for each line of text:

formatSpec = '%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(conc_file,'r');

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
conc = [dataArray{1:end-1}];
conc = conc(:,2);

%% Converts tracer concentration in precipitation water into SI units

conv_conc = 1; % conversion factor for kg/m³ -> kg/m³
conc = conc * conv_conc;


