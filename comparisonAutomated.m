% comparison script which is used to test the validity of Katies
% SpikeExtractionTool

% University of Melbourne
% Department of Biomedical Engineering
% Summer Research Experience Program - Spike Simulation Suite
% Christopher Ong & Chi Yung Darren Tan | 08 March 2021

%% How to use:

% The comparison tool will automate testing the accuracy of the extraction
% tool

% Initially, user must generate extraction and simulation pairs with files
% name being:
% 'simulation#' and 'extraction#' where # is replaced by a integer
% Simulation files need to be in a simulation folder
% Extraction files need to be in a extraction folder
% Adjust the directory pathway accordingly to where the Simulation and
% Extraction folders are with the corresponding pair in the correct folder

% User then can run this tool and it will automatically go through all
% simulation/extraction pairs and produce and excel spreadsheet with all
% data recorded.

%% Main
close all
clear

% Load in all files for both simulation and extraction
fileIDSimulation = struct2cell(dir(['/Users/darrentan/Documents/SREP/Simulations/Data collection/sim14/Simulations' '/*.mat']));
fileIDExtraction = struct2cell(dir(['/Users/darrentan/Documents/SREP/Simulations/Data collection/sim14/Extractions' '/*.mat']));

% Initialising the final struct
compiledTotal = {};

% Loop through every pair of simulation and extraction
for num = 1:size(dir(['/Users/darrentan/Documents/SREP/Simulations/Data collection/sim14/Simulations' '/*.mat']), 1)
    similarity = {};    % Creating a new similiarty cell array per pair
    
    % Loading in the simulation and extraction files.
    [simulated, extracted] = loadSimulations(fileIDSimulation(1:2,num), fileIDExtraction(1:2,num));
    % Plot the simulated sp locations (from chris code) for visualisation
    plotSpkLoc(simulated, extracted);

    % Note: if the buffer is tooooo large, there can be doubles in matchExtr,
    % and this will cause an ERROR
    buffer = 300;   % VARY THIS VALUE

    % Initialisation of similarity values
    similarity.totalExactMatched    = 0;
    similarity.totalBufferedMatched = 0;

    % Comparison code
    % Loop through every template
    for templ = 1:extracted.nTemplates
        % Loop through every family per template
        for fam = 1:extracted.nFamilies(templ)
            flag = 0;
            % Initialise the similarity arrays
            similarity.numMatched{templ}{fam}           = [];
            similarity.exactMatches{templ}{fam}         = [];
            similarity.exactMatcheswBuffer{templ}{fam}  = [];

            % Per template, per family, loop through the number of axons
            % simulated
            n = 1;  % Initialising n variable
            while n <= simulated.nAxons

                % Find where there is a match according to the initial buffer
                [matchSim, matchExtr]       = ismembertol(simulated.spkIndex{n}, extracted.spks{templ}{fam}, buffer, 'DataScale', 1);
                % Finds where in the simulation there is a match
                matchSimInd                 = find(matchSim == 1);

                % Find the values where the difference is non consecutive to
                % remove from the template/family array
                % Removes from the extraction index array where there is not
                % match between extracted and simulated datasets
                matchExtr(matchExtr == 0)   = [];

                % Creates a temp array with values from 1 to size of family
                % extracted in increments of 1. This is then compared to
                % matchExtr to find where the 'indexing' is not consecutive
                % (e.g. 1 2 3 4 6 7 is missing 5, hence 5 is the non
                % consecutive value/index).
                % The value corresponding to this index from the extracted
                % family array is then removed
                temp                        = [1:size(extracted.spks{templ}{fam}, 1)]';
                nonConsec                   = ismember(temp, matchExtr);
                extractedFamComp            = extracted.spks{templ}{fam};
                extractedFamComp(find(nonConsec == 0)) = [];

                % If there is a difference in size between the two arrays,
                % there will be an error printed. This will be due to the
                % initial buffer being too large so the warning message will
                % indicate to user to decrease the initial buffer value
                try
                    % Calculate the overall shift between extracted and simulated
                    shift = extractedFamComp - simulated.spkIndex{n}(matchSimInd);
                catch ME
                    if (strcmp(ME.identifier, 'MATLAB:dimagree'))
                        buffer  = buffer - 10;   % Reduce the buffer value
                        flag    = 1;
                        continue
                    end
                end  
                
                % First need to test if the buffer value has simply just
                % been increasing, if it has, reset it back to its max and
                % change flag value.
                % Intuitively, if buffer needs to be above its MAX value
                % (set as 300, can be varied), its most likely not a match
                % If buffer is less than 300 and flag value is 0, then will
                % increase the buffer by 10 to check if there has been a
                % deviation of shift
                if buffer > 300
                    buffer  = 200;
                    flag    = 1;
                elseif n == 1 && flag == 0
                    buffer = buffer + 10;
                    continue
                end

                % Depending on size of the array per template/family, there
                % will be different functions
                % Essentially, it will find the shift with the most occurances,
                % indicating that for this template, this is the overall shift.
                % For spikes that have this shift, it is an EXACT match. The
                % shiftBuffer is to give a small leaway for slight inaccuracies
                % that may occur.
                shiftBuffer = 2;
                % For arrays that are larger than 10 values, it will simply
                % find the shift value that occurs the most and this is the
                % offset
                if size(shift, 1) > 10  % NEEDS FINETUNING
                    shiftVal            = mode(shift);
                    shiftSame           = find(shift == shiftVal);
                    shiftSamewBuffer    = find(shift <= shiftVal + shiftBuffer & shift >= shiftVal - shiftBuffer);
                % For arrays that are smaller (between 3 to 10 values), 
                % This needs to be fixed, if max is only 1, indication that its
                % not picking up anything, will want to assume that there is no
                % match
                elseif (size(shift, 1) <= 10) && (size(shift, 1) >= 2)    % NEEDS FINETUNING
                    shiftUnique             = unique(shift);
                    shiftUniqueCount        = histc(shift(:), shiftUnique);
                    [~, ind]                = max(shiftUniqueCount);
                    if shiftUniqueCount(ind) == 1   % If the maximum repeat is only once, consider that this array is just full of random, non-matching spikes
                        shiftVal            = NaN;
                        shiftSame           = [];
                        shiftSamewBuffer    = [];
                    else
                        shiftVal            = shiftUnique(ind);
                        shiftSame           = find(shift == shiftVal);
                        shiftSamewBuffer    = find(shift <= shiftVal + shiftBuffer & shift >= shiftVal - shiftBuffer);
                    end
                % At this point, the array is probably empty, or just one
                % value. Ignore it basically as it is not significant.
                % However, there could be a way to check to see if the shift of
                % this single value corresponds to the shift in arrays > 10.
                else
                    shiftVal                = NaN;
                    shiftSame               = [];
                    shiftSamewBuffer        = [];
                end

                % Determines number of exact matches and buffered matches
                similarity.totalExactMatched                = similarity.totalExactMatched + size(shiftSame, 1);
                similarity.totalBufferedMatched             = similarity.totalBufferedMatched + size(shiftSamewBuffer, 1);
                similarity.exactMatches{templ}{fam}         = [similarity.exactMatches{templ}{fam}; size(shiftSame, 1)];
                similarity.exactMatcheswBuffer{templ}{fam}  = [similarity.exactMatcheswBuffer{templ}{fam}; size(shiftSamewBuffer, 1)];
                n = n + 1;
            end
            buffer = buffer + 10;
        end
    end

    % Prints out the report for overarching values
    report           = dispReport(extracted, simulated, similarity, num);

    % Compiles all the data into an array to export as an excel file
    compiledTotal    = arrangeIntoCells(extracted, simulated, similarity, report, num, compiledTotal);
    fprintf('\tTo access data, please open <strong>''compiled''</strong> in the workspace\n\n');
end

% Saves the file as an excel spreadsheet
saveFile(compiledTotal);

%% Functions

%% Load function

% Function will open the simulation files and extraction files and extract
% the required data into structs to be used by the comparison tool.
% Will automatically format arrays to the correct orientations.
function [simulated, extracted] = loadSimulations(simulation, extraction)

% Open simulation. Use if you want to manually select files
% [file, path] = uigetfile('*.mat');
% fullName = fullfile(path, file);
% load(fullName);

% Currently will automatically load the simulation files according to the
% file path
% load('/Users/darrentan/Documents/SREP/Simulations/Completed 1/simulation1.mat');
% load('/Users/darrentan/Documents/SREP/Simulations/Test 10/simulation10.mat');
fileName = fullfile(simulation(2), simulation(1));
load(fileName{1});

simulated.spks      = simulation.report.spks;                   % Extracts the spikes data
simulated.nAxons    = size(simulation.report.spks, 2);          % Number of axons that were simulated

for i = 1:simulated.nAxons
    simulated.spkIndex{i,:} = (find(simulated.spks(:,i) == 1)); % Loops through simulated axons and extracts the indexs where there were spikes
end

% Open extracted spikes. Use if you want to manually select files
% [file, path] = uigetfile('*.mat');
% fullName = fullfile(path, file);
% load(fullName);

% Currently will automatically load the extraction files according to the
% file path
% load('/Users/darrentan/Documents/SREP/Simulations/Completed 1/extraction1.mat');
% load('/Users/darrentan/Documents/SREP/Simulations/Test 10/extraction10.mat');
fileName = fullfile(extraction(2), extraction(1));
load(fileName{1});

extracted.dt               = extraction.dt;                         % Extracts the sampling rate

if size(extraction.APstimes, 2) > 1
    extraction.APstimes = extraction.APstimes';                     % Corrects formatting of cell if detected to be incorrect
end
extracted.nTemplates       = size(extraction.APstimes, 1);          % Number of templates

for i = 1:extracted.nTemplates
    extracted.nFamilies(i) = size(extraction.APstimes{i, 1}, 1);    % Number of families per template
end

extracted.spks             = indexConversion(extraction.APstimes, extracted);   % Determines the extracted indexing
end

%% Converts extracted from being indexed by time to index

% The extracted data from the Extraction tool originally will record spikes
% by simulation time, hence convertion of time to indexes is required
function index = indexConversion(times, extracted)

% Initialisng empty index cell
index = {};

for templ = 1:extracted.nTemplates          % Loop through every template
    for fam = 1:extracted.nFamilies(templ)  % Loop through every family per template
        % If statement is to account for different formatting of cell
        if size(times{templ,1}{fam, 1}, 2) > 1
            index{templ, 1}{fam, 1} = round(times{templ, 1}{fam, 1}'/extracted.dt);  % Converts times to indexes of the extracted data
        else
            index{templ, 1}{fam, 1} = round(times{templ, 1}{fam, 1}/extracted.dt);   % Converts times to indexes of the extracted data
        end
    end
end

end

%% Plotting the simulated spike and extracted locations

% Plotting the simulated and extracted spikes to visualise what and where
% is being matched and determining if there are errors
function plotSpkLoc(simulated, extracted)

figure;
hold on;

% Plotting the simulated
for nAxon = 1:simulated.nAxons
    nAxonTier = ((nAxon - 1)/length(simulated.spkIndex)) * 0.1;     % Scaling for seperate axons being simulated for visualisation
    stem(simulated.spkIndex{nAxon}, ones(size(simulated.spkIndex{nAxon})) - nAxonTier);
end

% Plotting the extracted
for templ = 1:extracted.nTemplates
    for fam = 1:extracted.nFamilies(templ)
        templTier   = (templ - 1) * 0.1;                            % Scaling for each template for visualisation
        famTier     = (fam / size(extracted.spks{templ}, 1)) * 0.1; % Scaling for each family for visualisation
        
        stem(extracted.spks{templ}{fam}(:,1), 0.5*ones(size(extracted.spks{templ}{fam}(:,1))) - famTier - templTier);
    end
end
end

%% Generate report
function report = dispReport(extracted, simulated, similarity, num)

% Find the total number of spikes simulated in the simulation file
for i = 1:simulated.nAxons
    report.numPerSimulatedAxon(i) = size(find(simulated.spks(:,i) == 1), 1);
end
report.numSimulatedSpks           = sum(report.numPerSimulatedAxon);

% Find the total number of spikes extracted by the extraction tool
for templ = 1:extracted.nTemplates
    sumTempl = 0;
    for fam = 1:extracted.nFamilies(templ)
        sumTempl = sumTempl + size(extracted.spks{templ, 1}{fam, 1}, 1);
    end
    report.numPerExtractedTemplate(templ) = sumTempl;
end
report.numExtractedSpks = sum(report.numPerExtractedTemplate);

% Print the number of spikes for each and include the percentage that the
% extractor was able to pick up
fprintf('\n\t<strong>Comparison number %d</strong>\n', num);
fprintf('\tNumber of Simulated spikes: %.2f\n', report.numSimulatedSpks);
fprintf('\tNumber of OBTAINED Extracted spikes: %.2f\n', report.numExtractedSpks);
fprintf('\tNumber of CORRECTLY Extracted spikes: %.2f\n', similarity.totalBufferedMatched);
fprintf('\t%.2f%% of spikes were extracted\n\n', (similarity.totalBufferedMatched/report.numSimulatedSpks)*100);   % This value is according to the buffered matched

end

%% Excel

% This function compiles all the data that has been determined into a cell
% array to be exported into excel
function compiledTotal = arrangeIntoCells(extracted, simulated, similarity, report, num, compiledTotal)

compiled = cell(1, 7);

compiled{3,1} = 'Comparison number:'; compiled{3,2} = num;

% Storing overall statistics
compiled{4,1} = 'Number of Simulated Spikes:'; compiled{4,2} = report.numSimulatedSpks;
compiled{4,3} = 'Number of OBTAINED Extracted Spikes:'; compiled{4,4} = report.numExtractedSpks;
compiled{4,5} = 'Number of CORRECTLY Extracted Spikes:'; compiled{4,6} = similarity.totalBufferedMatched;
% compiled{2,7} = '';

count = 5;  % Count value for formatting purposes

% Giving titles per column
compiled{5,1} = 'nAxon';                % Simulated axon number
compiled{5,2} = 'Template';             % Extracted template number
compiled{5,3} = 'Family';               % Extracted family number
compiled{5,4} = 'Exact Matched';        % Number of EXACT matches between extracted and simulated
compiled{5,5} = 'Buffered Matched';     % Number of BUFFERED matches between extracted and simulated
compiled{5,6} = 'Extracted Total';      % Total number of extracted spikes per template-family
compiled{5,7} = 'Certainty';            % Certainty that the template-family is from the simulated axon as a percentage (Buffered Matched value is used here)

for axon = 1:simulated.nAxons
    for templ = 1:extracted.nTemplates
        for fam = 1:extracted.nFamilies(templ)
            % Providing data according to the column title
            compiled{count + 1, 1} = axon;          % Simulation axon number
            compiled{count + 1, 2} = templ;         % Extracted template number
            compiled{count + 1, 3} = fam;           % Extracted family number
            compiled{count + 1, 4} = similarity.exactMatches{templ}{fam}(axon);         % Number of EXACT matches between extracted and simulated
            compiled{count + 1, 5} = similarity.exactMatcheswBuffer{templ}{fam}(axon);  % Number of BUFFERED matches between extracted and simulated
            compiled{count + 1, 6} = size(extracted.spks{templ}{fam}, 1);               % Total number of extracted spikes per template-family
            compiled{count + 1, 7} = [num2str(round((similarity.exactMatcheswBuffer{templ}{fam}(axon)...    % Certainty that the template-family is from the simulated axon as a percentage (Buffered Matched value is used here)
                /size(extracted.spks{templ}{fam}, 1))*100, 2)), 2 '%'];
            count                  = count + 1;
        end
    end
end

compiledTotal = vertcat(compiledTotal, compiled);

end

%% Safe file

% Function saves the compiled data as an .xlsx file to be accessed in excel
function saveFile(compiled)

[file, path] = uiputfile(['completedExtractionInformation' filesep 'completedExtraction.xlsx'], 'Save file name');
if file
    fileName = [path file];
    try
        writecell(compiled, fileName);
    catch E
        str = sprintf('\t Excel file couldn''t be saved. Please try again. \n\tError; %s\n\n', E.message);
        runtimeErrorHandler(E, 'message', str);
    end
else
    fprintf('\tUser didn''t choose a file location. The Excel file wasn''t saved.\n\n');
end

end