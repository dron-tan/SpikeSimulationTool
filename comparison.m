% comparison script which is used to test the validity of Katies
% SpikeExtractionTool

% University of Melbourne
% Department of Biomedical Engineering
% Summer Research Experience Program - Spike Simulation Suite
% Christopher Ong & Chi Yung Darren Tan | 08 March 2021

%%
% The comparison tool will automate testing the accuracy of the extraction
% tool

% Initially, user must generate extraction and simulation pairs with files
% name being:
% 'simulation#' and 'extraction#' where # is replaced by a integer

% comparison tool can then be run to test. The buffer will need to be
% adjusted

close all
clear

% Loading in the simulation and extraction files NEED TO AUTOMATE
[simulated, extracted] = loadSimulations;

% Plot the simulated sp locations (from chris code)
plotSpkLoc(simulated, extracted);

% Initialising values
% Note: if the buffer is tooooo large, there can be doubles in matchExtr,
% and this will cause an ERROR

buffer = 50;   % VARY THIS VALUE
fprintf('\n\t<strong>Buffer used = %d</strong>\n', buffer);

similarity.totalExactMatched = 0;
similarity.totalBufferedMatched = 0;


% Loop through every template
for templ = 1:extracted.nTemplates
    % Loop through every family per template
    for fam = 1:extracted.nFamilies(templ)
        % Initialise the similarity struct
        similarity.numMatched{templ}{fam} = [];
        similarity.exactMatches{templ}{fam} = [];
        similarity.exactMatcheswBuffer{templ}{fam} = [];
        
        % Per template, per family, loop through the number of axons
        % simulated
        for n = 1:simulated.nAxons
            % Find where there is a match according to the buffer
            [matchSim, matchExtr] = ismembertol(simulated.spkIndex{n}, extracted.spks{templ}{fam}, buffer, 'DataScale', 1);
            matchSimInd = find(matchSim == 1);
            
            % Find the values where the difference is non consecutive to
            % remove from the template/family array
            matchExtr(matchExtr == 0) = [];
%             matchExtrNonConsec = find(diff(matchExtr) > 1) + 1;    % Want to keep this value, this is the value that is not matched (or just keep count of amount not matched)
            
            temp = [1:size(extracted.spks{templ}{fam}, 1)]';
            nonConsec = ismember(temp, matchExtr);
            extractedFamComp = extracted.spks{templ}{fam};
            extractedFamComp(find(nonConsec == 0)) = [];
            
            % Calculate the overall shift between extracted and simulated
            try
                shift = extractedFamComp - simulated.spkIndex{n}(matchSimInd);
            catch ME
                if (strcmp(ME.identifier, 'MATLAB:dimagree'))
                    warning('There are overlaping spikes. You will need to decrease the initial buffer value.');
                end
                rethrow(ME)
            end
            
            % Categoris
            shiftBuffer = 2;
            if size(shift, 1) > 10  % THIS VALUE NEEDS TO BE ADJUSTED
                shiftVal = mode(shift);
                shiftSame = find(shift == shiftVal);
                shiftSamewBuffer = find(shift <= shiftVal + shiftBuffer & shift >= shiftVal - shiftBuffer);
            elseif (size(shift, 1) <= 10) && (size(shift, 1) > 2)    % THESE RANGES CAN BE CHANGED
                %??????? IDK WHAT TO DO!!!!
                % Need to check if it is empty first, if empty, can just
                % skip
                shiftUnique = unique(shift);
                shiftUniqueCount = histc(shift(:), shiftUnique);
                [~, ind] = max(shiftUniqueCount);
                shiftVal = shiftUnique(ind);
                shiftSame = find(shift == shiftVal);
                shiftSamewBuffer = find(shift <= shiftVal + shiftBuffer & shift >= shiftVal - shiftBuffer);
            else
                % If it is not empty, will need to check if there is a mode
                % and determine if it is significant
                shiftVal = NaN;
                shiftSame = [];
                shiftSamewBuffer = [];
            end
            
            % Determines number of exact matches and buffered matches
            similarity.totalExactMatched = similarity.totalExactMatched + size(shiftSame, 1);
            similarity.totalBufferedMatched = similarity.totalBufferedMatched + size(shiftSamewBuffer, 1);
            similarity.exactMatches{templ}{fam} = [similarity.exactMatches{templ}{fam}; size(shiftSame, 1)];
            similarity.exactMatcheswBuffer{templ}{fam} = [similarity.exactMatcheswBuffer{templ}{fam}; size(shiftSamewBuffer, 1)];
        end
    end
end

report = dispReport(extracted, simulated, similarity);
compiled = arrangeIntoCells(extracted, simulated, similarity, report);
fprintf('\tTo access data, please open ''compiled'' in the workspace\n');
saveFile(compiled);

%% Functions

%% Load function

function [simulated, extracted] = loadSimulations

% Open simulation
% [file, path] = uigetfile('*.mat');
% fullName = fullfile(path, file);
% load(fullName);

% load('/Users/darrentan/Documents/SREP/Simulations/Completed 1/simulation1.mat');
load('/Users/darrentan/Documents/SREP/Simulations/Test 10/simulation10.mat');
% simulated.simulated = simulation;
simulated.spks      = simulation.report.spks;
simulated.nAxons    = size(simulation.report.spks, 2);
for i = 1:simulated.nAxons
    simulated.spkIndex{i,:} = (find(simulated.spks(:,i) == 1));
end

% Open extracted spikes
% [file, path] = uigetfile('*.mat');
% fullName = fullfile(path, file);
% load(fullName);

% load('/Users/darrentan/Documents/SREP/Simulations/Completed 1/extraction1.mat');
load('/Users/darrentan/Documents/SREP/Simulations/Test 10/extraction10.mat');
% extracted.extracted  = extraction;
extracted.dt               = extraction.dt;
if size(extraction.APstimes, 2) > 1
    extraction.APstimes = extraction.APstimes';
end
        
extracted.nTemplates       = size(extraction.APstimes, 1);
for i = 1:extracted.nTemplates
    extracted.nFamilies(i) = size(extraction.APstimes{i, 1}, 1);
end
extracted.spks             = indexConversion(extraction.APstimes, extracted);
end

%% Converts extracted from being indexed by time to index
function index = indexConversion(times, extracted)

index = {};

for templ = 1:extracted.nTemplates
    for fam = 1:extracted.nFamilies(templ)
        if size(times{templ,1}{fam, 1}, 2) > 1
            index{templ, 1}{fam, 1} = round(times{templ, 1}{fam, 1}'/extracted.dt);  % THERE CAN BE VARYING MATRIX DIMENSIONS (ROW OR COLUMN MATRIX)
        else
            index{templ, 1}{fam, 1} = round(times{templ, 1}{fam, 1}/extracted.dt);  % THERE CAN BE VARYING MATRIX DIMENSIONS (ROW OR COLUMN MATRIX)
        end
    end
end

end

%% Plotting the simulated spike and extracted locations
function plotSpkLoc(simulated, extracted)

figure;
hold on;

% Plotting the simulated
for nAxon = 1:simulated.nAxons
    nAxonTier = ((nAxon - 1)/length(simulated.spkIndex)) * 0.1;
    stem(simulated.spkIndex{nAxon}, ones(size(simulated.spkIndex{nAxon})) - nAxonTier);
end

% PLotting the extracted
for templ = 1:extracted.nTemplates
    for fam = 1:extracted.nFamilies(templ)
        templTier   = (templ - 1) * 0.1;
        famTier     = (fam / size(extracted.spks{templ}, 1)) * 0.1;
        
        stem(extracted.spks{templ}{fam}(:,1), 0.5*ones(size(extracted.spks{templ}{fam}(:,1))) - famTier - templTier);
    end
end
end

%% Generate report
function report = dispReport(extracted, simulated, similarity)

% Find the number of spikes for simulated and extracted
for i = 1:simulated.nAxons
    report.numPerSimulatedAxon(i) = size(find(simulated.spks(:,i) == 1), 1);
end
report.numSimulatedSpks           = sum(report.numPerSimulatedAxon);

for templ = 1:extracted.nTemplates
    sumTempl = 0;
    for fam = 1:extracted.nFamilies(templ)
        sumTempl = sumTempl + size(extracted.spks{templ, 1}{fam, 1}, 1);    % THERE CAN BE VARYING MATRIX DIMENSIONS (ROW OR COLUMN MATRIX)
    end
    report.numPerExtractedTemplate(templ) = sumTempl;
end
report.numExtractedSpks = sum(report.numPerExtractedTemplate);

% Print the number of spikes for each and include the percentage that the
% extractor was able to pick up

fprintf('\n\tNumber of Simulated spikes: %.2f\n', report.numSimulatedSpks);
fprintf('\tNumber of OBTAINED Extracted spikes: %.2f\n', report.numExtractedSpks);
fprintf('\tNumber of CORRECTLY Extracted spikes: %.2f\n', similarity.totalExactMatched);
fprintf('\t%.2f%% of spikes were extracted\n\n', (report.numExtractedSpks/report.numSimulatedSpks)*100);

end

%% Excel
function compiled = arrangeIntoCells(extracted, simulated, similarity, report)

compiled = cell(1, 7);

compiled{1,1} = 'Number of Simulated Spikes:'; compiled{1,2} = report.numSimulatedSpks;
compiled{1,3} = 'Number of OBTAINED Extracted Spikes:'; compiled{1,4} = report.numExtractedSpks;
compiled{1,5} = 'Number of CORRECTLY Extracted Spikes:'; compiled{1,6} = similarity.totalBufferedMatched;
compiled{1,7} = '';

count = 3;

compiled{3,1} = 'nAxon';
compiled{3,2} = 'Template';
compiled{3,3} = 'Family';
compiled{3,4} = 'Exact Matched';
compiled{3,5} = 'Buffered Matched';
compiled{3,6} = 'Extracted Total';
compiled{3,7} = 'Certainty';

for axon = 1:simulated.nAxons
    
    for templ = 1:extracted.nTemplates
        for fam = 1:extracted.nFamilies(templ)
            compiled{count + 1, 1} = axon;
            compiled{count + 1, 2} = templ;
            compiled{count + 1, 3} = fam;
            compiled{count + 1, 4} = similarity.exactMatches{templ}{fam}(axon);
            compiled{count + 1, 5} = similarity.exactMatcheswBuffer{templ}{fam}(axon);
            compiled{count + 1, 6} = size(extracted.spks{templ}{fam}, 1);
            compiled{count + 1, 7} = [num2str(round((similarity.exactMatcheswBuffer{templ}{fam}(axon)/size(extracted.spks{templ}{fam}, 1))*100, 2)), 2 '%'];
            count = count + 1;
        end
    end
end

end

%% Safe file
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