% Script to compare extractor tool spike times with the simulated times

close all
clear

[simulated, extracted] = loadSimulations;

buffer = 30;
similarity.numMatched{extracted.nTemplates}{extracted.nFamilies} = [];
% similarity.matched{extracted.nTemplates}{extracted.nFamilies} = [];
matchedIndex{extracted.nTemplates}{extracted.nFamilies} = [];

for templ = 1:extracted.nTemplates
    for fam = 1:extracted.nFamilies(templ)
        for n = 1:simulated.nAxons  % Iterating through nAxons to find if there are matches
            match = ismembertol(simulated.spkIndex{n}, extracted.spks{templ}{fam}, buffer, 'DataScale', 1);
            matchedIndex{templ}{fam} = [matchedIndex{templ}{fam}; find(match == 1)];    % The index values for this axon will relate back to the SIMULATED index. FIX MATRIX DIMSNSIONS (TRANSFORM?)
        end
        % Determining the num of matches per temp-fam to nAxon
        similarity.numMatched{templ}{fam} = [similarity.numMatched{templ}{fam} size(matchedIndex{templ}{fam}, 2)];
        % Determining the % of matches per temp-fam to nAxon
%         similarity.matched{templ}{fam} = [similarity.matched{templ}{fam} (size(matchedIndex{templ}{fam}, 2)/size(extracted.spks{templ}{fam}, 1))];
    end
end

report = comparisonReport(simulated, extracted, similarity);

%% Functions

function [simulated, extracted] = loadSimulations

% Open simulation
% [file, path] = uigetfile('*.mat');
% fullName = fullfile(path, file);
% load(fullName);

load('/Users/darrentan/Documents/SREP/Simulations/Completed 1/simulation1.mat');
% simulated.simulated = simulation;
simulated.spks      = simulation.report.spks;
simulated.nAxons    = size(simulation.report.spks, 2);
for i = 1:simulated.nAxons
    simulated.spkIndex{i} = (find(simulated.spks(:,i) == 1))';
end

% Open extracted spikes
% [file, path] = uigetfile('*.mat');
% fullName = fullfile(path, file);
% load(fullName);

load('/Users/darrentan/Documents/SREP/Simulations/Completed 1/extraction1.mat');
% extracted.extracted  = extraction;
extracted.dt               = extraction.dt;
extracted.nTemplates       = size(extraction.APstimes, 1);
for i = 1:extracted.nTemplates
    extracted.nFamilies(i) = size(extraction.APstimes{i, 1}, 1);
end
extracted.spks             = indexConversion(extraction.APstimes, extracted);
end

function report = comparisonReport(simulated, extracted, similarity)

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

fprintf('\n Number of Simulated spikes: %.2f\n', report.numSimulatedSpks);
fprintf(' Number of Extracted spikes: %.2f\n', report.numExtractedSpks);
fprintf(' %.2f%% of spikes were extracted\n\n', (report.numExtractedSpks/report.numSimulatedSpks)*100);

% Need to report the similarity between the templates

fprintf('                  |'); fprintf(pad('Family', 14*extracted.nFamilies, 'both')); fprintf('\n');
for templ = 1:extracted.nTemplates
    str = [];
    fprintf(' <strong>-------------</strong>'); for fam = 1:extracted.nFamilies; fprintf('<strong>---------------</strong>'); end; fprintf('\n');
    fprintf('    Template');
    for fam = 1:extracted.nFamilies; fprintf('      |      %d', fam); end; fprintf('\n');
    for fam = 1:extracted.nFamilies; str = [str '|' pad(num2str(similarity.numMatched{templ}{fam}), 6, 'left') '/' pad(num2str(size(extracted.spks{templ}{fam}, 1)), 6, 'right')]; end
    fprintf('        %d         ', templ); disp(str);
end

end

function index = indexConversion(times, extracted)

index = {};

for templ = 1:extracted.nTemplates
    for fam = 1:extracted.nFamilies(templ)
        index{templ, 1}{fam, 1} = round(times{templ, 1}{fam, 1}/extracted.dt);  % THERE CAN BE VARYING MATRIX DIMENSIONS (ROW OR COLUMN MATRIX)
    end
end

end