% Script to compare extractor tool spike times with the simulated times

close all
clear

[simulated, extracted] = loadSimulations;

buffer = 50;
similarity = {};
matchedIndex{simulated.nAxons} = [];
totalCount = 0;

for templ = 1:extracted.nTemplates
    for fam = 1:extracted.nFamilies(templ)
        similarity{templ}{fam} = [];
        for n = 1:simulated.nAxons
            simCount = 0;
            for currentSpk = extracted.spks{templ}{fam}
%                 match = find((simulated.spkIndex{n}(:) <= extracted.spks{templ}{fam} + buffer) & (simulated.spkIndex{n}(:) >= extracted.spks{templ}{fam} - buffer));     % Matching to the simulation indexs
%                 match = find((extracted.spks{templ}{fam}(:) <= simulated.spkIndex{n} + buffer) & (extracted.spks{templ}{fam}(:) >= simulated.spkIndex{n} - buffer));
                match = find((currentSpk <= simulated.spkIndex{n} + buffer) & (currentSpk >= simulated.spkIndex{n} - buffer));
                if match; simCount = simCount + 1; totalCount = totalCount + 1; end
                matchedIndex{n} = [matchedIndex{n}; match];
%                 remove = ismember(simulated.spkIndex{n}, matched);
%                 simulated.spkIndex{n}(remove) = [];
            end
            % Calculates the similarity of the extracted to the simulated
            % axon
            similarity{templ}{fam} = [similarity{templ}{fam} (simCount/size(extracted.spks{templ}{fam}, 1))];
        end
    end
end
report = comparisonReport(simulated, extracted, similarity);

% buffer = 10;
% 
% % First 
% differences1 = abs(d_simulated - d_extracted(1));
% inToleranceIndexes1 = differences1 < buffer;
% withinTolerance1 = find(inToleranceIndexes1 == 1); %d_simulated(inToleranceIndexes1);
% 
% % Second
% differences2 = abs(d_simulated - d_extracted(2));
% inToleranceIndexes2 = differences2 < buffer;
% withinTolerance2 = find(inToleranceIndexes2 == 1);


% found_first = find(d_simulated == d_extracted(1), 1, 'first');
% found_second = find(d_simulated == d_extracted(2));



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
        sumTempl = sumTempl + size(extracted.spks{templ, 1}{fam, 1}, 1);
    end
    report.numPerExtractedTemplate(templ) = sumTempl;
end
report.numExtractedSpks = sum(report.numPerExtractedTemplate);

% Print the number of spikes for each and include the percentage that the
% extractor was able to pick up

fprintf('Number of Simulated spikes: %.2f\n', report.numSimulatedSpks);
fprintf('Number of Extracted spikes: %.2f\n', report.numExtractedSpks);
fprintf('%.2f%% of spikes were extracted\n', (report.numExtractedSpks/report.numSimulatedSpks)*100);

% Need to report the similarity between the templates




end

function index = indexConversion(times, extracted)

index = {};

for templ = 1:extracted.nTemplates
    for fam = 1:extracted.nFamilies(templ)
        index{templ, 1}{fam, 1} = round(times{templ, 1}{fam, 1}/extracted.dt);
    end
end

end