% SPIKESIMULATIONTOOL Script to simulate a spike train from single channel extracellular recording from various cells and drift*.
%
% *Drift is the growth of spike amplitude and noise level at a different
%  rate.
%
% University of Melbourne
% Department of Biomedical Engineering
% Artemio Soto-Breceda | 6/August/2019

%% Copyright 2019 Artemio Soto-Breceda
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.


function vsim = SpikeSimulationTool(varargin)
%% Options
Naxons      = 1;                    % Number of different templates
SNR         = 20;                   % Initial signal to noise ratio (it will change with drift)
growth      = [(1.1 + (2 - 1.1) * rand) (0.8 + (1.2 - 0.8) * rand)]; % [1.9 1.1];   % Growth of: [<spamp> <noise>]
total_time  = 20;                 % Seconds
fs          = 5000;                 % Sampling rate (current template file has this sampling rate, so it should stay like this unless the templates are fixed)
sr          = randi(10,1,Naxons)/2; % Spike rate
overlap     = false;                % If true, it allows spikes of diff axons to overlap
rpt_temp    = true;                 % True if more than one axon can have the same template
pw_locs     = round((total_time/4 + ((total_time*3/4) - (total_time/4)) * rand) * fs); % Locations of the change of drift. Only controls the growth of the spike amplitudes, not noise
pw_grow     = 1.2 + (1.5 - 1.2) * rand; % The new drift for each location. Only controls the growth of the spike amplitudes, not noise
has_drift   = true; 
has_noise   = true;
pre_noise   = true;                 % Append period of just noise at the start
do_filter   = true;                 % Bandpass filter the signal
passband    = [40 1200];%[80, 600]; % Passband
PLOT        = false;
templates   = [];                   % Initialize templates
templates_fn= 'templates_test_struct.mat'; % Default templates file

if nargin >= 1
   data       = varargin{1};
   
   Naxons     = data.Naxons;
   SNR        = data.SNR;
   total_time = data.total_time;
   pw_locs    = round((total_time/4 + ((total_time*3/4) - (total_time/4)) * rand) * fs); % Locations of the change of drift. Only controls the growth of the spike amplitudes, not noise
   fs         = data.fs;
   sr         = data.sr;
   overlap    = data.overlap;
   rpt_temp   = data.rpt_temp;
   has_drift  = data.has_drift;
   has_noise  = data.has_noise;
   pre_noise  = data.pre_noise;
   do_filter  = data.do_filter;
   passband   = data.passband;
   PLOT       = data.PLOT;
   templates_fn= data.templates_fn;
   
   swtch      = varargin{2};    % Pass in swtch struct for evnts parameters
   events     = varargin{3};    % Pass in events struct for evnts parameters
end

%% Events
evnts.inflammation_axons   = floor(0 + ((Naxons/2 - 0) * rand)); % Number of inflamed axons (increase the spike rate).
evnts.inflammation_onset   = round((total_time/4 + ((total_time*3/4) - (total_time/4)) * rand) * fs);  % High frequency at time
evnts.inflammation_tau     = 5e-3*fs;  % Time constant for increased spike rate to decay to spontaneous activity

% Natural
evnts.amplitude_nat_onset  = 500*fs;   % Change of amplitude in just some of axons
evnts.amplitude_nat_axons  = 0;        % Change of amplitude in just a couple of axons


evnts.amplitude_dist_onset = round((total_time/4 + ((total_time*3/4) - (total_time/4)) * rand) * fs);  % Change of amplitude in all axons
evnts.amplitude_dist_value = 0.5 + (1.5 - 0.5) * rand;      % Value of the new amplitude multiplier
evnts.amplitude_dist_prob  = 0.2; % Probability of having a change in the amplitude

evnts.prob_start           = floor(0 + ((Naxons/2 - 0) * rand)); % (Recruited) Number of axons that don't start at the beginning. They will randomly start somewhere along the recording.
evnts.prob_end             = floor(0 + ((Naxons/2 - 0) * rand)); % (Dismissed) Number of axons that don't last the whole recording. They will randomly end somewhere along the recording.

% If there is change in event value, input values
if nargin >=1
    for fn = fieldnames(events)'
        if swtch.(fn{1}) == 1
            evnts.(fn{1}) = events.(fn{1});
        end
    end
end

%% Run
% Load the templates matrix

if isempty(data.templates)
    % Find the directory path to templates folder
    [path,~,~] = fileparts(mfilename('fullpath'));
    load([path filesep 'templates' filesep templates_fn]);
else
    templates = data.templates; % Loads generated template from templatesApp.mlapp
end


dt = 1/fs;
try % Generate a train of extracellular spikes. There is no noise
   [v, vv, report] = gen_train(templates, Naxons, fs, total_time/dt,       ...
                               'SpikeRate',       sr,              ...
                               'Overlap',         overlap,         ...
                               'Recruited',       evnts.prob_start,...
                               'Dismissed',       evnts.prob_end,  ...
                               'Events',          evnts,           ...
                               'RepeatTemplates', rpt_temp         ...
                               );
catch E
   if strcmp('Manually stopped', E.message)
      fprintf(2,'\tManually stopped\n');
      return;
   else
      rethrow(E);
   end
end

% Add drift
if has_drift
   % Add the noise and drift
   v = add_drift(v,                           ...
                 'SNR',            SNR,       ...
                 'Noise',          has_noise, ...
                 'Growth',         growth,    ...
                 'PrecedingNoise', pre_noise, ...
                 'Linear',         false,     ...
                 'PwLocs',         pw_locs,   ...
                 'PwGrowth',       pw_grow    ...
                 );
elseif has_noise
   % Add only noise
   v = add_drift(v,                           ...
                 'SNR',            SNR,       ...
                 'Noise',          has_noise, ...
                 'Growth',         [1 1],     ...
                 'PrecedingNoise', false      ...
                 );
end

if do_filter && has_noise
   % Lowpass filter the signal
   try
    v = bandpass(v, passband, fs);
   catch E
       if strcmp('MATLAB:UndefinedFunction', E.identifier)
           str = sprintf('\tYou are using a Matlab version prior to 2018a, the simulation won''t be filtered.\n');
           cprintf('Errors', str);
       else
           rethrow(E);
       end
   end
end


%% Format recording for SpikeExtractionTool
vsim        = struct;
vsim.type   = 'voltage';
vsim.name   = 'vsim';
vsim.dt     = dt;
vsim.params = [];
vsim.data   = v;
vsim.time   = (dt:dt:length(vsim.data)*dt);
vsim.axons  = vv;
vsim.report = report;
[~, vsim.srt] = sort(max(vv), 'descend');
% Fix the dimensions if they are wrong
if size(vsim.data,1) < size(vsim.data,2), vsim.data = vsim.data'; end
if size(vsim.time,1) < size(vsim.time,2), vsim.time = vsim.time'; end

%% Plot Figure
if PLOT
   figure;
   plot(vsim.time, vsim.data);
   xlabel('Time (s)');
   ylabel('Amplitude');
   
   figure;
   vv = vsim.axons;
   max_vv = max(vv);
   %[~, vsim.srt] = sort(max_vv, 'descend');
   try
      if ~has_drift || ~pre_noise
         plot(vsim.time, vv(:,vsim.srt));
      else
         plot(vsim.time(101:end), vv(:,vsim.srt));
      end
   catch
      plot(vsim.time(101:end), vv(:,vsim.srt));
   end
   xlabel('Time (s)');
   ylabel('Amplitude');
   naxons = size(vsim.axons, 2);
   % Plot templates
   figure;
   for temp = 1:naxons
      nr = ceil(sqrt(naxons));
      nc = ceil(naxons/3);
      if naxons == 3
         nc = 1; nr = 3;
      end
      subplot(nc, nr, temp);
      plot(vv( max(vsim.report.locs{temp}(1) - 100, 1) : min(vsim.report.locs{temp}(1) + 100, size(vv,1)), temp));
   end
   clear vv
end

%% Print report
%printReport(report, dt);
%% Save
% Check size of the variable to save. If it is larger than 2Gb, remove the
% per axon field: vsim.axons

s = whos('vsim');
if s.bytes > 2e9
   vsim.axons = [];
   fprintf('\tThe file is too large, only the final recording will be saved, not the per-axon information.\n');
end
% 
% [file,path] = uiputfile(['simulations' filesep 'sim.mat'],'Save file name');
% if file
%    file_name = [path filesep file];
%    save(file_name, 'vsim');
% else
%    fprintf('\tUser didn''t chose a file location. The simulation wasn''t saved.\n');
% end

% sufix = 0;
% file = 'sim';
% valid_file = file;
% while exist([valid_file,'.mat'], 'file')
%    sufix = sufix + 1;
%    valid_file = [file num2str(sufix)];
% end
% file = valid_file;
% save(file, 'vsim');
