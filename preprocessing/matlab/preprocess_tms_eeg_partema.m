%% PREPROCESSING

clear all; close all; clc;

%% Set up paths
+
Using_Cluster = 0

if Using_Cluster

    FieldTrip_path = '/Volumes/home/rawad.alabboud/MATLAB/Matlab Software/fieldtrip-20231109';
    addpath(FieldTrip_path); ft_defaults;

    eeglab_path = '/Volumes/home/rawad.alabboud/MATLAB/Matlab Software/eeglab2025.0.0';
    addpath(eeglab_path); eeglab;

    Tesa_path = '/Volumes/home/rawad.alabboud/MATLAB/Matlab Software/TESA1.1.1';
    addpath(Tesa_path);

    External_functions = '/Volumes/home/rawad.alabboud/MATLAB/External_functions';
    addpath(External_functions);

    fastica_path = '/Volumes/home/rawad.alabboud/MATLAB/External_functions/FastICA_25';
    addpath(fastica_path);

else
    
    FieldTrip_path = '/Users/rawad.alabboud/Documents/MATLAB/Matlab Software/fieldtrip-20231109';
    addpath(FieldTrip_path); ft_defaults;

    eeglab_path = '//Users/rawad.alabboud/Documents/MATLAB/Matlab Software/eeglab2025.0.0';
    addpath(eeglab_path); eeglab;

    Tesa_path = '/Users/rawad.alabboud/Documents/MATLAB/Matlab Software/TESA1.1.1';
    addpath(Tesa_path);

    External_functions = '/Users/rawad.alabboud/Documents/MATLAB/External_functions';
    addpath(External_functions);

    fastica_path = '/Users/rawad.alabboud/Documents/MATLAB/External_functions/FastICA_25';
    addpath(fastica_path);
end
%% Experiment parameters

center = 'C1'; % C1: ETOILE / C2: BASTILLE / C3: APHP
subject = 'sub-001';
session = 'ses-1';
block = 'Stim'; % 'Stim' or 'PreStim_EC' or 'PreStim_EO' or 'PostStim_EC' or 'PostStim_EO'
stim = '01'; % 01 = 10 Hz left DLPFC

if Using_Cluster
    mainpath = '/Volumes/home/rawad.alabboud/DATA/PARTEMA Data/Raw Data/';
    save_path = ['/Volumes/home/rawad.alabboud/DATA/PARTEMA Data/Preprocessed Data/', subject];
else
    % Local GitHub paths for EEG data only
    mainpath = '/Users/rawad.alabboud/Documents/GitHub/partema-eeg-analysis/data/eeg/raw/';
    save_path = ['/Users/rawad.alabboud/Documents/GitHub/partema-eeg-analysis/data/eeg/preprocessed/', subject];
end
 
raw_data_path = [mainpath, subject, '/', session, '/eeg/'];
dataset = [center, '_', subject, '_', session, '_', block];
 
if ~exist(save_path, 'dir'); mkdir(save_path); end

extra_ICA = 0;
if extra_ICA
    fill_step = '_ICA3_';
else
    fill_step = '_';
end

triggercode = ['S TMS(1)'];   

n_pulses = 4;

%% Load EEG data from BrainVision .vhdr file
addpath(eeglab_path);

EEG = pop_loadbv([raw_data_path], [dataset,'.vhdr'], [], []); 

EEG = pop_chanedit(EEG, 'lookup', fullfile(eeglab_path, 'plugins', 'dipfit', 'standard_BESA', 'standard-10-5-cap385.elp'));
EEG = eeg_checkset(EEG);

eeglab redraw 

EEG_org = EEG; % save data in a copy struct

%% Remove specific channels by their names

% List of channel names to remove (example: {'Fpz', 'Fz', 'CP1', 'P5', 'AFz', 'Fp1'})
channels_to_remove = {'Fpz', 'Fz', 'CP1', 'P5', 'AFz', 'Fp1'};  

% Find indices of channels to remove
chan_idx_remove = find(ismember({EEG.chanlocs.labels}, channels_to_remove));

% Remove channels
EEG = pop_select(EEG, 'nochannel', chan_idx_remove);
EEG = eeg_checkset(EEG);

disp(['Removed channels: ', strjoin(channels_to_remove, ', ')]);

EEG = pop_select(EEG, 'nochannel', {'65', 'Sync', 'Gnd-Ref', 'Avg', 'Ref', 'Batt_0'});


%% Check for TMS pulses
for i = 1:length(EEG.event)
    if strcmp(EEG.event(i).type, 'S TMS(1)') && strcmp(EEG.event(i).code, 'Comment')
        EEG.event(i).code = 'Stimulus';
    end
end


TMSIndices = find(strcmp({EEG.event.type}, triggercode));
EEG_m = EEG.event(TMSIndices);
if strcmp(triggercode, 'S TMS(1)') && mod(length(EEG_m), 50) == 0  && all(strcmp({EEG_m.code}, 'Stimulus'))
    disp('No artefacts missing for TMS(1)');
else
    disp('Artefacts missing');
end

%% TMS Train detection and labeling (TMS1–TMS50)


if isempty(TMSIndices) || length(TMSIndices) < 2
    error('❌ Not enough TMS events to detect any trains.');
end

fs = EEG.srate;
max_isi_sec = 1;

% Convert latencies to time (s)
tms_latencies = [EEG.event(TMSIndices).latency] / fs;
dt = diff(tms_latencies);
train_starts = [1, find(dt > max_isi_sec) + 1];
train_ends   = [find(dt > max_isi_sec), length(tms_latencies)];

if max(train_starts) > length(TMSIndices) || max(train_ends) > length(TMSIndices)
    error('❌ Indexing error: train_start/end points exceed total TMS pulses.');
end

for i = 1:length(train_starts)
    pulse_range = TMSIndices(train_starts(i):train_ends(i));
    for p = 1:length(pulse_range)
        EEG.event(pulse_range(p)).type = ['TMS', num2str(p)];
    end
end

disp(['✅ ', num2str(length(train_starts)), ' trains labeled with TMS1–TMS50']);

%% Compute ISIs across 50 TMS pulses per train

nEpochs = length(train_starts);
TMS_lat = nan(nEpochs, 50);

for i = 1:nEpochs
    pulse_range = TMSIndices(train_starts(i):train_ends(i));
    latencies = [EEG.event(pulse_range).latency];
    TMS_lat(i, :) = latencies / fs * 1000;  % convert to ms
end

ISIs = diff(TMS_lat, 1, 2);  % [nEpochs x 49]

% Mimic original logic: use last ISI (between pulse 49 and 50)
TMS50_duration = mean(ISIs(:, end));
TMS50_sd = std(ISIs(:, end));

% Compute dist_cut as 3x TMS50 ISI + SD + buffer
buffer_ms = 100; % more suitable or high number of pulses... 50 and not 4 !
dist_cut = TMS50_duration * 49 + TMS50_sd + buffer_ms;
disp(['📏 Train cut length estimated from ISIs (TMS50*49 + SD + 40): ', num2str(dist_cut), ' ms']);

%% BASELINE CORRECTION AND VISUALIZATION

EEG = pop_epoch(EEG, {'TMS1'}, [-5.5 6], ...
    'newname', [dataset, '_TMS_train_epochs'], 'epochinfo', 'yes');

EEG = pop_rmbase(EEG, [-5000 -100]);

% Visual inspection
nChans = EEG.nbchan;
EEGLAB_plot_EEG_wrapper(EEG, [-499 5499], 20, 1:nChans);

clear data
data = squeeze(mean(EEG.data,3));

figure; pop_timtopo(EEG, [-100 100],  [], 'ERP data and scalp maps [-100 to 100 ms]');
figure; pop_timtopo(EEG, [-499 5499],  [], 'ERP data and scalp maps [-500 to 5500 ms]');
figure; pop_timtopo(EEG, [0 500],  [], 'ERP data and scalp maps [0 to 500 ms]');


%% MIRROR REPLACEMENT OF TMS TRAIN ARTIFACT (REVERSED MIRROR LOGIC)

train_len_ms = 5000;     % TMS train length
dist_cut = 5010;         % full length of contaminated segment
offset = dist_cut - train_len_ms;

% Define time windows
replace_start_time = -offset;                                 % e.g., -100
replace_end_time   = train_len_ms + offset;                   % e.g., 5100
replace_duration   = replace_end_time - replace_start_time;   % 5200 ms

mirror_start_time  = -train_len_ms - 3*offset;              % -5300
mirror_end_time    = -offset;                                 % -100
                                                              % duration is
                                                              % the same:
                                                              % 5200

% Validate EEG epoch covers full range
if EEG.times(1) > mirror_start_time || EEG.times(end) < replace_end_time
    error('⛔ Epoch time window does not cover required replacement or mirror range.');
end

% Sample indices
fs = EEG.srate;
time_vec = EEG.times;

idx_replace_start = find(time_vec >= replace_start_time, 1);
idx_replace_end   = find(time_vec >= replace_end_time, 1);
idx_mirror_start  = find(time_vec >= mirror_start_time, 1);
idx_mirror_end    = find(time_vec >= mirror_end_time, 1);

% Safe replace length
replace_len = min(idx_replace_end - idx_replace_start, idx_mirror_end - idx_mirror_start);
if replace_len <= 0
    error('⛔ Mirror or replacement range too short (length = %d)', replace_len);
end

% Extract and flip mirrored segment
mirror_segment = EEG.data(:, idx_mirror_start : idx_mirror_start + replace_len - 1, :);
flipped_segment = mirror_segment(:, end:-1:1, :);

% Apply mirrored segment to EEG
EEG.data(:, idx_replace_start : idx_replace_start + replace_len - 1, :) = flipped_segment;

% ✅ Log
fprintf('✅ Replaced signal from %.0f ms to %.0f ms with mirrored segment from %.0f ms to %.0f ms\n', ...
    replace_start_time, replace_start_time + replace_len/fs*1000, ...
    mirror_start_time, mirror_end_time);

% Fix times if mismatch
if length(EEG.times) ~= size(EEG.data, 2)
    EEG.times = linspace(EEG.xmin * 1000, EEG.xmax * 1000, size(EEG.data, 2));
    warning('⚠️ EEG.times was reset to match EEG.data (length = %d)', length(EEG.times));
end

% QA plot
% Butterfly plot (all channels overlaid)
figure;
plot(EEG.times, EEG.data(1:nChans, :, 1), 'Color', [0.2 0.2 0.2 0.3]);  % plot 1st epoch only
xlabel('Time (ms)');
ylabel('Amplitude (μV)');
title('Butterfly plot (all EEG channels, Epoch 1)');
xlim([-5500 5500]);
grid on;

hold on;
yl = ylim;
fill([replace_start_time replace_start_time + replace_len/fs*1000 ...
      replace_start_time + replace_len/fs*1000 replace_start_time], ...
     [yl(1) yl(1) yl(2) yl(2)], [1 0.9 0.9], 'EdgeColor', 'none');
uistack(findobj(gca, 'Type', 'patch'), 'bottom');

%% Remove TMS-evoked muscle activity (5001.5 to 5015 ms)... 15 ms after the last pulse in the train

% EEG_ptr = pop_tesa_removedata(EEG, [5001.5 5015], [], {'TMS1'});
EEG_ptr = EEG

% === Plot butterfly ERP ===
erp_avg = mean(EEG.data(1:nChans, :, :), 3);  % average over trials, exclude non-EEG chans
figure;
plot(EEG.times, erp_avg', 'Color', [0.2 0.2 0.2 0.3]);  % semi-transparent black
xlabel('Time (ms)'); ylabel('Amplitude (μV)');
title('Butterfly ERP with TMS1, TMS50, and replaced region');
xlim([-500 6000]); grid on; hold on;

% === Overlay event markers ===
% TMS1 is at 0 ms
xline(0, '--r', 'TMS1', 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'middle');

% TMS50: 49 ISIs after TMS1 → typically around 5000 ms
TMS50_time = 5000;  % Adjust if you measured ISIs dynamically
xline(TMS50_time, '--b', 'TMS50', 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'middle');

% === Shade mirror-replaced region ===
replace_start_time = -offset;       % e.g. -2
replace_end_time = 5001;            % you already did replacement up to here

yl = ylim;
fill([replace_start_time replace_end_time replace_end_time replace_start_time], ...
     [yl(1) yl(1) yl(2) yl(2)], [1 0.85 0.85], 'EdgeColor', 'none');
uistack(findobj(gca, 'Type', 'patch'), 'bottom');

% % === Optional: shade tail-cleaned area ===
% tail_start = 5001.5;
% tail_end   = 5020;
% fill([tail_start tail_end tail_end tail_start], ...
%      [yl(1) yl(1) yl(2) yl(2)], [0.85 1 0.85], 'EdgeColor', 'none');  % light green

EEG_first_rem  = EEG;

%% 
% === 1. Butterfly plot: averaged ERP across channels ===
erp_avg = squeeze(mean(EEG.data(1:nChans, :, :), 3));  % exclude non-EEG channels
figure;
plot(EEG.times, erp_avg', 'Color', [0.2 0.2 0.2 0.3]);  % semi-transparent lines
xlabel('Time (ms)'); ylabel('Amplitude (µV)');
title('Butterfly plot (averaged ERP across channels)');
xlim([-500 6000]); grid on;

% === 2. Subset preview of EEG data (raw signal shape) ===
data = squeeze(mean(EEG.data,3));  % [channels x time]
figure;
plot(data(:,100:10000)');
xlabel('Samples (subset)'); ylabel('Amplitude (µV)');
title('Channel-wise raw ERP preview (100–10000 samples)');

% === 3. Topographic ERP + scalp maps ===
% You can adjust time range to match your interest (e.g., post-train activity)
figure;
pop_timtopo(EEG, [dist_cut 5999], [], 'ERP + scalp maps (1000 ms post train)');


%%

% UNTIL HERE EVERYTHING IS GOOD AND I AM SURE OF IT


%% === Manual trial rejection using EEG_ptr ===

TMPREJ = [];  % clear any previous markings
idx_start = find(EEG_ptr.times >= dist_cut, 1);
data_post = EEG_ptr.data(:, idx_start:end, :);
time_post = EEG_ptr.times(idx_start:end);

% Visual inspection: mark segments to reject trials
eegplot(data_post, ...
    'srate', EEG_ptr.srate, ...
    'winlength', 5, ...
    'command', 'global TMPREJ; TMPREJ = TMPREJ;', ...  % ensures TMPREJ is captured
    'limits', [time_post(1) time_post(end)], ...
    'title', sprintf('Manual rejection: post-train window only (%.1f–%.1f ms)', time_post(1), time_post(end)));

disp('🕵️‍♀️ EEG_ptr displayed — mark segments to reject (ESC to close when done)');

%% === Convert TMPREJ → trial rejection vector
if exist('TMPREJ', 'var') && ~isempty(TMPREJ)
    [trialrej, ~] = eegplot2trial(TMPREJ(:,1:5), size(data_post,2), size(data_post,3));
    
    % Map back to full trial indices
    full_trialrej = false(1, size(EEG_ptr.data,3));
    full_trialrej(:) = trialrej;
else
    full_trialrej = false(1, size(EEG_ptr.data,3));
end

%%
EEG_ptr.BadTr = find(full_trialrej);

fprintf('🧹 Rejecting %d trials out of %d (from EEG_ptr)\n', ...
    length(EEG_ptr.BadTr), size(EEG_ptr.data, 3));

if ~isempty(EEG_ptr.BadTr)
    EEG_ptr = pop_rejepoch(EEG_ptr, EEG_ptr.BadTr, 0);
end


%% === Finalize: copy cleaned data to EEG
EEG = EEG_ptr;             % overwrite main EEG with cleaned copy
EEG_checkmist = EEG;       % backup for sanity tracking
disp('✅ Cleaned EEG_ptr assigned to EEG (ready for ICA or ERP)');
figure;
erp_avg = squeeze(mean(EEG.data(1:nChans, :, :), 3));
plot(EEG.times, erp_avg', 'k');
title('ERP after manual trial rejection');
xlim([-500 6000]);

figure; pop_timtopo(EEG, [-5499  5999],  [], 'ERP data and scalp maps of  resampled');

%% Check for recharging artifact
start_ms = 5000;
end_ms   = 5999.5;  % <- adjust upper bound to stay within EEG.times range

% Find sample indices
idx_start = find(EEG.times >= start_ms, 1);
idx_end   = find(EEG.times >= end_ms, 1, 'last');

% Check
if isempty(idx_start) || isempty(idx_end) || idx_end <= idx_start
    error('⛔ Invalid indices: EEG.times does not cover [%g–%g] ms', start_ms, end_ms);
end

% Slice data
data = squeeze(mean(EEG.data, 3));  % average over epochs
data_slice = data(:, idx_start:idx_end);
time_slice = EEG.times(idx_start:idx_end);

% Plot
figure;
plot(time_slice, data_slice(1:nChans,:)', 'Color', [0.2 0.2 0.2 0.3]);
xlabel('Time (ms)'); ylabel('Amplitude (µV)');
title(sprintf('ERP: %.1f–%.1f ms', start_ms, end_ms));
grid on;

%% Remove recharging artifact

% EEG = EEGLAB_remove_and_interpolate_recharging_artifact(EEG, 2); % 1 + single pulse, 2 = repetitive pulse
% 
% data = squeeze(mean(EEG.data,3));
% figure()
% plot(data(:,2000:22990)')
% figure; pop_timtopo(EEG, [-1190  1190],  [], 'ERP data and scalp maps of  resampled');
% figure; pop_timtopo(EEG, [-15 45],  [], 'ERP data and scalp maps of  resampled');
% figure; pop_timtopo(EEG, [150 210],  [], 'ERP data and scalp maps of  resampled');
% figure; pop_timtopo(EEG, [320 380],  [], 'ERP data and scalp maps of  resampled');
% figure; pop_timtopo(EEG, [485 550],  [], 'ERP data and scalp maps of  resampled');

%% Resample to 5Kz %if already at 5Kz, no need to do this

% EEG = pop_resample( EEG, 5000); %5000
EEG = pop_resample( EEG, 1000);

%% Prepare the data for ICA: First remove for a while the very worst channels!

EEG_evo = mean(EEG.data,3);
[~, sigmas] = DDWiener(EEG_evo);
figure;
plot(sigmas,'*')

%% labeling the very worst channels to not affect the ICA run
badC = find(sigmas > (median(sigmas) + 3*std(sigmas)));
goodC = setdiff(1:length(sigmas),badC);

EEG2ICA = pop_select( EEG, 'nochannel', badC);

%% Before applying ICA to reduce artifacts, check which is the maximum number of independent components via singular value decomposition


tmpdata=reshape(EEG2ICA.data,[size(EEG2ICA.data,1) size(EEG2ICA.data,2)*size(EEG2ICA.data,3)]); %%%
tmpdata=tmpdata-repmat(mean(tmpdata,2),[1,size(EEG2ICA.data,2)*size(EEG2ICA.data,3)]); %%%
tempCpca=svd(tmpdata); 
th=0.04;% this threshold is arbitrary but it works in most cases
figure,semilogy(tempCpca,'.-')% plot of the eigenvalues of the data matrix %%%
a=ylim;
hold on
Cpca=max(find(tempCpca>th));
line([Cpca Cpca],[a(1) a(2)],'Color','k','LineStyle','--');
hold off
'check whether the vertical dashed line correctly separates the eigenvalues close to zero from the ones much higher than zero, otherwise change the threshold th'

% if the vertical dotted line correctly separates the eigenvalues close to
% zero from the ones much higher than zero, Cpca correctly corresponds to
% the total number of independent components that can be estimated from the
% data matrix. if this is not the case, you have to change the threshold above or to manually correct Cpca

% %% 3. ICA - First round
% 
% Cpca = 20;
% 
% EEG.ICA_dims1 = Cpca;
% 
% EEG2ICA=ICA_analysis(EEG2ICA,EEG2ICA.data,Cpca);% run ICA using the runica method from EEGlab
% 
% EEGpart.ICA=EEG2ICA;
% EEG2ICA.comp2remove=[];
% 'ICA has been computed'
% 
% %% select and remove artefact-contaminated independent components (ICs)
% % this function opens a GUI that allows you to inspect single independent
% % components and to decide whether to remove or to keep them, looking at
% % the topography of ICs weights (bottom right), at the time course of the
% % ICs (bottom left), at the power spectrum of the ICs (center right), at
% % the single-trial time course of the ICs (center left), at the original
% % average TMS-evoked potentials (top left), at the average TMS-evoked
% % potentials after subtraction of the ICs that represent artefacts
% 
% %%
% 
% EEGviz = EEG2ICA;
% 
% % Get sample indices for 5000–5999.5 ms
% fs = EEGviz.srate;
% idx_start = find(EEGviz.times >= 5000, 1);
% idx_end   = find(EEGviz.times >= 5999.5, 1, 'last');
% 
% % Crop only the ERP display range (keep ICA weights)
% EEGviz.data  = EEGviz.data(:, idx_start:idx_end, :);
% EEGviz.pnts  = size(EEGviz.data, 2);
% EEGviz.times = EEGviz.times(idx_start:idx_end);
% EEGviz.xmin  = EEGviz.times(1)/1000;
% EEGviz.xmax  = EEGviz.times(end)/1000;
% 
% % Keep ICA info
% EEGviz.icaweights = EEG2ICA.icaweights;
% EEGviz.icasphere  = EEG2ICA.icasphere;
% EEGviz.icawinv    = EEG2ICA.icawinv;
% EEGviz.icaact     = [];  % let EEGLAB recompute if needed

%% === ICA on full epoch: -5000 to 5999.5 ms (with mirrored segment) ===

% === Step 1: Ensure full time window ===
t1 = -5000;
t2 = 5999.5;  % end of epoch window in ms

% Confirm epoch covers this range (should already be loaded)
% No cropping: we will apply ICA over the entire EEG dataset

% === Step 2: Perform baseline removal BEFORE stimulation ===
EEG = pop_rmbase(EEG, [-500 -100], []);  % baseline correction on pre-train window

% === Step 3: Optional resample for ICA speed ===
EEG_full = pop_resample(EEG, 1000);  % resample to 1000 Hz if desired

% === Step 4: Remove noisy channels before ICA ===
EEG_evo = mean(EEG_full.data, 3);
[~, sigmas] = DDWiener(EEG_evo);
badC = find(sigmas > (median(sigmas) + 3*std(sigmas)));
goodC = setdiff(1:length(sigmas), badC);
EEG2ICA = pop_select(EEG_full, 'nochannel', badC);

% === Step 5: Estimate ICA dimension via PCA ===
tmpdata = reshape(EEG2ICA.data, size(EEG2ICA.data,1), []);
tmpdata = tmpdata - mean(tmpdata,2);
tempCpca = svd(tmpdata);
th = 0.04;  % Threshold
% Cpca = 20

fprintf('⚙️ PCA estimate: Cpca = %d components\n', Cpca);

% === Step 6: Run ICA ===
EEG2ICA = ICA_analysis(EEG2ICA, EEG2ICA.data, Cpca);

% === Save ===
EEG_ica_full = EEG2ICA;

%% === Step 6: Review components ===
EEG_ica_full.comp2remove = [];

ICA_selectandremove(EEG_ica_full);

%%
EEG.data = EEG_ica_full.data;
%%
erp_avg = squeeze(mean(EEG.data, 3));
figure;
idx_start = find(EEG.times >= -5000, 1);
idx_end   = find(EEG.times >= 5999.5, 1, 'last');
time_slice = EEG.times(idx_start:idx_end);
plot(time_slice, erp_avg);  % erp_avg is [channels × timepoints]
xlabel('Time (ms)'); ylabel('Amplitude (µV)');
title('ERP after ICA component removal');
grid on;

%% Insert cleaned post-train window back into full EEG (only matching time window)
idx_start = find(EEG.times >= t1, 1);
idx_end = idx_start + size(EEG_ica_posttrain.data, 2) - 1;

if idx_end > size(EEG.data,2)
    error('⛔ The cleaned data window exceeds original EEG data size. Adjust indices.');
end

EEG.data(goodC, idx_start:idx_end, :) = EEG_ica_posttrain.data;
EEG.goodC = goodC;
EEG.badC = badC;

%% Step 8: Visualize updated data
EEGLAB_plot_EEG_wrapper(EEG, [-5000 5999], 20, 1:size(EEG.data,1));

%% Step 9: Perform baseline correction again
EEG = pop_rmbase(EEG, [-5000 -100], []);

plot_start = max(t1, EEG.times(1));
plot_end = min(t2, EEG.times(end));

if plot_start >= plot_end
    error('⛔ Selected plotting window [%.1f–%.1f ms] is invalid within EEG time range [%.1f–%.1f ms].', t1, t2, EEG.times(1), EEG.times(end));
end

figure;
pop_timtopo(EEG, [plot_start plot_end], [], sprintf('ERP data and scalp maps (%.0f–%.1f ms post-train)', plot_start, plot_end));

%% Step 10: Save current state
EEG_posttrain_cleaned = EEG;

%% === Save current EEG after ICA cleaning ===
EEG_presound = EEG;

%% === Integrate non-stationary sound + compute leadfield ===
EEG = pop_tesa_sound(EEG, 'lambdaValue', 0.1, 'iter', 10);

% Visualize averaged data
avg_data = squeeze(mean(EEG.data,3));
figure; pop_timtopo(EEG, [-5000 5999], [], 'ERP data and scalp maps (post sound integration)');
figure; pop_timtopo(EEG, [5250 5300], [], 'ERP data and scalp maps (zoomed window)');

EEG_post_sound = EEG;

%% === Apply band-pass filtering pre-SSPSIR ===
EEG = pop_tesa_filtbutter(EEG, 2, 200, 4, 'bandpass');

% Visualize after filtering
avg_data = squeeze(mean(EEG.data,3));
figure; plot(avg_data(:,100:2000)');

%% === Remove TMS artifact over slightly larger window ===
EEG = pop_tesa_removedata(EEG, [-5 dist_cut], [], {'TMS1'});

% Final visualization
figure; pop_timtopo(EEG, [5200 5999], [], 'ERP data and scalp maps (post artifact removal)');

%% === SSPSIR Step: Muscle Artifact Removal ===
EEG = pop_tesa_sspsir(EEG, 'artScale', 'automatic');  % automatic scaling to remove muscle artifacts

% Visualize cleaned data
EEGLAB_plot_EEG_wrapper(EEG, [5200 5999], 10, 1:size(EEG.data,1));

% Plot averaged signals after SSPSIR
clear avg_data;
avg_data = squeeze(mean(EEG.data,3));
figure; plot(avg_data(:,100:2000)');
figure; pop_timtopo(EEG, [5200 5999], [], 'ERP data and scalp maps (post SSPSIR muscle cleanup)');

%% === Final cleanup: remove extended TMS artifact window ===
EEG = pop_tesa_removedata(EEG, [-5 dist_cut], [], {'TMS1'});