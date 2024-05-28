clear
clc

%% MANUAL
patient = 'Patient010';

%% AUTO for a bit, but define bad elecs and trigger channels, etc
% Reset path for testing
restoredefaultpath

%% Load shared paths
addpath('/Volumes/Research/Epilepsy_ECOG/CodeBase/Matlab');
local_path = strcat('/Users/am4611/Dropbox/Research/ChickenSyntax/data/',patient,'/');

%% Read in data
data_pth = fullfile(local_path,strcat('data/',patient,'_ChickenSyntax_512.EDF'));
if patient == 'Patient001'
    [garbage,data]=edfread(data_pth);
    % Header file for this patient was wrong, must use updated
    load(fullfile(local_path, "data/Patient001_header_updated.mat")); % The header labels were off in the OG file, so they have to be imported separately
else
    [header,data]=edfread(data_pth);
    % Save header
    save(fullfile(local_path,'data',strcat(patient,'_header.mat')), 'header');
end

%% Read in reconstructed trigger log
opts = detectImportOptions(fullfile(local_path,strcat('log/final/',patient,'_log_512_aligned_with_ECoG_with_words_with_latencies.csv')));
word_log = readtable(fullfile(local_path,strcat('log/final/',patient,'_log_512_aligned_with_ECoG_with_words_with_latencies.csv')),opts);
word_log.task = [];

%% plot data
srate = mode(header.frequency);
eegplot(data(1:(find(strcmp([header.label],'EKGL'))-1),:),'srate',srate,'dispchans',60);

%% Define golbals
SJ = patient;
elecs = [1:find(strcmp([header.label],'EKGL'))-1]; 
display(find(strcmp([header.label],'EKGL'))-1) % print number of electrodes

%% Save raw data (pre-referenced) locally
raw_data = data(elecs,:);
%save(fullfile(local_path,'data',strcat(patient,'_raw_preCAR.mat')), 'raw_data', '-v7.3');

%% MANUALLY EDIT
bad_elecs = [33,41,49,64]; 

%% Store channel numbers for trigger, photo, mic
% (See preprocessing script "3 - downsample audio" for which channel is
% which)
trigger_channel = find(strcmp(header.label,'DC2'));
photo_channel = find(strcmp(header.label,'DC1'));
mic_channel = find(strcmp(header.label,'DC3'));

% Plot audio and diode triggers and audio recording
%plot(word_log.onset_trial,'Color','r');
hold on;
plot(data(trigger_channel,:)' .*.01,'Color','b');
%plot(data(photo_channel,:)' .*.04,'Color','r');
plot(data(mic_channel,:)' .*.1 + 10,'Color','c');
plot(word_log.onset_trial*10000 - 5000,'Color','g');
%plot(word_log.onset_chicken*10000 - 10000,'Color','r');
hold off;

%% Save some more stuff
% Get a list of the good electrodes
good_elecs = setdiff(elecs, bad_elecs);
% Data - just brain channels (no EKG, triggers, etc.)
gdat = data((1:elecs(length(elecs))),:);
Labels = header.label(1:elecs(end));
%% Save all data to network folder (TAKES LIKE 10 MINS)
%create_subj_globals(SJ,'ChickenSyntax',srate,srate, elecs, bad_elecs,data_pth);
%get_subj_globals(SJ,'ChickenSyntax');
%trigger = data(trigger_channel,:);
%mic = data(mic_channel,:);
%photo = data(photo_channel,:);  % phototdiode
%save([DTdir dlm 'gdat'],'gdat');
%save([DTdir dlm 'Labels'],'Labels');
%save([DTdir dlm 'trigger'] ,'trigger');
%save([DTdir dlm 'photo'] ,'photo');
%save([DTdir dlm 'mic'] ,'mic');

%% Save bad elecs locally
writecell(Labels(:,bad_elecs)',fullfile(local_path,'data',strcat(patient,'_bad_elecs.csv')));
writematrix(bad_elecs',fullfile(local_path,'data',strcat(patient,'_bad_elecs_order_1indexed.csv')));

%% Common Average Reference (CAR)
%get_subj_globals(SJ,'ChickenSyntax');
%create_CAR(SJ,'ChickenSyntax',elecs(end),elecs(end),[]);

%% Save CAR locally
%save(fullfile(local_path,'data',strcat(patient,'_gdat_CAR.mat')), 'gdat');

%% Adam's Common Average Reference
% Initialize storage - same size as data
car_adam = zeros(size(gdat));
% Loop through electrodes (rows) and subtract each channels' mean
for elec_loop = elecs
    car_adam(elec_loop,:) = gdat(elec_loop,:) - mean(gdat(elec_loop,:));    % remove mean before CAR    
end
% Now get the time-wise average of all GOOD centered channels 
car_adam = sum(car_adam(good_elecs,:),1) ./ length(good_elecs);

% Subtract CAR from all channels
data_car = gdat - car_adam;

% Save
%save(fullfile(local_path,'data',strcat(patient,'_car_data.mat')), 'data_car', '-v7.3');

%% Compare:
% OG data in red
eegplot(gdat,'srate',srate,'dispchans',50,'color',{'r'});
% CAR'ed data in blue
eegplot(data_car,'srate',srate,'dispchans',50,'color',{'b'});

%% Formerly Stage 7, now combined with 6:

%% Get lists of times for each word during character naming
onsetTimes = find((word_log.task_type=="name_character" | ...
    word_log.task_type=="produce_sentence" | ...
    word_log.task_type=="produce_list") & ...
    (word_log.onset_verb == 1 | ...
    word_log.onset_aux == 1 | ...
    word_log.onset_determiner == 1 | ...
    word_log.onset_by == 1 | ...
    word_log.onset_passive_be == 1 | ...
    word_log.onset_chicken==1 | ...
    word_log.onset_dog==1 | ...
    word_log.onset_dracula==1 | ...
    word_log.onset_frankenstein==1 | ...
    word_log.onset_ninja==1 | ...
    word_log.onset_nurse==1));
onsetTimeLabels_word = word_log.word(onsetTimes);
onsetTimeLabels_pos = word_log.pos(onsetTimes);
onsetTimeLabels_dp_or_np = word_log.dp_or_np(onsetTimes);
onsetTimeLabels_syntactic_role = word_log.syntactic_role(onsetTimes);
onsetTimeLabels_production_latency = word_log.samples_from_stimulus_onset(onsetTimes);
onsetTimeLabels_quick_image_duration = word_log.quick_image_duration(onsetTimes);
onsetTimeLabels_noun1 = word_log.noun1(onsetTimes);
onsetTimeLabels_noun2 = word_log.noun2(onsetTimes);
onsetTimeLabels_phase = word_log.phase(onsetTimes);
onsetTimeLabels_verb = word_log.verb(onsetTimes);
onsetTimeLabels_verb_voice = word_log.verb_voice(onsetTimes);
onsetTimes_stimLocked = onsetTimes - str2double(onsetTimeLabels_production_latency);
onsetTimes_substitution_error = word_log.substitution_error(onsetTimes);
onsetTimes_correct_noun = word_log.correct_noun(onsetTimes);
onsetTimes_wrong_noun = word_log.wrong_noun(onsetTimes);
% Added July 2023
onsetTimes_text_active_or_passive = word_log.active_or_passive(onsetTimes);
onsetTimes_arrow = word_log.image_arrow(onsetTimes);
onsetTimes_image_id = word_log.image_id(onsetTimes);

%% Remove trials with NaN onset time values
% Find NaNs in trial onset times and just get rid of those trials
remove_trials = isnan(onsetTimes_stimLocked);
% Shouldn't be more than a handful:
disp(strcat("Removing ",num2str(sum(remove_trials))," trials due to NaN onset times. (Shouldn't be more than a handful.)"))

% Define a double of logical values for which rows to keep (not NaN rows)
keep_trials = logical(1 - remove_trials);

% Subset to just the not NaN rows
onsetTimes = onsetTimes(keep_trials);
onsetTimeLabels_word = onsetTimeLabels_word(keep_trials);
onsetTimeLabels_pos = onsetTimeLabels_pos(keep_trials);
onsetTimeLabels_dp_or_np = onsetTimeLabels_dp_or_np(keep_trials);
onsetTimeLabels_syntactic_role = onsetTimeLabels_syntactic_role(keep_trials);
onsetTimeLabels_noun1 = onsetTimeLabels_noun1(keep_trials);
onsetTimeLabels_noun2 = onsetTimeLabels_noun2(keep_trials);
onsetTimeLabels_verb = onsetTimeLabels_verb(keep_trials);
onsetTimeLabels_verb_voice = onsetTimeLabels_verb_voice(keep_trials);
onsetTimeLabels_quick_image_duration = onsetTimeLabels_quick_image_duration(keep_trials);
onsetTimeLabels_production_latency = onsetTimeLabels_production_latency(keep_trials);
onsetTimeLabels_phase = onsetTimeLabels_phase(keep_trials);
onsetTimes_stimLocked = onsetTimes_stimLocked(keep_trials);
onsetTimes_substitution_error = onsetTimes_substitution_error(keep_trials);
onsetTimes_correct_noun = onsetTimes_correct_noun(keep_trials);
onsetTimes_wrong_noun = onsetTimes_wrong_noun(keep_trials);
% Added July 2023
onsetTimes_text_active_or_passive = onsetTimes_text_active_or_passive(keep_trials);
onsetTimes_arrow = onsetTimes_arrow(keep_trials);
onsetTimes_image_id = onsetTimes_image_id(keep_trials);

%% Turn linguistic information for each trial into a table
trial_labels_table = cell2table([onsetTimeLabels_word, ...
    onsetTimeLabels_syntactic_role, ...
    onsetTimeLabels_dp_or_np, ...
    onsetTimeLabels_pos, ...
    onsetTimeLabels_noun1, ...
    onsetTimeLabels_noun2, ...
    onsetTimeLabels_verb, ...
    onsetTimeLabels_verb_voice, ...
    onsetTimeLabels_quick_image_duration, ...
    onsetTimeLabels_production_latency, ...
    onsetTimes_substitution_error, ...
    onsetTimes_correct_noun, ...
    onsetTimes_wrong_noun, ...
    onsetTimes_text_active_or_passive, ...
    onsetTimes_arrow, ...
    onsetTimes_image_id, ...
    onsetTimeLabels_phase], ...
    'VariableNames',{'word' ...
    'case' ...
    'dp_or_np' ...
    'pos' ...
    'noun1' ...
    'noun2' ...
    'verb' ...
    'verb_voice' ...
    'quick_image_duration' ...
    'production_latency' ...
    'substitution_error' ...
    'correct_noun' ...
    'wrong_noun' ...
    'question_voice' ...
    'arrow' ....
    'image_id' ...
    'phase'});

%% Save
% If the save directory doesn't exist, create it
if ~exist(fullfile(local_path, "data/epoched/"), 'dir')
       mkdir(fullfile(local_path, "data/epoched/"))
end

% Write linguistic data:
writetable(trial_labels_table, fullfile(local_path, strcat("data/epoched/",patient,"_trial_labels_table.csv")));

% How many trials?
nTrials = size(onsetTimes,1);
writematrix(nTrials, fullfile(local_path, strcat("data/epoched/",patient,"_nTrials_before_trial_exclusions.csv")));

%% Filter high-gamma (Hilbert transform)
srate = mode(header.frequency);
%hgdat = abs(my_hilbert(data_car,srate,70,150,1));

%% Average logarithmically-spaced frequency bands within high gamma
n_frequency_bands_hg = 3;
f1_hg = 70;
f2_hg = 150;
freq_bands_hg = logspace(log10(f1_hg),log10(f2_hg),n_frequency_bands_hg+1);
% Exclude 120
freq_bands_hg = [freq_bands_hg(1:2) 118 122 freq_bands_hg(4)];
band_avg_hg = zeros(size(data_car));
for freq_loop = 1:length(freq_bands_hg)-1
    if freq_bands_hg(freq_loop) ~= 118 % Skip the 118 to 122 Hz band
        disp(strcat("Adding HG frequencies from ",string(freq_bands_hg(freq_loop))," to ",string(freq_bands_hg(freq_loop + 1))," Hz."))
        temp = abs(my_hilbert(data_car,srate,...
            freq_bands_hg(freq_loop),freq_bands_hg(freq_loop+1)));
        band_avg_hg = band_avg_hg + temp./mean(temp,2);
    end
end
hgdat = band_avg_hg ./ n_frequency_bands_hg;
disp('multiband - analytic amplitude - high gamma - done!')
clear temp

%% Also get beta, single frequency band since so narrow
n_frequency_bands_beta = 1;
f1_beta = 12;
f2_beta = 30;
freq_bands_beta = logspace(log10(f1_beta),log10(f2_beta),n_frequency_bands_beta+1);
band_avg_beta = zeros(size(data_car));
for freq_loop = 1:length(freq_bands_beta)-1
    temp = abs(my_hilbert(data_car,srate,...
        freq_bands_beta(freq_loop),freq_bands_beta(freq_loop+1)));
    band_avg_beta = band_avg_beta + temp./mean(temp,2);
end
beta_data = band_avg_beta ./ n_frequency_bands_beta;
disp('multiband - analytic amplitude - beta - done!')
clear temp

%% Define bad electrodes to exclude
% Read in old bad elecs
% Many of these were bad due to 60 Hz (line) noise and will be OK now
bad_elecs = table2array(readtable(strcat(local_path,'data/',patient,'_bad_elecs.csv'),'ReadVariableNames',false));
bad_elecs_numeric = table2array(readtable(strcat(local_path,'data/',patient,'_bad_elecs_order_1indexed.csv'),'ReadVariableNames',false));

%% Double check
% In red: original CAR'ed data: should have lots of low-freq stuff
eegplot(data_car,'srate',srate,'dispchans',20,'color',{'r'});
% In blue: HGA, should have only high freq stuff
eegplot(hgdat,'srate',srate,'dispchans',20,'color',{'b'});
% In green: Beta amplitude
eegplot(beta_data,'srate',srate,'dispchans',20,'color',{'g'});

%% Bad in red, good in blue
eegplot(hgdat(bad_elecs_numeric,:),'srate',srate,'dispchans',50,'color',{'r'});
eegplot(hgdat(setdiff(1:size(data_car,1),bad_elecs_numeric),:),'srate',srate,'dispchans',50,'color',{'b'});

%% Remove in black
% Remove electrodes (just define which ones)
% should be a subset of bad_elecs.
% but not all because bad_elecs was defined for CAR.
% some electrodes are noisy and shouldn't be included for CAR but may be
% fine after.
remove_elecs = [bad_elecs_numeric([1,2,4])']; 
eegplot(hgdat(remove_elecs,:),'srate',srate);

%% Automatic - All below 

%% Epoch dimensions
% How much time to epoch?
epochDuration_preOnset = 2200; % in ms
epochDuration_postOnset = 1200; % in ms
nEpochSamples_preOnset = round((epochDuration_preOnset) * srate / 1000);
nEpochSamples_postOnset = round((epochDuration_postOnset) * srate / 1000);
% Flip this for stimulus-locked epochs:
nEpochSamples_preOnset_stimLocked = round((epochDuration_postOnset) * srate / 1000);
nEpochSamples_postOnset_stimLocked = round((epochDuration_preOnset) * srate / 1000);

% Window over which to baseline
baselineDuration = 200; % in ms
nBaselineSamples = round((baselineDuration) * srate / 1000);

% Production-locked trial sample labels
sample_labels = [strcat("sample_neg",flip(string(1:nEpochSamples_preOnset))), "sample_0", strcat("sample_",string(1:nEpochSamples_postOnset))];
% Write trial sample labels
writematrix(sample_labels, fullfile(local_path, strcat("data/epoched/",patient,"_sample_labels.csv")));

% Stimulus-locked trial sample labels
sample_labels_stimLocked = [strcat("sample_neg",flip(string(1:nEpochSamples_preOnset_stimLocked))), "sample_0", strcat("sample_",string(1:nEpochSamples_postOnset_stimLocked))];
% Write trial sample labels
writematrix(sample_labels_stimLocked, fullfile(local_path, strcat("data/epoched/",patient,"_sample_labels_stimLocked.csv")));

% Baseline sample labels
baseline_sample_labels = strcat("sample_neg",flip(string(1:nBaselineSamples)));
% Write baseline sample labels
writematrix(baseline_sample_labels, fullfile(local_path, strcat("data/epoched/",patient,"_baseline_sample_labels.csv")));

%% Set up smoothing
%n_smooth = table2array(readtable(strcat(local_path,'../','n_pre_and_post_data_smoothing_samples.csv'),'ReadVariableNames',false)); % number of pre & post smoothing samples
% (actual number of samples smoothed will be [(2 * n_smooth) + 1])

%% Set up save directories
% Create 1 table for each electrode: nTrials rows and nEpochSamples columns
% Also create a table of labels: (chicken, dog, ...) and POS/case tags
% (NP, DP, ...; subject, object, ...)

% If the save directory doesn't exist, create it
smoothing_loop = "multiband_analytic_amplitude";
% For locked-to-production-onset data:
if ~exist(fullfile(local_path, strcat("data/epoched/electrode_data/",smoothing_loop,"/locked_to_production_onset/z_scored/")), 'dir')
       mkdir(fullfile(local_path, strcat("data/epoched/electrode_data/",smoothing_loop,"/locked_to_production_onset/z_scored/")))
end
if ~exist(fullfile(local_path, strcat("data/epoched/electrode_data/",smoothing_loop,"/locked_to_production_onset/hgp/")), 'dir')
       mkdir(fullfile(local_path, strcat("data/epoched/electrode_data/",smoothing_loop,"/locked_to_production_onset/hgp/")))
end
if ~exist(fullfile(local_path, strcat("data/epoched/electrode_data/",smoothing_loop,"/locked_to_production_onset/percent_change_baseline/")), 'dir')
       mkdir(fullfile(local_path, strcat("data/epoched/electrode_data/",smoothing_loop,"/locked_to_production_onset/percent_change_baseline/")))
end
if ~exist(fullfile(local_path, strcat("data/epoched/electrode_data/",smoothing_loop,"/locked_to_production_onset/raw_potentials_postCAR/")), 'dir')
       mkdir(fullfile(local_path, strcat("data/epoched/electrode_data/",smoothing_loop,"/locked_to_production_onset/raw_potentials_postCAR/")))
end
if ~exist(fullfile(local_path, strcat("data/epoched/electrode_data/",smoothing_loop,"/locked_to_production_onset/z_scored_beta/")), 'dir')
       mkdir(fullfile(local_path, strcat("data/epoched/electrode_data/",smoothing_loop,"/locked_to_production_onset/z_scored_beta/")))
end

% And for locked-to-stimulus-onset-data:
if ~exist(fullfile(local_path, strcat("data/epoched/electrode_data/",smoothing_loop,"/locked_to_stimulus_onset/z_scored/")), 'dir')
       mkdir(fullfile(local_path, strcat("data/epoched/electrode_data/",smoothing_loop,"/locked_to_stimulus_onset/z_scored/")))
end
if ~exist(fullfile(local_path, strcat("data/epoched/electrode_data/",smoothing_loop,"/locked_to_stimulus_onset/hgp/")), 'dir')
       mkdir(fullfile(local_path, strcat("data/epoched/electrode_data/",smoothing_loop,"/locked_to_stimulus_onset/hgp/")))
end
if ~exist(fullfile(local_path, strcat("data/epoched/electrode_data/",smoothing_loop,"/locked_to_stimulus_onset/percent_change_baseline/")), 'dir')
       mkdir(fullfile(local_path, strcat("data/epoched/electrode_data/",smoothing_loop,"/locked_to_stimulus_onset/percent_change_baseline/")))
end
if ~exist(fullfile(local_path, strcat("data/epoched/electrode_data/",smoothing_loop,"/locked_to_stimulus_onset/raw_potentials_postCAR/")), 'dir')
       mkdir(fullfile(local_path, strcat("data/epoched/electrode_data/",smoothing_loop,"/locked_to_stimulus_onset/raw_potentials_postCAR/")))
end
if ~exist(fullfile(local_path, strcat("data/epoched/electrode_data/",smoothing_loop,"/locked_to_stimulus_onset/z_scored_beta/")), 'dir')
       mkdir(fullfile(local_path, strcat("data/epoched/electrode_data/",smoothing_loop,"/locked_to_stimulus_onset/z_scored_beta/")))
end

% And for baseline
if ~exist(fullfile(local_path, strcat("data/epoched/electrode_data/",smoothing_loop,"/baseline_",int2str(baselineDuration),"ms_pre_stimulus/")), 'dir')
       mkdir(fullfile(local_path, strcat("data/epoched/electrode_data/",smoothing_loop,"/baseline_",int2str(baselineDuration),"ms_pre_stimulus/")))
end

%% Epoch: Make electrode matrices
smoothing_loop = "multiband_analytic_amplitude";
disp(strcat("Beginning ",smoothing_loop," loop."))

% Include bad electrodes, for now anyway
if length(remove_elecs)==0
    remove_elecs=["no_bad_elecs"]
end

% Open 7 parallel workers
parpool(4)

% Loop thru electrodes in parallel
parfor electrode = 1:length(Labels) % Cycle thru electrodes
    % electrode = 1; % (for troubleshooting)

    % Only for electrodes not listed as "bad"
    if any(remove_elecs == electrode)
        disp(strcat("Excluding bad electrode ",string(Labels(electrode))))

    else
        disp(strcat("Electrode ",string(electrode)," of ",string(length(Labels))))
        
        % Set up matrices for labels (chicken, dog, ...) and prod-locked ECoG data
        %elec_data = nan(nTrials, nEpochSamples_preOnset + 1 + nEpochSamples_postOnset);
        elec_data_z = nan(nTrials, nEpochSamples_preOnset + 1 + nEpochSamples_postOnset);
        elec_data_z_beta = nan(nTrials, nEpochSamples_preOnset + 1 + nEpochSamples_postOnset);
        %elec_data_percent = nan(nTrials, nEpochSamples_preOnset + 1 + nEpochSamples_postOnset);
        elec_data_potential = nan(nTrials, nEpochSamples_preOnset + 1 + nEpochSamples_postOnset);

        % Set up matrices for labels (chicken, dog, ...) and stim-locked ECoG data
        %elec_data_stimLocked = nan(nTrials, nEpochSamples_preOnset_stimLocked + 1 + nEpochSamples_postOnset_stimLocked);
        elec_data_z_stimLocked = nan(nTrials, nEpochSamples_preOnset_stimLocked + 1 + nEpochSamples_postOnset_stimLocked);
        elec_data_z_stimLocked_beta = nan(nTrials, nEpochSamples_preOnset_stimLocked + 1 + nEpochSamples_postOnset_stimLocked);
        %elec_data_percent_stimLocked = nan(nTrials, nEpochSamples_preOnset_stimLocked + 1 + nEpochSamples_postOnset_stimLocked);
        elec_data_potential_stimLocked = nan(nTrials, nEpochSamples_preOnset_stimLocked + 1 + nEpochSamples_postOnset_stimLocked);

        % Set up matrix for baseline data
        baseline_data = nan(nTrials, nBaselineSamples);
        baseline_data_beta = nan(nTrials, nBaselineSamples);

        % Define slices of data outside for loop so that parfor will work
        current_hgdat = hgdat(electrode, :)
        current_data_car = data_car(electrode, :)
        current_beta_data = beta_data(electrode, :)

        % Add one row at a time of high gamma power
        for trial_loop = 1:nTrials
            % trial_loop = 1; % (for troubleshooting)
            % Get onset indices
            temp_prodLocked_onset = onsetTimes(trial_loop);
            temp_stimLocked_onset = onsetTimes_stimLocked(trial_loop);

            % Find the 200ms window pre-stimulus onset for baseline
            temp_baseline_end = temp_prodLocked_onset - str2num(cell2mat(onsetTimeLabels_production_latency(trial_loop))) - 1;
            temp_baseline_start = temp_baseline_end - nBaselineSamples + 1;

            % Find start and end indices for production-locked epochs
            temp_trial_start = temp_prodLocked_onset - nEpochSamples_preOnset;
            temp_trial_end = temp_prodLocked_onset + nEpochSamples_postOnset;

            % Find start and end indices for stimulus-locked epochs
            temp_trial_start_stimLocked= temp_stimLocked_onset - nEpochSamples_preOnset_stimLocked;
            temp_trial_end_stimLocked = temp_stimLocked_onset + nEpochSamples_postOnset_stimLocked;

            % Get data
            % Get this trial's baseline data and stats
            baseline = current_hgdat(:, temp_baseline_start:temp_baseline_end);
            baseline_mean = mean(baseline);
            baseline_sd = std(baseline);
            % Get this trial's prod-locked data
            trial = current_hgdat(:, temp_trial_start:temp_trial_end);
            trial_potential = current_data_car(:, temp_trial_start:temp_trial_end);
            % Get this trial's stim-locked data
            trial_stimLocked = current_hgdat(:, temp_trial_start_stimLocked:temp_trial_end_stimLocked);
            trial_potential_stimLocked = current_data_car(:, temp_trial_start_stimLocked:temp_trial_end_stimLocked);

            % Add to matrices
            baseline_data(trial_loop,:) = baseline;
            % Prod locked
            %elec_data(trial_loop,:) = trial;
            elec_data_z(trial_loop,:) = (trial - baseline_mean) / baseline_sd;
            %elec_data_percent(trial_loop,:) = (trial / baseline_mean) - 1;
            elec_data_potential(trial_loop,:) = trial_potential;
            % Stim locked
            %elec_data_stimLocked(trial_loop,:) = trial_stimLocked;
            elec_data_z_stimLocked(trial_loop,:) = (trial_stimLocked - baseline_mean) / baseline_sd;
            %elec_data_percent_stimLocked(trial_loop,:) = (trial_stimLocked / baseline_mean) - 1;
            elec_data_potential_stimLocked(trial_loop,:) = trial_potential_stimLocked;

            % Get beta data
            % Get this trial's baseline data and stats
            baseline_beta = current_beta_data(:, temp_baseline_start:temp_baseline_end);
            baseline_mean_beta = mean(baseline_beta);
            baseline_sd_beta = std(baseline_beta);
            % Get this trial's prod-locked data
            trial_beta = current_beta_data(:, temp_trial_start:temp_trial_end);
            % Get this trial's stim-locked data
            trial_stimLocked_beta = current_beta_data(:, temp_trial_start_stimLocked:temp_trial_end_stimLocked);

            % Add to matrices
            baseline_data_beta(trial_loop,:) = baseline_beta;
            % Prod locked
            elec_data_z_beta(trial_loop,:) = (trial_beta - baseline_mean_beta) / baseline_sd_beta;
            % Stim locked
            elec_data_z_stimLocked_beta(trial_loop,:) = (trial_stimLocked_beta - baseline_mean_beta) / baseline_sd_beta;

        end
        % Write CSV of potentials at 512 Hz
        % Production-locked data:
        %writematrix(elec_data, fullfile(local_path, strcat("data/epoched/electrode_data/",smoothing_loop,"/locked_to_production_onset/hgp/"), strcat(patient,"_",string(Labels(electrode)),".csv")))
        writematrix(elec_data_z, fullfile(local_path, strcat("data/epoched/electrode_data/",smoothing_loop,"/locked_to_production_onset/z_scored/"), strcat(patient,"_",string(Labels(electrode)),".csv")))
        writematrix(elec_data_z_beta, fullfile(local_path, strcat("data/epoched/electrode_data/",smoothing_loop,"/locked_to_production_onset/z_scored_beta/"), strcat(patient,"_",string(Labels(electrode)),".csv")))
        %writematrix(elec_data_percent, fullfile(local_path, strcat("data/epoched/electrode_data/",smoothing_loop,"/locked_to_production_onset/percent_change_baseline/"), strcat(patient,"_",string(Labels(electrode)),".csv")))
        writematrix(elec_data_potential, fullfile(local_path, strcat("data/epoched/electrode_data/",smoothing_loop,"/locked_to_production_onset/raw_potentials_postCAR/"), strcat(patient,"_",string(Labels(electrode)),".csv")))
        % Stimulus-locked data:
        %writematrix(elec_data_stimLocked, fullfile(local_path, strcat("data/epoched/electrode_data/",smoothing_loop,"/locked_to_stimulus_onset/hgp/"), strcat(patient,"_",string(Labels(electrode)),".csv")))
        writematrix(elec_data_z_stimLocked, fullfile(local_path, strcat("data/epoched/electrode_data/",smoothing_loop,"/locked_to_stimulus_onset/z_scored/"), strcat(patient,"_",string(Labels(electrode)),".csv")))
        writematrix(elec_data_z_stimLocked_beta, fullfile(local_path, strcat("data/epoched/electrode_data/",smoothing_loop,"/locked_to_stimulus_onset/z_scored_beta/"), strcat(patient,"_",string(Labels(electrode)),".csv")))
        %writematrix(elec_data_percent_stimLocked, fullfile(local_path, strcat("data/epoched/electrode_data/",smoothing_loop,"/locked_to_stimulus_onset/percent_change_baseline/"), strcat(patient,"_",string(Labels(electrode)),".csv")))
        writematrix(elec_data_potential_stimLocked, fullfile(local_path, strcat("data/epoched/electrode_data/",smoothing_loop,"/locked_to_stimulus_onset/raw_potentials_postCAR/"), strcat(patient,"_",string(Labels(electrode)),".csv")))
        % Baselines:
        %writematrix(baseline_data, fullfile(local_path, strcat("data/epoched/electrode_data/",smoothing_loop,"/baseline_",int2str(baselineDuration),"ms_pre_stimulus/"), strcat(patient,"_",string(Labels(electrode)),".csv")))
    end
end

% Close parallel workers
delete(gcp('nocreate'))

%% Finish up
disp("Done")
