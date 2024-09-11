%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:  mass-univariate analysis with cluster-based permutation tests   %
%         of event-related potentials                                     %
% Author: Mareike J. Hülsemann                                            %
% Date:   06/2020                                                         %
% Note:                                                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BSD 3-Clause License
% 
% Copyright (c) 2020, Mareike Huelsemann
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
% 
% 3. Neither the name of the copyright holder nor the names of its
%    contributors may be used to endorse or promote products derived from
%    this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


%% Load Raw Data and Calculate Mean EEG for t-Tests

% initialize random number generater (otherwise the same random numbers are generated with each start of MATLAB)
rng(datenum(datetime));

% start EEGLAB
eeglab_path = '/your/path/to/eeglab/'; % EEGLAB toolbox v2019.1 was used in original analysis
addpath(eeglab_path); 
eeglab; 
pop_editoptions('option_single', 0);


PATH_IN       = '/your/path/to/preprocessed/data/';
allSets       = dir([PATH_IN '*.set']);
PATH_StimList = '/path/to/trial/list/.../rejincl_stimulus_list/';

% define path for secure-saving data
PATH_OUT = '/path/where/to/save/data/';
note     = 'same amount of trials';

load('/path/to/EEG.chanlocs/structure/')
AMOUNT_CHANS   = length(EEG_chanlocs);
AMOUNT_SAMPLES = 2000;

MEANS   = {'corr_cgrt_vol1', 'corr_cgrt_vol2', 'corr_cgrt_vol3', 'corr_cgrt_vol4', ...
    'corr_icgr_vol1', 'corr_icgr_vol2', 'corr_icgr_vol3', 'corr_icgr_vol4'};

% initialise variables
chan_order_control = cell(length(allSets),AMOUNT_CHANS,2);
nperm              = 1000;
trial_amount       = nan(length(allSets),numel(MEANS)/2,numel(MEANS)/4);

% initialize MeanActivty Data for t-tests
idx_elec = 8:61; % electrode indexes to use for statistical testing
stimes = 400:100:1700;
etimes = 499:100:1799;  % originally with etimes = [500:100:1700 1799];
MeanActivity = nan(length(allSets),numel(idx_elec),numel(stimes),numel(MEANS));

% loop over all subjects to calculate original EEG mean activity
for n = 1:length(allSets)
    
    disp([allSets(n).name(1:9) ' started: ' datestr(datetime)])
    
    % load data set, which is preprocessed 
    EEG = pop_loadset('filename',allSets(n).name, 'filepath', PATH_IN);
    
    % adjust channel order
    new_order = nan(length(EEG_chanlocs), 1);
    for ci=1:length(EEG_chanlocs)
        new_order(ci) = find(strcmp({EEG.chanlocs.labels}, EEG_chanlocs(ci).labels));
    end
    EEG.data = EEG.data(new_order,:,:);
    chan_order_control(n,:,1) = {EEG.chanlocs.labels}';
    chan_order_control(n,:,2) = {EEG.chanlocs(new_order).labels}';
    % % % % validity check
    % % % a(:,1) = {EEG_chanlocs.labels}';
    % % % a(:,2) = {EEG.chanlocs(new_order).labels}';
    % % % a(:,3) = {EEG.chanlocs.labels}';
    % % % a(:,4) = num2cell(EEG.data(:,1,1));
    % % % a(:,5) = num2cell(EEG.data(new_order,1,1));

    % load stimulus list
    load([PATH_StimList allSets(n).name(1:9) '_StimulusList.mat'])
    stim_list = stim_list(2:end,:);
    
    % define trial indices
    %%%position without noise trials
    idx_pos_L = find(strcmp(stim_list(:, 9),'Liegen') & ~strcmp(stim_list(:,10),'Rauschen'));
    idx_pos_S = find(strcmp(stim_list(:, 9),'Stehen') & ~strcmp(stim_list(:,10),'Rauschen'));
    %%%congruency
    idx_cgr_0 = find([stim_list{:,11}]'==0);
    idx_cgr_1 = find([stim_list{:,11}]'==1);
    
    %%%volume
    idx_vol_1 = find([stim_list{:,12}]'==1 & ~strcmp(stim_list(:,10),'Rauschen'));
    idx_vol_2 = find([stim_list{:,12}]'==2 & ~strcmp(stim_list(:,10),'Rauschen'));
    idx_vol_3 = find([stim_list{:,12}]'==3 & ~strcmp(stim_list(:,10),'Rauschen'));
    idx_vol_4 = find([stim_list{:,12}]'==4 & ~strcmp(stim_list(:,10),'Rauschen'));
    
    %%%stimulus amount 
    % 20 different stimuli, 4 volumes, 6 repetitions (in 2 positions) = 480 trials
    % 18 different stimuli, 4 volumes, 5 repetitions (in 2 positions) = 360 trials
    idx_stima = find([stim_list{:, 8}]'>4 & ~strcmp(stim_list(:,10),'Rauschen'));
    %%%artefact free trial
    idx_noart = find([stim_list{:,18}]'==1);
    
    %%%correct response
    idx_ACC_1 = find([stim_list{:,14}]'==1);
    %%%no response given
    idx_norep = find(isnan([stim_list{:,13}])');
    
    % definition of averaged trials
    trials{1,1} = intersect(intersect(intersect(intersect(idx_cgr_1, idx_vol_1), idx_stima), idx_noart), idx_ACC_1);
    trials{1,2} = intersect(intersect(intersect(intersect(idx_cgr_1, idx_vol_2), idx_stima), idx_noart), idx_ACC_1);
    trials{1,3} = intersect(intersect(intersect(intersect(idx_cgr_1, idx_vol_3), idx_stima), idx_noart), idx_ACC_1);
    trials{1,4} = intersect(intersect(intersect(intersect(idx_cgr_1, idx_vol_4), idx_stima), idx_noart), idx_ACC_1);
    trials{2,1} = intersect(intersect(intersect(intersect(idx_cgr_0, idx_vol_1), idx_stima), idx_noart), idx_ACC_1);
    trials{2,2} = intersect(intersect(intersect(intersect(idx_cgr_0, idx_vol_2), idx_stima), idx_noart), idx_ACC_1);
    trials{2,3} = intersect(intersect(intersect(intersect(idx_cgr_0, idx_vol_3), idx_stima), idx_noart), idx_ACC_1);
    trials{2,4} = intersect(intersect(intersect(intersect(idx_cgr_0, idx_vol_4), idx_stima), idx_noart), idx_ACC_1);
    
    % choose random trials for condition with larger amount of trials
    % find condition with less amount of trials (a-trials)
    for vi=1:4 % 4 volume levels
        [min_trials,index] = min([numel(trials{1,vi}) numel(trials{2,vi})]);
        if index==1 % condition 1 contains less trials
            % define a-trials of condition 2
            rp = randperm(numel(trials{2,vi}),numel(trials{1,vi}));
            % only use these for later averaging
            trials{2,vi} = trials{2,vi}(rp);
        else % condition 2 contains less trials
            % define a-trials of condition 1
            rp = randperm(numel(trials{1,vi}),numel(trials{2,vi}));
            % only use these for later averaging
            trials{1,vi} = trials{1,vi}(rp);
        end
    end
    
    % save amount trials
    trial_amount(n,1,:) = [numel(trials{1,1}) numel(trials{2,1})];
    trial_amount(n,2,:) = [numel(trials{1,2}) numel(trials{2,2})];
    trial_amount(n,3,:) = [numel(trials{1,3}) numel(trials{2,3})];
    trial_amount(n,4,:) = [numel(trials{1,4}) numel(trials{2,4})];
    
    % save single-trial data per condition
    single_trial_cgrt_vol1 = EEG.data(:,:,trials{1,1}); save([PATH_OUT allSets(n).name(1:9) '_cgrt_vol1.mat'],'single_trial_cgrt_vol1');
    single_trial_cgrt_vol2 = EEG.data(:,:,trials{1,2}); save([PATH_OUT allSets(n).name(1:9) '_cgrt_vol2.mat'],'single_trial_cgrt_vol2');
    single_trial_cgrt_vol3 = EEG.data(:,:,trials{1,3}); save([PATH_OUT allSets(n).name(1:9) '_cgrt_vol3.mat'],'single_trial_cgrt_vol3');
    single_trial_cgrt_vol4 = EEG.data(:,:,trials{1,4}); save([PATH_OUT allSets(n).name(1:9) '_cgrt_vol4.mat'],'single_trial_cgrt_vol4');
    single_trial_icgr_vol1 = EEG.data(:,:,trials{2,1}); save([PATH_OUT allSets(n).name(1:9) '_icgr_vol1.mat'],'single_trial_icgr_vol1');
    single_trial_icgr_vol2 = EEG.data(:,:,trials{2,2}); save([PATH_OUT allSets(n).name(1:9) '_icgr_vol2.mat'],'single_trial_icgr_vol2');
    single_trial_icgr_vol3 = EEG.data(:,:,trials{2,3}); save([PATH_OUT allSets(n).name(1:9) '_icgr_vol3.mat'],'single_trial_icgr_vol3');
    single_trial_icgr_vol4 = EEG.data(:,:,trials{2,4}); save([PATH_OUT allSets(n).name(1:9) '_icgr_vol4.mat'],'single_trial_icgr_vol4');
    
    % calculate mean for t-Tests: first 4 stimuli excluded, correct response trials, CGR
    for i=1:numel(stimes)
        col=0;    
        col=col+1; MeanActivity(n,:,i,col) = nanmean(nanmean(EEG.data(idx_elec,find(EEG.times ==  stimes(i)):find(EEG.times ==  etimes(i)),trials{1,1}),3),2);
        col=col+1; MeanActivity(n,:,i,col) = nanmean(nanmean(EEG.data(idx_elec,find(EEG.times ==  stimes(i)):find(EEG.times ==  etimes(i)),trials{1,2}),3),2);
        col=col+1; MeanActivity(n,:,i,col) = nanmean(nanmean(EEG.data(idx_elec,find(EEG.times ==  stimes(i)):find(EEG.times ==  etimes(i)),trials{1,3}),3),2);
        col=col+1; MeanActivity(n,:,i,col) = nanmean(nanmean(EEG.data(idx_elec,find(EEG.times ==  stimes(i)):find(EEG.times ==  etimes(i)),trials{1,4}),3),2);
        col=col+1; MeanActivity(n,:,i,col) = nanmean(nanmean(EEG.data(idx_elec,find(EEG.times ==  stimes(i)):find(EEG.times ==  etimes(i)),trials{2,1}),3),2);
        col=col+1; MeanActivity(n,:,i,col) = nanmean(nanmean(EEG.data(idx_elec,find(EEG.times ==  stimes(i)):find(EEG.times ==  etimes(i)),trials{2,2}),3),2);
        col=col+1; MeanActivity(n,:,i,col) = nanmean(nanmean(EEG.data(idx_elec,find(EEG.times ==  stimes(i)):find(EEG.times ==  etimes(i)),trials{2,3}),3),2);
        col=col+1; MeanActivity(n,:,i,col) = nanmean(nanmean(EEG.data(idx_elec,find(EEG.times ==  stimes(i)):find(EEG.times ==  etimes(i)),trials{2,4}),3),2);
    end
end

% variables to secure save at the end:
save([PATH_OUT allSets(n).name(1:6) 'chan_order_control.mat'],'chan_order_control');
save([PATH_OUT allSets(n).name(1:6) 'CGR_trial_amount.mat'],'trial_amount');
save([PATH_OUT allSets(n).name(1:6) 'CGR_MeanActivity.mat'],'MeanActivity', 'trial_amount', 'note');


%% Execute t-Tests for original data

TTEST = nan(numel(idx_elec),numel(stimes),4,3); % 4 for 4 volume levels, 2 for two test statistics T, p, omega)
nobs = size(allSets,1)*2;
for ti=1:numel(stimes)
    for ei=1:numel(idx_elec)
        for vi=1:4 % volume levels
            
            [H,P,CI,STATS] = ttest(squeeze(MeanActivity(:,ei,ti,vi)), squeeze(MeanActivity(:,ei,ti,4+vi)));
            
            % calculate effect size omega^2
            Tval = STATS.tstat;
            fs=((Tval^2)-1)/(nobs);
            wval = (fs)/(1+fs);
            
            TTEST(ei,ti,vi,1) = Tval; % save Temp
            TTEST(ei,ti,vi,2) = P;    % save p-val
            TTEST(ei,ti,vi,3) = wval; % save omega^2           
            
            
        end
    end
end

% variables to secure save at the end:
save([PATH_OUT allSets(n).name(1:6) 'CGR_TTEST.mat'],'TTEST', 'trial_amount', 'note');


