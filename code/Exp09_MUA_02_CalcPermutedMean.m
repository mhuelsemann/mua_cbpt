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


%% Load Raw Data and Calculate Permuted Mean EEG for t-Tests

% initialize random number generater (otherwise the same random numbers are generated with each start of MATLAB)
rng('shuffle');

% define DATA PATH
VOLUMES  = {'vol1','vol2','vol3','vol4'};
PATH_IN  = '/your/path/to/preprocessed/data/';
% note     = 'different amount of trials';
note     = 'same amount of trials';
allSets  = dir([PATH_IN 'Exp09_*_cgrt_' VOLUMES{1} '.mat']);

% initialise variables
idx_elec             = 8:61;
stimes               = 400:100:1700;
etimes               = 499:100:1799;
nperm                = 10000; % run programm with less permuations multiple times and combine matrices later on
MeanActivity_permuted = nan(size(allSets,1),numel(idx_elec),numel(stimes),2*numel(VOLUMES),nperm);
EEG_times = -200:1799;



% loop over all subjects
for n=1:size(allSets,1)
    disp(['Started VP' allSets(n).name(7:9) ': ' datestr(datetime)])
    % loop over all volumes
    for vi=1:4
        VOLUME = VOLUMES{vi};
        
        % load data
        single_trial_cgrt = importdata([PATH_IN allSets(n).name(1:9) '_cgrt_' VOLUME '.mat']);
        single_trial_icgr = importdata([PATH_IN allSets(n).name(1:9) '_icgr_' VOLUME '.mat']);
        
        trials_cgrt = size(single_trial_cgrt,3);
        trials_icgr = size(single_trial_icgr,3);
        trials_sum  = trials_cgrt + trials_icgr;
        
        single_trial(:,:, 1             : trials_cgrt) = single_trial_cgrt;
        single_trial(:,:, trials_cgrt+1 : trials_sum ) = single_trial_icgr;
        
        
        % permutation of trials (conserving the amount of trials per condition)
        for pi=1:nperm
            % generate random trial order
            rp = randperm(trials_sum);
            
            % take first n(trials_cgrt) for permuted condition 1
            perm_cgrt = single_trial(:,:,rp(1             : trials_cgrt));
            % take the last n(trials_icgr) for permuted condition 2
            perm_icgr = single_trial(:,:,rp(trials_cgrt+1 : trials_sum ));
            
            % calculate perm_mean for t-Tests: first 4 stimuli excluded, correct response trials, CGR
            for ti=1:numel(stimes)
                MeanActivity_permuted(n,:,ti,  vi,pi) = nanmean(nanmean(perm_cgrt(idx_elec,find(EEG_times ==  stimes(ti)):find(EEG_times ==  etimes(ti)),:),3),2);
                MeanActivity_permuted(n,:,ti,4+vi,pi) = nanmean(nanmean(perm_icgr(idx_elec,find(EEG_times ==  stimes(ti)):find(EEG_times ==  etimes(ti)),:),3),2);
            end
        end
    end
    
    
end

%%%% DO NOT FORGET TO SAVE! %%%%