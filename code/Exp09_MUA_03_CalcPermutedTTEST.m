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


%% Execute t-Tests for permuted data

% define settings
% note = 'different_amount_trials';
note = 'same_amount_trials';
suffix = '_vol1'; % for FILE_IN and _OUT

% load permuted MeanActivty data
PATH_IN = ['/your/path/to/preprocessed/data/' note '/'];
load([PATH_IN 'Exp09_CGR_MeanActivity_permuted' suffix '.mat'])
load([PATH_IN 'Exp09_CGR_trial_amount.mat'])
PATH_OUT = PATH_IN;
FILE_OUT = ['Exp09_CGR_TTEST_permuted' suffix '.mat'];

% initialize
nperm          = size(MeanActivity_permuted,5);
idx_elec       = 8:61;
stimes         = 400:100:1700;
etimes         = 499:100:1799;
TTEST_permuted = nan(nperm, numel(idx_elec), numel(stimes), 4, 3); % 4 for 4 volume levels, 2 for two test statistics T, p, omega)
nobs           = size(MeanActivity_permuted,1)*2;

for pi=1:nperm
    disp(num2str(pi))
    for ti=1:numel(stimes)
        for ei=1:numel(idx_elec)
            for vi=1:1%4 % volume levels
                
                [H,P,CI,STATS] = ttest(squeeze(MeanActivity_permuted(:,ei,ti,vi,pi)), squeeze(MeanActivity_permuted(:,ei,ti,4+vi,pi)));
                
                % calculate effect site omega^2
                Tval = STATS.tstat;
                fs=((Tval^2)-1)/(nobs);
                wval = (fs)/(1+fs);
                
                TTEST_permuted(pi,ei,ti,vi,1) = Tval; % save Temp
                TTEST_permuted(pi,ei,ti,vi,2) = P;    % save p-val
                TTEST_permuted(pi,ei,ti,vi,3) = wval; % save omega^2
                
                
            end
        end
    end
end

% variables to secure save at the end:
save([PATH_OUT  FILE_OUT],'TTEST_permuted', 'trial_amount', 'note');


