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


%% Definition of Embodiment Cluster

% load empirical t-vals
PATH_IN = '/your/path/to/preprocessed/data/';
load([PATH_IN 'Exp09_CGR_TTEST.mat'])
tvals = squeeze(TTEST(:,:,:,1));

% define critical t-val
tails = 'twosided';
alpha = .05;
df    = 65;
if strcmp(tails,'twosided'), alpha = alpha/2; end
tcrit = tinv(1-alpha,df); % critical t-value

% load chanloc structure
load('/path/to/EEG.chanlocs/structure/', 'EEG_chanlocs54');
chanlocs = EEG_chanlocs54;
% choose max. neighbour distance:
% 40 corresponds to 4 direct neighbours per index electrode (above, below, left, right)
%                   0 1 0
%                   1 X 1
%                   0 1 0
% 55 corresponds to 8 direct neighbours per index electrode (above, below, left, right, + diagonal neighbours)
%                   1 1 1
%                   1 X 1
%                   1 1 1
max_elec_dist = 40;
if     max_elec_dist==40, note2 = 'calculated with 4 neighbours (dist=40)';
elseif max_elec_dist==55, note2 = 'calculated with 8 neighbours (dist=55)';
else,                     note2 = 'neighbour distance neither 40 nor 55';
end

for vi=1:size(tvals,3) % 4 volume levels
    
    tmp_tvals = squeeze(tvals(:,:,vi));
    
    %% define (positive) significant t-vals
    tvals_psig = zeros(size(tmp_tvals));
    tvals_psig(tmp_tvals > tcrit) = 1;
    
    
    %% define (positive) clusters within each electrode
    elec_pclusters = cell(1,numel(chanlocs));
    for ei=1:numel(chanlocs)
        ci=0;
        for ti=1:size(tvals_psig,2)
            if ti==1 % begin of array
                if tvals_psig(ei,ti) == 1
                    % begin of a cluster
                    ci=ci+1;
                    clusterstart = ti;
                else % tvals_sig(ti) == 0
                    % do nothing
                end
            else % ti!=1
                if tvals_psig(ei,ti) == 1 && tvals_psig(ei,ti-1) == 0
                    % begin of a cluster
                    ci=ci+1;
                    clusterstart = ti;
                elseif tvals_psig(ei,ti) == 1 && tvals_psig(ei,ti-1) == 1
                    % continuation of a cluster
                    % do nothing
                elseif tvals_psig(ei,ti) == 0 && tvals_psig(ei,ti-1) == 1
                    % end of cluster one time point before
                    clusterend = ti-1;
                    elec_pclusters{ei}{ci} = clusterstart:clusterend;
                else % tvals_sig(ti) == 0 && tvals_sig(ti-1) == 0
                    % do nothing
                end
            end
            if ti==size(tvals_psig,2) && tvals_psig(ei,ti)==1 % end of array and cluster end
                clusterend = ti;
                elec_pclusters{ei}{ci} = clusterstart:clusterend;
            end
        end
    end
    clear clusterstart clusterend ci ei ti
    
    
    %% define (positive) clusters across electrodes
    tvals_psig(tvals_psig==1) = -1;
    clusternr = 0;
    for ei=1:numel(chanlocs) % loop over all electrodes
        for ci=1:numel(elec_pclusters{ei}) % loop over all elec_clusters of the current electrode
            
            % Is cluster already named?
            if tvals_psig(ei,elec_pclusters{ei}{ci}(1)) == -1 % no, not yet named
                clusternr = clusternr +1; % increase by one for new cluster
                tvals_psig(ei,elec_pclusters{ei}{ci}) = clusternr;
                %imagesc(tvals_psig(elec_order,:)); colorbar; set(gca,'YTick',1:54,'YTickLabels',{EEG_chanlocs54(elec_order).labels}); drawnow; pause(1);
            end
            
            % Are there neighbours with overlapping clusters in time?
            for ni=1:numel(chanlocs) % loop over all other electrodes
                % find potential neighbours
                if ei==ni % do nothing, as this is the same electrode
                    % verify, whether electrode is direct neighbour
                elseif sqrt((chanlocs(ei).X - chanlocs(ni).X)^2 + ...
                        (chanlocs(ei).Y - chanlocs(ni).Y)^2 + ...
                        (chanlocs(ei).Z - chanlocs(ni).Z)^2) < max_elec_dist
                    
                    for nci=1:numel(elec_pclusters{ni}) % loop over all clusters of the neighbour
                        
                        % Does cluster overlap with the index cluster?
                        if intersect(elec_pclusters{ei}{ci},elec_pclusters{ni}{nci}) % yes, overlapping
                            % Is this neighbour-cluster already named?
                            if tvals_psig(ni,elec_pclusters{ni}{nci}(1)) == -1 % no, not yet named
                                tvals_psig(ni,elec_pclusters{ni}{nci}) = tvals_psig(ei,elec_pclusters{ei}{ci}(1)); % name it according to index cluster
                                %imagesc(tvals_psig(elec_order,:)); colorbar; set(gca,'YTick',1:54,'YTickLabels',{EEG_chanlocs54(elec_order).labels}); drawnow; pause(1);
                            else % yes, already named. Is this name different from the index cluster?
                                if tvals_psig(ei,elec_pclusters{ei}{ci}(1)) ~= tvals_psig(ni,elec_pclusters{ni}{nci}(1)) % if cluster number is different
                                    lower_cluster_nr  = min(tvals_psig(ei,elec_pclusters{ei}{ci}(1)), tvals_psig(ni,elec_pclusters{ni}{nci}(1)));
                                    higher_cluster_nr = max(tvals_psig(ei,elec_pclusters{ei}{ci}(1)), tvals_psig(ni,elec_pclusters{ni}{nci}(1)));
                                    tvals_psig(tvals_psig==higher_cluster_nr) = lower_cluster_nr;
                                    %imagesc(tvals_psig(elec_order,:)); colorbar; set(gca,'YTick',1:54,'YTickLabels',{EEG_chanlocs54(elec_order).labels}); drawnow; pause(1);
                                end
                            end
                        else % no, clusters do not overlap, do nothing
                        end
                    end
                else % electrode is not a neighbour, do nothing
                end
            end
        end
    end
    clear ei ci ni nci elec_pclusters
    
    
    %% define (negative) significant t-vals
    tvals_nsig = zeros(size(tmp_tvals));
    tvals_nsig(tmp_tvals < (-tcrit)) = 1;
    
    
    %% define (negative) clusters within each electrode
    elec_nclusters = cell(1,numel(chanlocs));
    for ei=1:numel(chanlocs)
        ci=0;
        for ti=1:size(tvals_nsig,2)
            if ti==1 % begin of array
                if tvals_nsig(ei,ti) == 1
                    % begin of a cluster
                    ci=ci+1;
                    clusterstart = ti;
                else % tvals_sig(ti) == 0
                    % do nothing
                end
            else % ti!=1
                if tvals_nsig(ei,ti) == 1 && tvals_nsig(ei,ti-1) == 0
                    % begin of a cluster
                    ci=ci+1;
                    clusterstart = ti;
                elseif tvals_nsig(ei,ti) == 1 && tvals_nsig(ei,ti-1) == 1
                    % continuation of a cluster
                    % do nothing
                elseif tvals_nsig(ei,ti) == 0 && tvals_nsig(ei,ti-1) == 1
                    % end of cluster one time point before
                    clusterend = ti-1;
                    elec_nclusters{ei}{ci} = clusterstart:clusterend;
                else % tvals_sig(ti) == 0 && tvals_sig(ti-1) == 0
                    % do nothing
                end
            end
            if ti==size(tvals_nsig,2) && tvals_nsig(ei,ti)==1 % end of array and cluster end
                clusterend = ti;
                elec_nclusters{ei}{ci} = clusterstart:clusterend;
            end
        end
    end
    clear clusterstart clusterend ci ei ti
    
    
    %% define (negative) clusters across electrodes
    tvals_nsig(tvals_nsig==1) = -1;
    % do not set clusternr = 0 because positive and negative clusters should have unique names
    for ei=1:numel(chanlocs) % loop over all electrodes
        for ci=1:numel(elec_nclusters{ei}) % loop over all elec_clusters of the current electrode
            
            % Is cluster already named?
            if tvals_nsig(ei,elec_nclusters{ei}{ci}(1)) == -1 % no, not yet named
                clusternr = clusternr +1; % increase by one for new cluster
                tvals_nsig(ei,elec_nclusters{ei}{ci}) = clusternr;
                %imagesc(tvals_nsig(elec_order,:)); colorbar; set(gca,'YTick',1:54,'YTickLabels',{EEG_chanlocs54(elec_order).labels}); drawnow; pause(1);
            end
            
            % Are there neighbours with overlapping clusters in time?
            for ni=1:numel(chanlocs) % loop over all other electrodes
                % find potential neighbours
                if ei==ni % do nothing, as this is the same electrode
                    % verify, whether electrode is direct neighbour
                elseif sqrt((chanlocs(ei).X - chanlocs(ni).X)^2 + ...
                        (chanlocs(ei).Y - chanlocs(ni).Y)^2 + ...
                        (chanlocs(ei).Z - chanlocs(ni).Z)^2) < max_elec_dist
                    
                    for nci=1:numel(elec_nclusters{ni}) % loop over all clusters of the neighbour
                        
                        % Does cluster overlap with the index cluster?
                        if intersect(elec_nclusters{ei}{ci},elec_nclusters{ni}{nci}) % yes, overlapping
                            % Is this neighbour-cluster already named?
                            if tvals_nsig(ni,elec_nclusters{ni}{nci}(1)) == -1 % no, not yet named
                                tvals_nsig(ni,elec_nclusters{ni}{nci}) = tvals_nsig(ei,elec_nclusters{ei}{ci}(1)); % name it according to index cluster
                                %imagesc(tvals_nsig(elec_order,:)); colorbar; set(gca,'YTick',1:54,'YTickLabels',{EEG_chanlocs54(elec_order).labels}); drawnow; pause(1);
                            else % yes, already named. Is this name different from the index cluster?
                                if tvals_nsig(ei,elec_nclusters{ei}{ci}(1)) ~= tvals_nsig(ni,elec_nclusters{ni}{nci}(1)) % if cluster number is different
                                    lower_cluster_nr  = min(tvals_nsig(ei,elec_nclusters{ei}{ci}(1)), tvals_nsig(ni,elec_nclusters{ni}{nci}(1)));
                                    higher_cluster_nr = max(tvals_nsig(ei,elec_nclusters{ei}{ci}(1)), tvals_nsig(ni,elec_nclusters{ni}{nci}(1)));
                                    tvals_nsig(tvals_nsig==higher_cluster_nr) = lower_cluster_nr;
                                    %imagesc(tvals_nsig(elec_order,:)); colorbar; set(gca,'YTick',1:54,'YTickLabels',{EEG_chanlocs54(elec_order).labels}); drawnow; pause(1);
                                end
                            end
                        else % no, clusters do not overlap, do nothing
                        end
                    end
                else % electrode is not a neighbour, do nothing
                end
            end
        end
    end
    clear ei ci ni nci elec_nclusters
    clear clusternr higher_cluster_nr lower_cluster_nr
    
        
    %% extract cluster statistic
    % merge positive and negative cluster
    tvals_sig = tvals_psig + tvals_nsig;
    clear tvals_psig tvals_nsig
    
    % change cluster numbers to 1-by-1 increasing cluster numbers
    tvals_sig = tvals_sig * (-1);
    clusternames = unique(tvals_sig);
    clusternames(clusternames==0) = [];
    for i=1:numel(clusternames)
        tvals_sig(tvals_sig==clusternames(i)) = i;
    end
    clear clusternames
    
    % define all clusters
    amount_cluster = unique(tvals_sig);
    amount_cluster(amount_cluster==0) = [];
    
    % define cluster sizes
    tmax   = nan(size(amount_cluster));
    clsize = nan(size(amount_cluster));
    for ci=1:numel(amount_cluster) % loop through all clusters
        tmax(ci)   = sum(tmp_tvals(tvals_sig==amount_cluster(ci)));
        clsize(ci) = numel(tmp_tvals(tvals_sig==amount_cluster(ci)));
    end
    
    
    %% save all relevant information
    cluster_result.amount_cluster = amount_cluster;
    cluster_result.tvals_sig      = tvals_sig;
    cluster_result.tmax           = tmax;
    cluster_result.clsize         = clsize;
    
    info.chanlocs = chanlocs;
    info.note     = note;
    info.note2    = note2;
    info.alpha    = alpha;
    info.tails    = tails;
    info.tcrit    = tcrit;
    info.volume   = vi;
    cluster_result.info = info;
    
    save([PATH_IN 'Exp09_CGR_ClusterResult_dist' num2str(max_elec_dist) '_vol' num2str(vi) '.mat'],'cluster_result')
    
end

