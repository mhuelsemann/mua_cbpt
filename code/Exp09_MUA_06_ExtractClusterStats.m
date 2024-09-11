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


%% Extraction of Embodiment Cluster Statistics for Publication

% load original ANOVA and cluster stats
PATH_IN = '';
trialamount = 'same';
effectnames = {'CGR-N400','CGRxVOL-N400', ...
    'CGR-N180','CGRxVOL-N180', 'CGR-P280','CGRxVOL-P280', ...
    'VOL-N180','VOL-P280','VOL-N400'};
nperm = 10000;
alpha = .05;
timews = {400:100:1700;400:100:1700; 130; 130; 230; 230; 130; 230; 400:100:1700};

for ef=1:numel(effectnames) % loop over all effects
    timew = timews{ef,:};
    
    % load ANOVA and cluster statistics
    load([PATH_IN trialamount '\Exp09_ANOVA_OrigCluster_' trialamount 'AmountTrials_' effectnames{ef} '.mat'])
    load([PATH_IN trialamount '\Exp09_ANOVA_PermCluster_' trialamount 'AmountTrials_' effectnames{ef} '.mat'])
    
    % calculate p-val
    pvals = nan(cluster.amount_cluster,2);
    for i=1:cluster.amount_cluster
        pvals(i,1) = sum(cluster_permuted.Fmax   > abs(cluster.Fmax(i))  )/nperm;
        pvals(i,2) = sum(cluster_permuted.clsize > abs(cluster.clsize(i)))/nperm;
    end
    cluster.pvals_Fmax   = pvals(:,1);
    cluster.pvals_clsize = pvals(:,2);
    save([PATH_IN 'Exp09_ANOVA_OrigCluster_' trialamount 'AmountTrials_' effectnames{ef} '_v2.mat'],'cluster')
    
    % save figure of permuted cluster distribution
    figure('Units','Centimeters','Position',[5 5 20 10]);
    subplot(131); histogram(cluster_permuted.amount_cluster); xlabel('Amount Cluster'); ylabel('Count'); title(effectnames{ef}); ylim([0 5000]);
    subplot(132); histogram(cluster_permuted.Fmax          ); xlabel('F-Mass');         title(effectnames{ef}); ylim([0 5000]);
    subplot(133); histogram(cluster_permuted.clsize        ); xlabel('Cluster Size');   title(effectnames{ef}); ylim([0 5000]);
    saveas(gcf, [PATH_IN 'Exp09_ANOVA_PermCluster_' trialamount 'AmountTrials_' effectnames{ef} '.png'])
    close(gcf)
    
    % print out stats in command window for each cluster
    % open text file to save results for publication
    fileID = fopen([PATH_IN 'Exp09_ANOVA_OrigCluster_' trialamount 'AmountTrials_' effectnames{ef} '.txt'],'w');
    for i=1:cluster.amount_cluster
        
        % test whether current cluster is (marginally) significant
        if cluster.pvals_Fmax(i)>=2*alpha && cluster.pvals_clsize(i)>=2*alpha
            fprintf(fileID, 'Cluster %d is not significant.\r\n\r\n',i);
        else
            fprintf(fileID, 'Cluster %d is (marginally) significant.\r\n',i);
            fprintf(fileID, 'F-mass: %0.2f, p=%0.3f;\r\ncluster size: %d, p=%0.3f.\r\n',cluster.Fmax(i),cluster.pvals_Fmax(i),cluster.clsize(i),cluster.pvals_clsize(i));
            fprintf(fileID, 'The following electrodes contribute to the cluster:\r\n');
            
            [contr_elecs(:,1), contr_elecs(:,2)] = (find(cluster.Fvals_sig(:,:) == i)); 
            Fvals = nan(size(contr_elecs,1), 1); wvals = nan(size(contr_elecs,1), 1); 
            for j=1:size(contr_elecs,1)
                Fvals(j,1) = cluster.stats(contr_elecs(j,1), 1, contr_elecs(j,2));
                wvals(j,1) = cluster.stats(contr_elecs(j,1), 6, contr_elecs(j,2));
                fprintf(fileID, '%s, %s-%sms: F(%d,%d) = %0.2f, p = %0.3f, omega = %0.2f (e=%0.2f, eta=%0.2f)\r\n', ...
                    cluster.info.chanlocs(contr_elecs(j,1)).labels, ... % elec
                    num2str(timew(contr_elecs(j,2))), ... % timewindow start
                    num2str(timew(contr_elecs(j,2))+100), ... % timewindow end
                    cluster.stats(contr_elecs(j,1), 2, contr_elecs(j,2)), ... % dfZ
                    cluster.stats(contr_elecs(j,1), 3, contr_elecs(j,2)), ... % dfN
                    cluster.stats(contr_elecs(j,1), 1, contr_elecs(j,2)), ... % F
                    cluster.stats(contr_elecs(j,1), 5, contr_elecs(j,2)), ... % pval
                    cluster.stats(contr_elecs(j,1), 6, contr_elecs(j,2)), ... % wval
                    cluster.stats(contr_elecs(j,1), 4, contr_elecs(j,2)), ... % epsilon
                    cluster.stats(contr_elecs(j,1), 7, contr_elecs(j,2)));    % eta
            end            
            
            max_elec = find(Fvals==max(Fvals));
            fprintf(fileID, 'The electrode with the largest effect is %s at %s-%sms.\r\n', ...
                cluster.info.chanlocs(contr_elecs(max_elec,1)).labels, ...
                num2str(timew(contr_elecs(max_elec,2))), ... % timewindow start
                num2str(timew(contr_elecs(max_elec,2))+100)); % timewindow end
                
            fprintf(fileID, 'Effect sizes for the contributing electrodes ranged between %0.2f < omega < %0.2f, with an average (plusminus SD) of omega = %0.2f plusminus %0.2f.\r\n\r\n', ...
                min(wvals), max(wvals), mean(wvals), std(wvals));
            
            clear contr_elecs Fvals wvals max_elec
        end
    end
    fclose(fileID);
    
    clear cluster* pvals i fileID contr_elecs ti stats j max_elec Fvals wvals
end
