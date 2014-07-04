function results = yale_mvpa_display_process_results( yale_mvpa_config, subs, results ) %#ok<INUSL>

%results starts out as n_subs x 1 cell array
% each cell is a struct (of size 1 x nits) with fields:
%  - n_features: single double value
%  - acts: nconds x ntrials output array
%  - testtargs: nconds x ntrials output array (same as acts)

% n.b.: need to figure out what to do with by_timepoint (or by_channel) type analyses
% - presumably timepoints/channels will turn the struct from (1 x nits) to (channel x nits)


n_subs =                                length( results );
% wta_means =                             zeros( n_subs, n_timepoints );
% tiedrank_means =                        zeros( n_subs, n_timepoints );
% auc_means =                             zeros( n_subs, n_timepoints );
% svm_means =                             zeros( n_subs, n_timepoints );

n_timepoints_or_channels =              size( results{1}, 1 );
nits =                                  size( results{1}, 2 );
if n_timepoints_or_channels>1
    error('Not supported yet!');
end

wta_means =                             zeros(n_subs, 1);
tiedrank_means =                        zeros(n_subs, 1);
auc_means =                             zeros(n_subs, 1);
% svm_means =                             zeros(n_subs, 1);


% do some calculations on results
for i=1:n_subs
%     for j=1:n_timepoints
    for j = 1:nits
        results{i}(j).wta =                 yale_mvpa_results_wta( results{i}(j) );
        results{i}(j).tiedrank =            yale_mvpa_results_tiedrank( results{i}(j) );
        results{i}(j).auc =                 yale_mvpa_results_auc( results{i}(j) );

        wta_means(i,j) =                    mean(mean(squeeze(results{i}(j).wta)));
        tiedrank_means(i,j) =               mean(mean(squeeze(results{i}(j).tiedrank)));
        auc_means(i,j) =                    mean(mean(squeeze(results{i}(j).auc)));
        
        results{i}(j).wta_mean =            wta_means(i,j);
        results{i}(j).tiedrank_mean =       tiedrank_means(i,j);
        results{i}(j).auc_mean =            auc_means(i,j);
    end
%     end
end

wta_means = mean(wta_means,2);
tiedrank_means = mean(tiedrank_means,2);
auc_means = mean(auc_means,2);

if yale_mvpa_config.general.results_verbose
    disp(' '); disp(' '); disp(' ');
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    disp('% RESULTS');
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    disp(' ');
    disp('WTA:');
    disp(wta_means);
    disp(' ');
    disp('Tiedrank:');
    disp(tiedrank_means);
    disp(' ');
    disp('AUC:');
    disp(auc_means);
    disp(' ');
end


% %-----------------------------------------------------------------
% %-----------------------------------------------------------------
% function svm_out = mrj_eeg_classify_process_results_svm( results )
% % svm_out should be (# of 2-condition combinations) x (# of trials per condition in test set)
% %  values: 1 or 0, indicating whether SVM was right for each binary classifier
% %  (chance is 50%, if you average all the 1's and 0's together)
% %  unlike wta and tiedrank, we'll go ahead and average over iterations here to keep things to a 2-D matrix
% 
% nits =                                  length(results.acts);
% svm_out =                               zeros(size( results.svm_orig_accs{1} ));
% 
% for i=1:nits
%     svm_out =                           svm_out + results.svm_orig_accs{i};
% end
% 
% svm_out =                               svm_out / nits;
