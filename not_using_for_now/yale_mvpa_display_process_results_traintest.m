function results = mrj_eeg_classify_display_process_results_traintest( subs, results )

%OK, actually just modifying this whole thing for Img_Perc fMRI
%copying most of it over from main display_process_results file, and using some cheap hacks

do_svm =                                isfield( results{1}, 'svm_orig_accs' );
n_subs =                                length( results );
n_timepoints =                          length( results{1} );
wta_means =                             zeros( n_subs, n_timepoints );
tiedrank_means =                        zeros( n_subs, n_timepoints );
auc_means =                             zeros( n_subs, n_timepoints );
svm_means =                             zeros( n_subs, n_timepoints );

% do some calculations on results
for i=1:n_subs
    results{i}.testtargs{1} = kron(eye(4),ones(1,size(subs(i).trainset.classifier_data{1},2))); %main cheap hack for Img_Perc fMRI
    for j=1:n_timepoints
        results{i}(j).wta =                 mrj_eeg_classify_process_results_wta( results{i}(j) ); % indexed timepoint, iteration, trial
        results{i}(j).tiedrank =            mrj_eeg_classify_process_results_tiedrank( results{i}(j) ); % ditto
        results{i}(j).auc =                 mrj_eeg_classify_process_results_auc( results{i}(j) ); % ditto

        wta_means(i,j) =                    mean(mean(squeeze(results{i}(j).wta)));
        tiedrank_means(i,j) =               mean(mean(squeeze(results{i}(j).tiedrank)));
        auc_means(i,j) =                    mean(mean(squeeze(results{i}(j).auc)));
        
%         if do_svm
%             results{i}(j).svm_accs =        mrj_eeg_classify_process_results_svm( results{i}(j) ); % indexed timepoint, combination, trial
%             svm_means(i,j) =                mean(mean(squeeze(results{i}(j).svm_accs)));
%         end
    end
end

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
if do_svm
    disp('SVM:');
    disp(svm_means);
    disp(' ');
end
disp('AUC:');
disp(auc_means);
disp(' ');

% keyboard;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KEY SUB-FUNCTION(S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-----------------------------------------------------------------
%-----------------------------------------------------------------
function wta_out = mrj_eeg_classify_process_results_wta( results )
% wta_out will be nits x ntrials matrix of 1's (winner == correct) and 0's
%  (chance is 25%, if you average all the 1's and 0's together)

nits =                                  length(results.acts);
wta_out =                               zeros( nits, size(results.acts{1},2) );

for i=1:nits
    acts =                              results.acts{i};
    testtargs =                         results.testtargs{i};
    acts_maxes =                        max(acts);
    acts_tmp =                          zeros(size(acts));
    for j=1:length(acts_maxes)
        acts_tmp(:,j) =                 (acts(:,j)==acts_maxes(j));
    end
    acts_tmp2 =                         (acts_tmp & testtargs);
    acts_accs_wta =                     sum(acts_tmp2);
    acts_accs_wta(acts_accs_wta>0) =    1./sum(acts_tmp(:,acts_accs_wta>0)); %adjusts for ties -- trials with ties get fractional accuracies
    wta_out(i,:) =                      acts_accs_wta;
end


%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
function tiedrank_out = mrj_eeg_classify_process_results_tiedrank( results )
% tiedrank_out should be nits x ntrials, containing rank score for the correct category on each trial
%  (scale is 1-4, with 4 being best; chance is 2.5)

nits =                                  length(results.acts);
tiedrank_out =                          zeros( nits, size(results.acts{1},2) );

for i=1:nits
    acts =                              results.acts{i};
    testtargs =                         logical(results.testtargs{i});
    tiedrank_tmp =                      tiedrank(acts);
    tiedrank_out(i,:) =                 tiedrank_tmp(testtargs);
end


%-----------------------------------------------------------------
%-----------------------------------------------------------------
function svm_out = mrj_eeg_classify_process_results_svm( results )
% svm_out should be (# of 2-condition combinations) x (# of trials per condition in test set)
%  values: 1 or 0, indicating whether SVM was right for each binary classifier
%  (chance is 50%, if you average all the 1's and 0's together)
%  unlike wta and tiedrank, we'll go ahead and average over iterations here to keep things to a 2-D matrix

nits =                                  length(results.acts);
svm_out =                               zeros(size( results.svm_orig_accs{1} ));

for i=1:nits
    svm_out =                           svm_out + results.svm_orig_accs{i};
end

svm_out =                               svm_out / nits;


%-----------------------------------------------------------------
%-----------------------------------------------------------------
function auc_out = mrj_eeg_classify_process_results_auc( results )
% auc_out will be n_conds x nits matrix of AUC values (should be > .5 if better than chance)

nits =                                  length( results.acts );
n_conds =                               size( results.acts{1},1 );
auc_out =                               zeros( n_conds, nits );

for i=1:nits
    acts =                              results.acts{i};
    testtargs =                         results.testtargs{i};
    
    for j=1:n_conds
        [tp fp] =                       roc( testtargs(j,:)', acts(j,:)' );
        auc_out(j,i) =                  auroc( tp, fp );
    end
end


% for i=1:length(subs) %n_subs
%     %changed for Img_Perc fMRI
% %     sub_ntpc(i,:)=subs(i).testset.n_trials_per_cond; %#ok<AGROW>
%     sub_ntpc(i,:) = ones(subs(i).testset(1).n_conds,1)  * size(subs(i).testset(1).classifier_data{1},2); %#ok<AGROW>
% end
% 
% for i=1:length(subs)
%     end_inds=cumsum(sub_ntpc(i,:));
%     start_inds=[1,(end_inds(1:end-1)+1)];
%     
%     for j=1:length(subs(i).trainset) %n train sets
%         these_acts=[];
%         for k=1:5 %classifier nits
%             these_acts(k,:,:)=results{i}(j).acts{k}; %#ok<AGROW>
%         end
%         these_acts=squeeze(mean(these_acts,1)); %now just n_train_conds (for this particular training set) x n_trials
% 
%         for k=1:length(start_inds) %n test conditions
%             this_testcond_acts=these_acts(:,start_inds(k):end_inds(k));
%             results_tt(i,j,k,:)=mean(this_testcond_acts,2); %#ok<AGROW,NASGU>
%         end
%     end
% 
% end
% 
% conds=subs(1).testset.cond_names'; %#ok<NASGU>
% 
% keyboard;
