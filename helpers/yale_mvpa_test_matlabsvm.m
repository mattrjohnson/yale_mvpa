function [acts s]=test_matlabsvm_mrj( testpats, testtargs, svmstruct )

% testpats is nfeatures x ntrials (with multi-class, includes trials from all classes)
% testtargs is nclasses x ntrials (logical matrix)
% svmstruct is a struct of structs, with fields
%  svmstruct(i).struct -- what is returned by svmtrain for one pairwise SVM
%  svmstruct(i).conds  -- condition indices for the pairwise SVM, i.e. conditions 1 and 3 (work as indices into the rows of traintargs)

% acts will be nclasses x ntrials and represents "activations" for each class in each trial, for a winner-take-all type analysis
% s is included as an output value to match the other MVPA toolbox functions
%  - to match the toolbox, it should probably return a copy of svmstruct
%  - however, mrj_eeg_classify doesn't use it at all right now, so really we can stick whatever we want in it (that's what she said)
s=svmstruct; %for the moment; may want to return other stuff at some point

%not going to do a lot of error-checking to make sure # of classes match, but we'll do a tiny bit here
n_classes=size(testtargs,1);
if n_classes<2
    error('Need at least 2 classes!');
end
class_pairs=sortrows(combnk(1:n_classes,2));
if size(class_pairs,1) ~= length(svmstruct)
    error('Something is wrong -- not the same # of classes as training, maybe?');
end

acts=zeros(size(testtargs));
% acts_times_used=acts;
for i=1:length(svmstruct)
    these_conds=svmstruct(i).conds;
    
    trials_to_use = testtargs(these_conds(1),:) | testtargs(these_conds(2),:);
%     these_sample=(testpats(:,trials_to_use))';
    these_groups=(testtargs(these_conds(1),trials_to_use))'; %first class listed will be the 1's class, other will be the 0's class
    
%     [group_preds rawresults]=svmclassify(svmstruct(i).struct,these_sample);
    [group_preds rawresults]=svmclassify(svmstruct(i).struct,testpats');
    
%     acts(these_conds(1),trials_to_use)=acts(these_conds(1),trials_to_use)+rawresults';
%     acts(these_conds(2),trials_to_use)=acts(these_conds(2),trials_to_use)-rawresults';
    acts(these_conds(1),:)=acts(these_conds(1),:)-rawresults';
    acts(these_conds(2),:)=acts(these_conds(2),:)+rawresults';
    
%     acts_times_used(these_conds(1),trials_to_use)=acts_times_used(these_conds(1),trials_to_use)+1;
%     acts_times_used(these_conds(2),trials_to_use)=acts_times_used(these_conds(2),trials_to_use)+1;
    
%     s(i).accs=(these_groups==group_preds);
    s(i).accs=(these_groups==group_preds(trials_to_use));
end

% acts = acts ./ acts_times_used;