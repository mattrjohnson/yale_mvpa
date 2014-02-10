function svmstruct = yale_mvpa_train_matlabsvm( trainpats, traintargs, args )

% trainpats is nfeatures x ntrials (with multi-class, includes trials from all classes)
% traintargs is nclasses x ntrials (logical matrix)

% args: may not be needed yet?

% svmstruct will be a struct of structs, with fields
%  svmstruct(i).struct -- what is returned by svmtrain for one pairwise SVM
%  svmstruct(i).conds  -- condition indices for the pairwise SVM, i.e. conditions 1 and 3 (work as indices into the rows of traintargs)



%first: break into pairs of classes/conditions
n_classes=size(traintargs,1);
if n_classes<2
    error('Need at least 2 classes!');
end
class_pairs=sortrows(combnk(1:n_classes,2));

for i=1:size(class_pairs,1)
    these_conds=class_pairs(i,:);
    svmstruct(i).conds=these_conds; %#ok<AGROW>
    
    trials_to_use = traintargs(these_conds(1),:) | traintargs(these_conds(2),:);
    these_training = (trainpats(:,trials_to_use))';
    these_groups = (traintargs(these_conds(1),trials_to_use))'; %first class listed will be the 1's class, other will be the 0's class
    
    svmstruct(i).struct = svmtrain( these_training, these_groups ); %#ok<AGROW>
end


% function [svm_struct, svIndex] = svmtrain(training, groupnames, varargin)
%training is trials x features
%groupnames is trials x 1, can be anything but for simplicity we'll use 1's and 0's