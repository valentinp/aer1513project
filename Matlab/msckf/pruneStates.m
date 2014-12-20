function [ prunedMsckfState, deletedCamStates ] = pruneStates( msckfState )
%PRUNESTATES Prunes any states that have no tracked features and updates
%covariances
    
    prunedMsckfState.imuState = msckfState.imuState;
    prunedMsckfState.imuCovar = msckfState.imuCovar;
    
    %Find all camStates with no tracked landmarks    
    deleteIdx = [];
    for c_i = 1:length(msckfState.camStates)
        if isempty(msckfState.camStates{c_i}.trackedFeatureIds)
            deleteIdx(end+1) = c_i;
        end
    end
    
    %Prune the damn states!
    
    deletedCamStates = msckfState.camStates(deleteIdx);
    prunedMsckfState.camStates = removeCells(msckfState.camStates, deleteIdx);
    
    keepStatesIdx = 1:size(msckfState.camCovar,1);
    keepStatesMask = true(1, numel(keepStatesIdx));
    for dIdx = deleteIdx
        keepStatesMask(6*dIdx - 5:6*dIdx) = false(6,1);
    end
    
    keepStatesIdx = keepStatesIdx(keepStatesMask);
    
    prunedMsckfState.camCovar = msckfState.camCovar(keepStatesIdx, keepStatesIdx);
    %Keep rows, prune columns of upper right covariance matrix
    prunedMsckfState.imuCamCovar = msckfState.imuCamCovar(:, keepStatesIdx);
     
end

