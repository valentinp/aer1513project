function [ prunedMsckfState ] = pruneStates( msckfState )
%PRUNESTATES Prunes any states that have no tracked features and updates
%covariances
    
    prunedMsckfState.imuState = msckfState.imuState;
    prunedMsckfState.imuCovar = msckfState.imuCovar;
    
    %Find all camStates with no tracked landmarks    
    pruneStateBinaryFlags = zeros(1, length(msckfState.camStates));
    
    for c_i = 1:length(msckfState.camStates)
        if isempty(msckfState.camStates{c_i}.trackedFeatureIds)
            pruneStateBinaryFlags(c_i) = 1;
        end
    end
    
    %Prune the damn states!
    keepStatesIdx = (pruneStateBinaryFlags == 0);
    
    prunedMsckfState.camStates = msckfState.camStates(keepStatesIdx);
    prunedMsckfState.camCovar = msckfState.camCovar(keepStatesIdx, keepStatesIdx);
    %Keep rows, prune columns of upper right covariance matrix
    prunedMsckfState.imuCamCovar = msckfState.imuCamCovar(:, keepStatesIdx);
     
end

