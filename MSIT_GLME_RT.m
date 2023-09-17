function rtPredictModel = MSIT_GLME_RT(groupedData,groupedTable,start_times)
%written by Alexander B Herman 2022
%inputs:
% groupedData: tables of neural features for each trial (spike counts)
% aggregated across cells and participants
% 
% %groupedTable: table of covariates for all trials and all participants, one cell entry for each time-window to analyzewith variables RTall are the raw reaction times, fr_precue_all is the z-scored spike count in the window of interest, zLogRT is zscored reaction time,
% zLogLastRT is zscored last trial reaction time, prevConflictType is the
%fr_precue_all_meanVec is the mean firing rate across trials for each cell
% conflict type on the last trial, subjVec are the participant numbers and
% cellCatVec is a categorical vector of cell numbers, eriksenTrials is a
% dummy variable = 1 if the trial had eriksen conflict, simonTrials = 1 for
% the presence of simon conflict
% with one cell entry for each time-window to analyze




linkFunc='identity'; %link function
dummyvarcoding='effects'; %coding for categorical variables
excludeIncor=0; % set to 1 if want only correct trials


glmeFormula_RT_1='RTall ~  eriksenTrials*fr_precue_all*cellCatVec + simonTrials*fr_precue_all*cellCatVec + fr_precue_all_meanVec + prevConflictType + zLogLastRT - eriksenTrials:cellCatVec - simonTrials:cellCatVec - cellCatVec + (1|subjVec)'; %model to test for single cell effects in predicting RT


glmeFormula_RT_2='RTall ~  eriksenTrials*fr_precue_all + simonTrials*fr_precue_all + fr_precue_all_meanVec + prevConflictType + zLogLastRT + trialStartTime + (1|subjVec)'; %this equation tests for effects of time into experiment on a region level.



for tw=1:length(start_times)
    rtPredictModel{tw}= MSIT_GLME(glmeFormula_RT_1,'normal','log',groupedData{tw},groupedTable{tw},'effects',excludeIncor);

end
end
