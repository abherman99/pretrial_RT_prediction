function [erdPredictModel_noConf, FStat_noConf, p_perm_noConf, erdPredictModel_confBin, FStat_confBin, p_perm_ConfBin, erdPredictModel_confType, FStat_ConfType, p_perm_confType] = MSIT_GLME_ERD(groupedData,groupedTable,start_times)
%written by Alexander B Herman 2022
%inputs:
% groupedData: tables of neural features for each trial (spike counts)
% aggregated across cells and participants, 
%groupedTable: table of covariates for all trials and all participants, one cell entry for each time-window to analyze
% with variables fr_precue_all_raw is the spike count in the window of interest, zLogRT is zscored reaction time,
% zLogLastRT is zscored last trial reaction time, prevConflictType is the
% conflict type on the last trial, subjVec are the participant numbers and
% cellNumVec are the cell numbers for random effects.
%participants, with one cell entry for each time-window to analyze


%define the models
glmeFormulaFR1 = 'fr_precue_all_raw ~ zLogRT + zLogLastRT + prevConflictType  + (1|subjVec) + (1|cellNumVec:subjVec)';  %No conflict term
glmeFormulaFR2 = 'fr_precue_all_raw ~ ConflictBin*zLogRT  + zLogLastRT + prevConflictBin + (1|subjVec) + (1|cellNumVec:subjVec)'; %binary term for any kind of conflict
glmeFormulaFR3 = 'fr_precue_all_raw ~ conflict*zLogRT + zLogLastRT + prevConflictType + (1|subjVec) + (1|cellNumVec:subjVec)'; %conflict specific term



for tw=1:length(start_times)

    [erdPredictModel_noConf{tw},FStat_noConf(:,tw),p_perm_noConf(:,tw)] = MSIT_GLME_permTest2(glmeFormulaFR1,'poisson','log',groupedData{tw},groupedTable{tw},'effects',excludeIncor,1000);
    [erdPredictModel_confBin{tw},FStat_confBin(:,tw),p_perm_confBin(:,tw)] = MSIT_GLME_permTest2(glmeFormulaFR2,'poisson','log',groupedData{tw},groupedTable{tw},'effects',excludeIncor,1000);
    [erdPredictModel_confType{tw},FStat_confType(:,tw),p_perm_confType(:,tw)] = MSIT_GLME_permTest2(glmeFormulaFR3,'poisson','log',groupedData{tw},groupedTable{tw},'effects',excludeIncor,1000);

end
end
