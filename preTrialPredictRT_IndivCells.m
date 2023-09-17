function preTrialPredictRT_IndivCells(datafile,datavarname,window_type,window_size,epoch_starts,epoch_end,epoch_postcue_starts,epoch_postcue_end,dist,linkfunc,loglastRT,zscore_lastRT,cor_uncor_flag,stim_type,exclude_prob)
zscore_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));

%%written by Alexander B Herman 2022
%%This function takes as inputs a fieltrip-style cell array of individual
%%neurons with spike counts and trial by trial behavioral variables. For
%%each cell it computes a model that predicts reaction times from firing
%%rates, and then it computes the angle and correlations between the
%%coefficient vectors. 
%%example" preTrialPredictRT_makeData('cellData.mat','growing_pre',50,[800],1000,1000,1300,'normal','log',0,1,1,'both',0.005)

correct_flag=0;
Zscore_flag=0;
logRTFlag=0;
trialThresh=5;

dbstop if error

load(datafile);


load(datafile);
switch datavarname
    case 'neuralData'
        data=neuralData;
        clear neuralData
    case 'data_cing_all'
        data=data_cing_all;
        clear data_cing_all
    case 'alldata_cing'
        % data=data_cing_all;
        data=alldata_cing;
        clear alldata_cing
end



%data=datavarname;
numPerms=5000;
%epoch_end=1000;
init_end=1300;
init_start = init_end-window_size;
epoch_length=50;
step=2;
stepsize=5;
% if window_size==50
%     max_steps = 250;
% end
MAD_Factor = 3;
max_steps=100;
excludeFactor=100;
Beta_thresh=2;
Bonf_p_value=0.025;
switch window_type
    %If the goal is to examine a "growing window"
    case 'growing_pre'
        % max_steps = 2;
        %for jj=1:max_steps
        for jj=1:length(epoch_starts)
            
            % epoch_start=init_start-(jj-1)*stepsize;
            %epoch_end=init_end-(jj-1)*stepsize;
            epoch_analyze{jj}=epoch_starts(jj):epoch_end;
            %epoch_analyze_postcue=epoch_postcue_starts:epoch_postcue_end;
            epoch_analyze_postcue{jj}=epoch_postcue_starts:epoch_postcue_end;
            %  epoch_analyze_all{c,i}=epoch_analyze;
        end
        num_epochs=length(epoch_starts);
    
            case 'growing_post'
        % max_steps = 2;
        %for jj=1:max_steps
        for jj=1:length(epoch_postcue_end)
            
            % epoch_start=init_start-(jj-1)*stepsize;
            %epoch_end=init_end-(jj-1)*stepsize;
            epoch_analyze{jj}=epoch_starts:epoch_end;
            %epoch_analyze_postcue=epoch_postcue_starts:epoch_postcue_end;
            epoch_analyze_postcue{jj}=epoch_postcue_starts(jj):epoch_postcue_end(jj);
            %  epoch_analyze_all{c,i}=epoch_analyze;
        end
        num_epochs=length(epoch_postcue_starts);
        
        %if the goal is to examine a sliding window of fixed size
    case 'sliding'
        % max_steps = 2;
        for jj=1:max_steps
            epoch_start=init_start-(jj-1)*stepsize;
            epoch_end=init_end-(jj-1)*stepsize;
            epoch_analyze{jj}=epoch_start:epoch_end;
            epoch_postcue_start=init_end + (jj-1)*stepsize;
            epoch_postcue_end=init_end + (jj-1)*stepsize + window_size;
            epoch_analyze_postcue{jj}=epoch_postcue_start:epoch_postcue_end;
            %  epoch_analyze_all{init_end + (jj-1)*stepsize;c,i}=epoch_analyze;
        end
        num_epochs=max_steps;
end
num_cells=length(data);

for c = 1:num_cells
  % for  c = 1:4
    %disp(c)
    warning('off','all')
    subj{c}=data{c}.filename(1:7);
    conflict = data{c}.vars(:,9);
    RT = data{c}.vars(:,15);
    RTorig=RT;
    lastRTorig=[NaN; RTorig(1:end-1)];
    zLogRT=zscore(log(RT));
    zLogLastRT=[NaN; zLogRT(1:end-1)];
    switch logRTFlag
        case 1
            RT=log(RT);
            switch Zscore_flag
                case 1
                    RT=zscore_xnan(RT);
                case 0
            end
            lastRT=[NaN; RT(1:end-1)];
            nextRT=[RT(2:end); NaN];
            lastRTorig=[NaN; RTorig(1:end-1)];
        case 0
            if loglastRT==1
                temp_RT=log(RT);
                
            else
                temp_RT=RT;
            end
            if zscore_lastRT==1
                temp_RT=zscore_xnan(RT);
                lastRT=[NaN; temp_RT(1:end-1)];
                lastRTorig=[NaN; RTorig(1:end-1)];
                % nextRT=[temp_RT(2:end); NaN];
            else
                lastRT=[NaN; temp_RT(1:end-1)];
                %nextRT=[temp_RT(2:end); NaN];
                
            end
            
            %
            %             switch Zscore_flag
            %                 case 1
            %                     RT=zscore(RT);
            %                 case 0
            %             end
    end
    
        temp_RT=zscore_xnan(log(RT));
        logRT=log(RTorig);
        logLastRT=log(lastRTorig);
    zLogRT=temp_RT;
    zLogLastRT=[NaN; temp_RT(1:end-1)];
    RTdiff=zscore_xnan(log(RT)-log(lastRTorig));
    CorrectTrials = data{c}.vars(:,16);
    IncorrectTrials = ~CorrectTrials;
    numCorrectTrials=sum(CorrectTrials);
    numIncorrecTrials=sum(IncorrectTrials);
    simonResponse = data{c}.vars(:,12);
    eriksenResponse = data{c}.vars(:,13);
    simonTrials = ~isnan(simonResponse);
    eriksenTrials = ~isnan(eriksenResponse);
    simonTrialsExc = ~isnan(simonResponse) & isnan(eriksenResponse);
    eriksenTrialsExc = ~isnan(eriksenResponse) & isnan(simonResponse);
    bothTrials = ~isnan(eriksenResponse) & ~isnan(simonResponse);
    noneTrials = isnan(eriksenResponse) & isnan(simonResponse);
    correctSimonRate(c) = sum(CorrectTrials & simonTrials)/sum(simonTrials);
    correctEriksenRate(c) = sum(CorrectTrials & eriksenTrials)/sum(eriksenTrials);
    %    noneTrials = isnan(simonResponse) & (eriksenResponse);
    prevEriksenTrials=[0; eriksenTrials(1:end-1)];
    prevSimonTrials=[0; simonTrials(1:end-1)];
    prevCorrect=[0; CorrectTrials(1:end-1)];
    psth_cue = data{c}.psth_cue;
    psth_response = data{c}.psth_response;
    psth_fixation = data{c}.psth_fixation;
    psth_cue_all{c}=psth_cue;
    psth_fixation_all{c}=psth_fixation;
    % normalize
    % npsth_fix = (psth_fix - mean(p:))) ./ std(psth_fix(:));
    %   npsth_cue = (psth_cue - mean(psth_cue(:))) sth_fix(./ std(psth_cue(:));
    %   npsth_cue_allCells{c}=npsth_cue;
    
    
    
    lastConflict = [NaN; conflict(1:end-1)];
    prevCorrect=[NaN; CorrectTrials(1:end-1)];
    prevIncorrect = [NaN; ~CorrectTrials(1:end-1)];
    numPrevIncor(c)=nansum(prevIncorrect);
    % ExcludeCellsRT=~(RT<mean(RT)+excludeFactor*std(RT)) & (RT>(mean(RT)-excludeFactor*std(RT)));
    ExcludeCellsRT=[];
    
    ExcludeCellsAll{c}=~(RT<mean(RT)+excludeFactor*std(RT)) & (RT>(mean(RT)-excludeFactor*std(RT)));
    numExcludeCells{c}=sum(~(RT<mean(RT)+excludeFactor*std(RT)) & (RT>(mean(RT)-excludeFactor*std(RT))));
    AllTrialsFinal= ~isnan(lastConflict) & ~isnan(lastRT) & (RT<(mean(RT)+excludeFactor*std(RT)) & (RT>mean(RT)-excludeFactor*std(RT)));
    %SimonTrialsFinalAll{c}= ~isnan(lastConflict) & ~isnan(lastRT) & simonTrials & (RT<(mean(RT)+excludeFactor*std(RT)) & (RT>mean(RT)-excludeFactor*std(RT)));
    % EriksenTrialsFinalAll{c}= ~isnan(lastConflict) & ~isnan(lastRT) & eriksenTrials & (RT<(mean(RT)+excludeFactor*std(RT))) & (RT>mean(RT)-excludeFactor*std(RT));
    % AllTrialsFinalAll{c}= ~isnan(lastConflict) & ~isnan(lastRT) & (RT<(mean(RT)+excludeFactor*std(RT))) & (RT>mean(RT)-excludeFactor*std(RT));
    SimonTrialsFinal= ~isnan(lastConflict) & ~isnan(lastRT) & simonTrials & (RT<(mean(RT)+excludeFactor*std(RT))) & (RT>mean(RT)-excludeFactor*std(RT));
    SimonTrialsFinalExc= ~isnan(lastConflict) & ~isnan(lastRT) & simonTrials & ~eriksenTrials & (RT<(mean(RT)+excludeFactor*std(RT))) & (RT>mean(RT)-excludeFactor*std(RT));
    
    % SimonTrialsFinalNEW= ~isnan(lastConflict) & ~isnan(lastRT) & simonTrialsNEW & (RT<(mean(RT)+excludeFactor*std(RT))) & (RT>mean(RT)-excludeFactor*std(RT));
    % SimonTrialsFinalNEWAll{c}= ~isnan(lastConflict) & ~isnan(lastRT) & simonTrialsNEW & (RT<(mean(RT)+excludeFactor*std(RT))) & (RT>mean(RT)-excludeFactor*std(RT));
    
    EriksenTrialsFinal= ~isnan(lastConflict) & ~isnan(lastRT) & eriksenTrials & (RT<(mean(RT)+excludeFactor*std(RT))) & (RT>mean(RT)-excludeFactor*std(RT));
    EriksenTrialsFinalExc= ~isnan(lastConflict) & ~isnan(lastRT) & eriksenTrials &~simonTrials & (RT<(mean(RT)+excludeFactor*std(RT))) & (RT>mean(RT)-excludeFactor*std(RT));
    %EriksenTrialsFinalNEW= ~isnan(lastConflict) & ~isnan(lastRT) & eriksenTrialsNEW & (RT<(mean(RT)+excludeFactor*std(RT))) & (RT>mean(RT)-excludeFactor*std(RT));
    % EriksenTrialsFinalNEWAll{c}= ~isnan(lastConflict) & ~isnan(lastRT) & eriksenTrialsNEW & (RT<(mean(RT)+excludeFactor*std(RT))) & (RT>mean(RT)-excludeFactor*std(RT));
    
    NoneTrialsFinal=isnan(eriksenResponse) & isnan(simonResponse) & ~isnan(lastConflict) & ~isnan(lastRT) & (RT<(mean(RT)+excludeFactor*std(RT))) & (RT>mean(RT)-excludeFactor*std(RT));
    % NoneTrialsFinalNEW= ~isnan(lastConflict) & ~isnan(lastRT) & ~simonTrialsNEW & ~eriksenTrialsNEW & (RT<(mean(RT)+excludeFactor*std(RT))) & (RT>mean(RT)-excludeFactor*std(RT));
    
    BothTrialsFinal=~isnan(eriksenResponse) & ~isnan(simonResponse) & ~isnan(lastConflict) & ~isnan(lastRT) & (RT<(mean(RT)+excludeFactor*std(RT))) & (RT>mean(RT)-excludeFactor*std(RT));
    EitherTrials=~isnan(eriksenResponse) | ~isnan(simonResponse);
    EitherTrialsFinal=~isnan(eriksenResponse) | ~isnan(simonResponse) & ~isnan(lastConflict) & ~isnan(lastRT) & (RT<(mean(RT)+excludeFactor*std(RT))) & (RT>mean(RT)-excludeFactor*std(RT));
    EitherTrialsFinallAll{c}=EitherTrialsFinal;
    AllTrialsFinalAll{c}=AllTrialsFinal;
    SimonTrialsFinalAll{c}=SimonTrialsFinal;
    SimonTrialsFinaExclAll{c}=SimonTrialsFinalExc;
    EriksenTrialsFinalAll{c}=EriksenTrialsFinal;
    EriksenTrialsFinalExcAll{c}=EriksenTrialsFinalExc;
    NoneTrialsFinalAll{c} = NoneTrialsFinal;
    BothTrialsFinalAll{c} = BothTrialsFinal;
    RTnone=RT(NoneTrialsFinal);
    
    CorrectAll=CorrectTrials(AllTrialsFinal);
    CorrectSimon=CorrectTrials(SimonTrialsFinal);
    CorrectEriksen=CorrectTrials(EriksenTrialsFinal);
    CorrectBoth=CorrectTrials(BothTrialsFinal);
    CorrectNone=CorrectTrials(NoneTrialsFinal);
    
    prevEriksenTrialsAll=prevEriksenTrials(AllTrialsFinal);
    prevSimonTrialsAll=prevSimonTrials(AllTrialsFinal);
    prevConflictBinary = prevEriksenTrialsAll | prevSimonTrialsAll;
    %prevEriksenTrials(1)=0;
    %prevSimonTrials(1)=0;
    prevConflictBinary = prevEriksenTrials | prevSimonTrials;
    prevConflict= prevEriksenTrials | prevSimonTrials;
    noPrevConflict = ~prevEriksenTrials & ~ prevSimonTrials;
    %prevConflictCatFinal = lastConflict(AllTrialsFinal);
    prevConflictCatFinal = (prevEriksenTrials*1 + prevSimonTrials*2 +1);
    prevConflictCat= (prevEriksenTrials*1 + prevSimonTrials*2 +1);
    ConflictCatFinal = (EriksenTrialsFinal*1 + SimonTrialsFinal*2 + 1);
    ConflictCat = eriksenTrials*1 + simonTrials*2 +1;
    ConflictBin = eriksenTrials | simonTrials;
    NoConflictBin = ~eriksenTrials & ~simonTrials;
    ConflictBinFinal = EriksenTrialsFinal | SimonTrialsFinal;
    %ConflictCatFinal = conflict(AllTrialsFinal);
    
    CongTrialsBin=((prevSimonTrials | prevEriksenTrials) & (simonTrials | eriksenTrials)) | ((~prevSimonTrials & ~prevEriksenTrials) & (~simonTrials & ~eriksenTrials));
    InCongTrialsBin=((prevSimonTrials | prevEriksenTrials) & (~simonTrials & ~eriksenTrials)) | ((~prevSimonTrials & ~prevEriksenTrials) & (simonTrials | eriksenTrials));
    

    RTall=RT;
    lastRTall=lastRT;
    nextRTall=[RTall(2:end); NaN];
    NextEriksenTrials=[eriksenTrials(2:end); NaN];
    NextSimonTrials=[simonTrials(2:end); NaN];
    % end
    meanRTEriksen=mean(RTall(EriksenTrialsFinalExc));
    meanRTSimon=mean(RTall(SimonTrialsFinalExc));
    meanRTBoth=mean(RTall(BothTrialsFinal));
    meanRTNone=mean(RTall(NoneTrialsFinal));
    
    
    RTAllSimon=RTall.*simonTrials;
    RTAllEriksen=RTall.*eriksenTrials;
    RTAllNoConf=RTall.*noneTrials;
    RTAllSimonExc=RTall.*(simonTrials & ~eriksenTrials);
    RTAllEriksenExc=RTall.*(eriksenTrials & ~ simonTrials);
    RTAllEither=RTall.*(simonTrials | eriksenTrials);
    RTAllNone=RTall.*(~simonTrials & ~eriksenTrials);
    RTAllCong=RTall.*CongTrialsBin;
    RTAllInCong=RTall.*InCongTrialsBin;

    RTAllBoth=RTall.*(simonTrials & eriksenTrials);
    
    
    RTAllall{c}=RTall;
    RTAllSimonAll{c}=RTAllSimon;
    RTAllEriksenAll{c}=RTAllEriksen;
    RTAllEitherAll{c}=RTAllEither;
    RTAllNoneAll{c}=RTAllNone;
    RTAllSimonExcAll{c}=RTAllSimonExc;
    RTAllEriksenExcAll{c}=RTAllEriksenExc;
    RTAllBothAll{c}=RTAllBoth;
    
    RTRawall{c}=RTorig;
    RTRawSimon=RTorig.*SimonTrialsFinal;
    RTRawEriksen=RTorig.*EriksenTrialsFinal;
    RTRawEither=RTorig.*(SimonTrialsFinal | EriksenTrialsFinal);
    RTRawNone=RTorig.*(~SimonTrialsFinal & ~EriksenTrialsFinal);
    RTRawSimonExc=RTorig.*SimonTrialsFinalExc;
    RTRawEriksenExc=RTorig.*EriksenTrialsFinalExc;
    RTRawBoth=RTorig.*BothTrialsFinal;
    
    RTRawSimonAll{c}=RTRawSimon;
    RTRawEriksenAll{c}=RTRawEriksen;
    RTRawEitherAll{c}=RTRawEither;
    RTRawNoneAll{c}=RTRawNone;
    RTRawSimonExcAll{c}=RTRawSimon;
    RTRawEriksenExcAll{c}=RTRawEriksen;
    RTRawBothAll{c}=RTRawBoth;
    
    
    
    npsth_cue = (psth_cue - mean(psth_cue(:))) ./ std(psth_cue(:));
    
    psth_cue_Simon{c} = psth_cue.*SimonTrialsFinal;
    psth_cue_Eriksen{c} = psth_cue.*EriksenTrialsFinal;

    psth_cue_SimonExc{c} = psth_cue.*SimonTrialsFinalExc;
    psth_cue_EriksenExc{c} = psth_cue.*EriksenTrialsFinalExc;
    
    for ii=1:num_epochs
        switch window_type
            case 'growing_pre'
        epoch_analyze_loc=cell2mat(epoch_analyze(ii));
        %        epoch_analyze_postcue_loc=cell2mat(epoch_analyze_postcue);
         epoch_analyze_postcue_loc= cell2mat(epoch_analyze_postcue(1));
            case 'growing_post'
        epoch_analyze_postcue_loc=cell2mat(epoch_analyze_postcue(ii));
        epoch_analyze_loc=cell2mat(epoch_analyze(1));
        end
        
        epoch_precue=epoch_analyze_loc;
        %SpikeCount_precue_all = sum(psth_cue(AllTrialsFinal,epoch_analyze_loc),2);
        SpikeCount_precue_all = nansum(psth_cue(:,epoch_analyze_loc),2);
        %fr_precue_all = nanmean(npsth_cue(:,epoch_analyze_loc),2);
        switch stim_type
            case 'cue'
                fr_precue_all =zscore(nanmean(psth_cue(:,epoch_analyze_loc),2));
                fr_precue_all_raw=nansum(psth_cue(:,epoch_analyze_loc),2);
                fr_precue_all_ms=nanmean(psth_cue(:,epoch_analyze_loc),2)-nanmean(nanmean(psth_cue(:,epoch_analyze_loc),2));
                fr_precue_all_raw_2=nanmean(psth_cue(:,epoch_analyze_loc),2);
                
                fr_postcue_all =zscore(nanmean(psth_cue(:,epoch_analyze_postcue_loc),2));
                fr_postcue_all_ms=nanmean(psth_cue(:,epoch_analyze_postcue_loc),2)-nanmean(nanmean(psth_cue(:,epoch_analyze_postcue_loc),2));
                fr_postcue_all_raw=nansum(psth_cue(:,epoch_analyze_postcue_loc),2);
                fr_postcue_all_raw_2=nanmean(psth_cue(:,epoch_analyze_postcue_loc),2);
            case 'response'
                fr_precue_all =zscore(nanmean(psth_response(:,epoch_analyze_loc),2));
                fr_precue_all_raw=nansum(psth_response(:,epoch_analyze_loc),2);
                fr_precue_all_raw_2=nanmean(psth_response(:,epoch_analyze_loc),2);
                
                fr_postcue_all =zscore(nanmean(psth_response(:,epoch_analyze_postcue_loc),2));
                fr_postcue_all_raw=nansum(psth_response(:,epoch_analyze_postcue_loc),2);
                fr_postcue_all_raw_2=nanmean(psth_response(:,epoch_analyze_postcue_loc),2);
                fr_postcue_all_ms=nanmean(psth_response(:,epoch_analyze_loc),2)-nanmean(nanmean(psth_response(:,epoch_analyze_loc),2));

                
            case 'both'
                
                fr_postcue_all =zscore(nanmean(psth_response(:,epoch_analyze_postcue_loc),2));
                fr_postcue_all =zscore(nanmean(psth_response(:,epoch_analyze_postcue_loc),2));
                fr_postcue_all_raw=nansum(psth_response(:,epoch_analyze_postcue_loc),2);
                fr_postcue_all_raw_2=nanmean(psth_response(:,epoch_analyze_postcue_loc),2);
                fr_postcue_all_ms=nanmean(psth_response(:,epoch_analyze_loc),2)-nanmean(nanmean(psth_response(:,epoch_analyze_loc),2));

                
                fr_precue_all =zscore(nanmean(psth_cue(:,epoch_analyze_loc),2));
                fr_precue_all_raw=nansum(psth_cue(:,epoch_analyze_loc),2);
                fr_precue_all_raw_2=nanmean(psth_cue(:,epoch_analyze_loc),2);
                fr_precue_all_ms=nanmean(psth_cue(:,epoch_analyze_loc),2)-nanmean(nanmean(psth_cue(:,epoch_analyze_loc),2));

        end
        %SpikeCount_precue_all_allCells{c,ii}=SpikeCount_precue_all;
        %precue_all_ds_PredFR=table(SpikeCount_precue_all, prevEriksenTrialsAll, prevSimonTrialsAll, lastRTall, RTall);
        %precue_all_ds_PredFR=table(SpikeCount_precue_all, prevEriksenTrials, prevSimonTrials, lastRTall, RTall);
        fr_precue_all_mean_epoch(c)=nanmean(nansum(psth_cue(:,800:1000),2))/2;
        fr_precure_all_all{c}=fr_precue_all;
        fr_precure_all_all_raw{c}=fr_precue_all_raw;
        %fr_precue_Simon=fr_precue_all.*SimonTrialsFinal;
        fr_precue_Simon=fr_precue_all.*simonTrials;
        fr_precue_Simon_ms=fr_precue_all_ms.*simonTrials;

        fr_postcue_Simon=fr_postcue_all.*simonTrials;
        fr_postcue_Simon_ms=fr_postcue_all_ms.*simonTrials;

        fr_postcue_NextSimon=fr_postcue_all.*NextSimonTrials;
        fr_precue_SimonRaw=fr_precue_all_raw_2.*simonTrials;
        %fr_precue_Eriksen=fr_precue_all.*EriksenTrialsFinal;
        fr_precue_Eriksen=fr_precue_all.*eriksenTrials;
        fr_precue_Eriksen_ms=fr_precue_all_ms.*eriksenTrials;

        fr_postcue_Eriksen=fr_postcue_all.*eriksenTrials;
        fr_postcue_Eriksen_ms=fr_postcue_all_ms.*eriksenTrials;

        fr_postcue_NextEriksen=fr_postcue_all.*NextEriksenTrials;
        
        fr_precue_EriksenRaw=fr_precue_all_raw_2.*eriksenTrials;
        
        % fr_precue_Either=fr_precue_all.*(SimonTrialsFinal | EriksenTrialsFinal);
        fr_precue_Either=fr_precue_all.*(simonTrials | eriksenTrials);
        fr_precue_Either_ms=fr_precue_all_ms.*(simonTrials | eriksenTrials);


        fr_postcue_Either_ms=fr_postcue_all_ms.*(simonTrials | eriksenTrials);

        fr_precue_None=fr_precue_all.*(~SimonTrialsFinal & ~EriksenTrialsFinal);
        fr_precue_SimonExc=fr_precue_all.*simonTrialsExc;
        fr_precue_EriksenExc=fr_precue_all.*eriksenTrialsExc;
        fr_precue_Both=fr_precue_all.*bothTrials;
        fr_precue_Both_quad=(fr_precue_all.*BothTrialsFinal).*(fr_precue_all.*BothTrialsFinal).*sign(fr_precue_all);
        fr_precue_preConf=fr_precue_all(prevConflict);
        
        fr_precue_NoPreConf=fr_precue_all(~prevConflictBinary);
        
        %[h_preConf,p_preConf]=ttest2(fr_precue_preConf(2:end),fr_precue_NoPreConf(2:end));
        [p_preConf,h_preConf]=ranksum(fr_precue_preConf(2:end),fr_precue_NoPreConf(2:end));
        h_preConf_All(c,ii)=h_preConf;
        p_preConf_All(c,ii)=p_preConf;
        
        % precue_all_ds_PredFR_all{c,ii}=precue_all_ds_PredFR;
        
        %precue_all_ds_PredFR_one=table(SpikeCount_precue_all, RTAllSimon, RTAllEriksen, RTAllNone, RTAllSimonExc, RTAllEriksenExc, RTAllBoth,prevEriksenTrials, prevSimonTrials, prevConflictBinary, lastRTall, RTall, RTAllEither, prevCorrect);
        % precue_all_ds_PredFR_one_all{c,ii}=precue_all_ds_PredFR_one;
        
        
                       RTeither=RTall.*(simonTrials | eriksenTrials);
        fr_precue_Either=fr_precue_all.*(simonTrials | eriksenTrials);
                 RTsimon=RTall.*(simonTrials);
        fr_precue_simon=fr_precue_all.*(simonTrials);
                        RTeriksen=RTall.*(eriksenTrials);
        fr_precue_eriksen=fr_precue_all.*(eriksenTrials);
        
        %subjVec=repmat(subj,numTrials,1);
        
        precue_all_ds_PredRT_one=table(RTeither,RTsimon,fr_precue_simon,fr_precue_Simon_ms,fr_precue_Either_ms,fr_postcue_Either_ms,RTeriksen,fr_precue_eriksen,fr_precue_Eriksen_ms,fr_postcue_NextEriksen,fr_postcue_NextSimon,NextEriksenTrials,NextSimonTrials,logRT,logLastRT,zLogRT,zLogLastRT,RTall, RTdiff, nextRTall, fr_postcue_Eriksen,fr_postcue_Simon, fr_postcue_all,fr_postcue_Eriksen_ms,fr_postcue_Simon_ms, fr_postcue_all_ms,fr_precue_all_raw,fr_postcue_all_raw, fr_precue_all_raw_2, EitherTrials, bothTrials, eriksenTrialsExc, prevIncorrect,simonTrialsExc, fr_precue_SimonExc, fr_precue_EriksenExc, fr_precue_all, fr_precue_all_ms,fr_precue_Simon, ConflictBinFinal, CongTrialsBin, InCongTrialsBin, prevConflict, NoConflictBin, CorrectTrials, IncorrectTrials, noPrevConflict, ConflictBin, simonTrials, eriksenTrials, SimonTrialsFinal, EriksenTrialsFinal, fr_precue_Eriksen, fr_precue_Either, fr_precue_None, fr_precue_SimonExc, fr_precue_EriksenExc, fr_precue_Both,fr_precue_Both_quad, prevEriksenTrials, prevSimonTrials, prevConflictBinary, lastRTall, RTAllEither, prevCorrect,prevConflictCatFinal,prevConflictCat,ConflictCatFinal,ConflictCat,prevConflictBinary);
       
        precue_all_ds_PredRT_AllCells{c}=precue_all_ds_PredRT_one;
        
        
        
        % fastTrialsErik = RTAllEriksen <= quantile(RTAllEriksen,0.4);
        % fastTrialsSimon = RTAllSimon <= quantile(RTAllSimon,0.4);
        
        fastTrials = RTall <= quantile(RTall,0.5);
        slowTrials = RTall > quantile(RTall,0.5);
        
        
        RTAllEriksen(RTAllEriksen==0)=NaN;
        RTAllSimon(RTAllSimon==0)=NaN;
        RTAllNoConf(RTAllNoConf==0)=NaN;
        fastTrialsErik = RTAllEriksen <= quantile(RTAllEriksen,0.5);
        fastTrialsSimon = RTAllSimon <= quantile(RTAllSimon,0.5);
        fastTrialsNoConf = RTAllNoConf <= quantile(RTAllNoConf,0.5);

        slowTrialsNoConf = RTAllNoConf > quantile(RTAllNoConf,0.5);
        slowTrialsErik = RTAllEriksen > quantile(RTAllEriksen,0.5);
        slowTrialsSimon = RTAllSimon > quantile(RTAllSimon,0.5);
        
        fastTrialsErikExc = RTAllEriksenExc <= quantile(RTAllEriksenExc,0.5);
        fastTrialsSimonExc = RTAllSimonExc <= quantile(RTAllSimonExc,0.5);
        
        slowTrialsErikExc = RTAllEriksenExc > quantile(RTAllEriksenExc,0.5);
        slowTrialsSimonExc = RTAllSimonExc > quantile(RTAllSimonExc,0.5);
        
        %  slowTrialsErik = RTAllEriksen > quantile(RTAllEriksen,0.6);
        % slowTrialsSimon = RTAllSimon > quantile(RTAllSimon,0.6);
        
        fastTrialsAnyConf = RTAllEither <= quantile(RTAllEither,0.5);
        fastTrialsNoConf = RTAllNone <= quantile(RTAllNone,0.5);
        
        
        %slowTrialsAnyConf = RTAllEither > quantile(RTAllEither,0.66);
        % slowTrialsNoConf = RTAllNone > quantile(RTAllNone,0.66);
        
        slowTrialsAnyConf = RTAllEither > quantile(RTAllEither,0.5);
        slowTrialsNoConf = RTAllNone > quantile(RTAllNone,0.5);
        
        RTAllCong(RTAllCong==0)=NaN;
        RTAllInCong(RTAllInCong==0)=NaN;
        fastTrialsCong = RTAllCong <= quantile(RTAllCong,0.5);
        fastTrialsInCong = RTAllInCong <= quantile(RTAllInCong,0.5);
        
        
        %slowTrialsAnyConf = RTAllEither > quantile(RTAllEither,0.66);
        % slowTrialsNoConf = RTAllNone > quantile(RTAllNone,0.66);
        
        slowTrialsCong = RTAllEither > quantile(RTAllCong,0.5);
        slowTrialsInCong = RTAllNone > quantile(RTAllInCong,0.5);
        
        [p_FastSlow_fix_RS(c),h_FastSlow_fix_RS(c),stats_FastSlow_fix_RS{c}]=ranksum(fr_precue_all(fastTrials),fr_precue_all(slowTrials));
        [h_FastSlow_fix(c),p_FastSlow_fix(c),~,stats_FastSlow_fix{c}]=ttest2(fr_precue_all(fastTrials),fr_precue_all(slowTrials));
        
        fr_prefix_fast{c}=fr_precue_all_raw(fastTrials);
        fr_prefix_slow{c}=fr_precue_all_raw(slowTrials);
        
        
        fr_precue_all_ms_ErikFast{c} = fr_precue_all_ms(fastTrialsErik);
           fr_precue_all_ms_ErikSlow{c} = fr_precue_all_ms(slowTrialsErik);
       fr_precue_all_ms_SimonFast{c} = fr_precue_all_ms(fastTrialsSimon);
           fr_precue_all_ms_SimonSlow{c} = fr_precue_all_ms(slowTrialsSimon);

        [h_ErikFastSlow(c),p_ErikFastSlow(c),~,stats_ErikFastSlow{c}]=ttest2(fr_precue_all(fastTrialsErik),fr_precue_all(slowTrialsErik));
        [h_SimonFastSlow(c),p_SimonFastSlow(c),~,stats_SimonFastSlow{c}]=ttest2(fr_precue_all(fastTrialsSimon),fr_precue_all(slowTrialsSimon));
        [p_ErikFastSlow_RS(c),h_ErikFastSlow_RS(c),stats_ErikFastSlow_RS{c}]=ranksum(nanmean(psth_cue(fastTrialsErik,epoch_analyze_loc),2),nanmean(psth_cue(slowTrialsErik,epoch_analyze_loc),2));
        [p_SimonFastSlow_RS(c),h_ErikFastSlow_RS(c),stats_ErikFastSlow_RS{c}]=ranksum(nanmean(psth_cue(fastTrialsSimon,epoch_analyze_loc),2),nanmean(psth_cue(slowTrialsSimon,epoch_analyze_loc),2));
        
        [h_ErikSimonFast(c),p_ErikSimonFast(c),~,stats_ErikSimonFast{c}]=ttest2(fr_precue_all(fastTrialsErik),fr_precue_all(fastTrialsSimon));
        [h_ErikSimonSlow(c),p_ErikSimonSlow(c),~,stats_ErikSimonSlow{c}]=ttest2(fr_precue_all(slowTrialsErik),fr_precue_all(slowTrialsSimon));



        [p_ErikFastSlow_RS2(c),h_ErikFastSlow_RS2(c),stats_ErikFastSlow_RS2{c}]=ranksum(fr_precue_all(fastTrialsErik),fr_precue_all(slowTrialsErik));
        [p_SimonFastSlow_R2(c),h_SimonFastSlow_RS2(c),stats_SimonFastSlow_RS2{c}]=ranksum(fr_precue_all(fastTrialsSimon),fr_precue_all(slowTrialsSimon));
        
        
        [h_ErikFastSlowPost(c),p_ErikFastSlowPost(c)]=ttest2(fr_postcue_all(fastTrialsErik),fr_postcue_all(slowTrialsErik,:));
        [h_SimonFastSlowPost(c),p_SimonFastSlowPost(c)]=ttest2(fr_postcue_all(fastTrialsSimon),fr_postcue_all(slowTrialsSimon));
        [p_ErikFastSlow_RS2_Post(c),h_ErikFastSlow_RS2_Post(c)]=ranksum(fr_postcue_all(fastTrialsErik),fr_postcue_all(slowTrialsErik));
        [p_SimonFastSlow_R2_Post(c),h_SimonFastSlow_RS2_Post(c)]=ranksum(fr_postcue_all(fastTrialsSimon),fr_postcue_all(slowTrialsSimon));
        
        
        try
            [p_FR_prevCorvIncor_RS(c),h_FR_prevCorvIncor_RS(c),FR_prevCorvIncor_RS_stats{c}]=ranksum(fr_precue_all(logical(prevCorrect(2:end))),fr_precue_all(logical(prevIncorrect(2:end))));
            [h_FR_prevCorvIncor(c),p_FR_prevCorvIncor(c)]=ttest2(fr_precue_all(logical(prevCorrect(2:end))),fr_precue_all(logical(prevIncorrect(2:end))));
            [p_RT_prevCorvIncor_RS(c),h_RT_prevCorvIncor_RS(c),RT_prevCorvIncor_RS_stats{c}]=ranksum(RTall(logical(prevCorrect(2:end))),RTall(logical(prevIncorrect(2:end))));
            [p_RTdiff_prevCorvIncor_RS(c),h_RTdiff_prevCorvIncor_RS(c),RTdiff_prevCorvIncor_RS_stats{c}]=ranksum(RTdiff(logical(prevCorrect(2:end))),RTdiff(logical(prevIncorrect(2:end))));
        catch
            [p_FR_prevCorvIncor_RS(c),h_FR_prevCorvIncor_RS(c)]=deal(NaN,0);
            [h_FR_prevCorvIncor(c),p_FR_prevCorvIncor(c)]=deal(0,NaN);
            [p_RT_prevCorvIncor_RS(c),h_RT_prevCorvIncor_RS(c)]=deal(NaN,0);
            [p_RTdiff_prevCorvIncor_RS(c),h_RTdiff_prevCorvIncor_RS(c)]=deal(NaN,0);
        end
        
        
        
        [p_ErikFastSlow_raw_RS(c),h_ErikFastSlow_raw_RS(c)]=ranksum(fr_precue_all_raw_2(fastTrialsErik),fr_precue_all_raw_2(slowTrialsErik));
        [p_SimonFastSlow_raw_RS(c),h_SimonFastSlow_raw_RS(c)]=ranksum(fr_precue_all_raw_2(fastTrialsSimon),fr_precue_all_raw_2(slowTrialsSimon));
        
        [h_ErikFastSlowExc(c),p_ErikFastSlowExc(c)]=ttest2(fr_precue_all(fastTrialsErikExc),fr_precue_all(slowTrialsErikExc));
        [h_SimonFastSlowExc(c),p_SimonFastSlowExc(c)]=ttest2(fr_precue_all(fastTrialsSimonExc),fr_precue_all(slowTrialsSimonExc));
        [p_ErikFastSlowExc_RS(c),h_ErikFastSlowExc_RS(c)]=ranksum(fr_precue_all(fastTrialsErikExc),fr_precue_all(slowTrialsErikExc));
        [p_SimonFastSlowExc_RS(c),h_SimonFastSlowExc_RS(c)]=ranksum(fr_precue_all(fastTrialsSimonExc),fr_precue_all(slowTrialsSimonExc));
        
        
        [h_AnyConfFastSlow(c),p_AnyConfFastSlow(c)]=ttest2(fr_precue_all(fastTrialsAnyConf),fr_precue_all(slowTrialsAnyConf));
        [h_NoConfFastSlow(c),p_NoConfFastSlow(c)]=ttest2(fr_precue_all(fastTrialsNoConf),fr_precue_all(slowTrialsNoConf));
        [h_CongFastSlow(c),p_CongFastSlow(c)]=ttest2(fr_precue_all(fastTrialsCong),fr_precue_all(slowTrialsCong));
        [h_InCongFastSlow(c),p_InCongFastSlow(c)]=ttest2(fr_precue_all(fastTrialsInCong),fr_precue_all(slowTrialsInCong));
        
        fr_precue_raw=nanmean(psth_cue(:,epoch_analyze_loc),2);
        fr_ErikFast{c,:}=fr_precue_raw(fastTrialsErik);
        fr_ErikSlow{c,:}=fr_precue_raw(slowTrialsErik);
        fr_SimonFast{c,:}=fr_precue_raw(fastTrialsSimon);
        fr_SimonSlow{c,:}=fr_precue_raw(slowTrialsSimon);
        
        zpsth_cue=zscore(psth_cue);
        
        
        psth_Fast_loc=psth_fixation(fastTrials,:);
        RTs_Fast_loc=RT(fastTrials);
        psth_Fast_fix{c}(:,:)=psth_Fast_loc;
        % zpsth_Fast_fix{c}(:,:)=zpsth_Fast_loc;
        RTs_Fast{c}(:,:)=RTs_Fast_loc;
        
        psth_Slow_loc=psth_fixation(slowTrials,:);
        RTs_Slow_loc=RT(slowTrials);
        psth_Slow_fix{c}(:,:)=psth_Slow_loc;
        %zpsth_Slow_fix{c}(:,:)=zpsth_Slow_loc;
        RTs_Slow{c}(:,:)=RTs_Slow_loc;
        
        
        psth_NoConffast_loc=psth_cue(fastTrialsNoConf,:);
        zpsth_NoConffast_loc=zpsth_cue(fastTrialsNoConf,:);
        RTs_NoConffast_loc=RT(fastTrialsNoConf);
        [RTs_NoConffast_loc,I]=sort(RTs_NoConffast_loc);
        psth_NoConffast_loc=psth_NoConffast_loc(I,:);
        zpsth_NoConffast_loc=zpsth_NoConffast_loc(I,:);
        psth_NoConfFast{c}(:,:)=psth_NoConffast_loc;
        zpsth_NoConfFast{c}(:,:)=zpsth_NoConffast_loc;
        RTs_NoConfFast{c}(:,:)=RTs_NoConffast_loc;


        psth_NoConfslow_loc=psth_cue(slowTrialsNoConf,:);
        zpsth_NoConfslow_loc=zpsth_cue(slowTrialsNoConf,:);
        RTs_NoConfslow_loc=RT(slowTrialsNoConf);
        [RTs_NoConfslow_loc,I]=sort(RTs_NoConfslow_loc);
        psth_NoConfslow_loc=psth_NoConfslow_loc(I,:);
        zpsth_NoConfslow_loc=zpsth_NoConfslow_loc(I,:);
        psth_NoConfSlow{c}(:,:)=psth_NoConfslow_loc;
        zpsth_NoConfSlow{c}(:,:)=zpsth_NoConfslow_loc;
        RTs_NoConfSlow{c}(:,:)=RTs_NoConfslow_loc;
        
        psth_ErikFast_loc=psth_cue(fastTrialsErik,:);
        zpsth_ErikFast_loc=zpsth_cue(fastTrialsErik,:);
        RTs_ErikFast_loc=RT(fastTrialsErik);
        [RTs_ErikFast_loc,I]=sort(RTs_ErikFast_loc);
        psth_ErikFast_loc=psth_ErikFast_loc(I,:);
        zpsth_ErikFast_loc=zpsth_ErikFast_loc(I,:);
        psth_ErikFast{c}(:,:)=psth_ErikFast_loc;
        zpsth_ErikFast{c}(:,:)=zpsth_ErikFast_loc;
        RTs_ErikFast{c}(:,:)=RTs_ErikFast_loc;
        
        psth_SimonFast_loc=psth_cue(fastTrialsSimon,:);
        zpsth_SimonFast_loc=zpsth_cue(fastTrialsSimon,:);
        RTs_SimonFast_loc=RT(fastTrialsSimon);
        [RTs_SimonFast_loc,I]=sort(RTs_SimonFast_loc);
        psth_SimonFast_loc=psth_SimonFast_loc(I,:);
        zpsth_SimonFast_loc=zpsth_SimonFast_loc(I,:);
        zpsth_SimonFast{c}(:,:)=zpsth_SimonFast_loc;
        psth_SimonFast{c}(:,:)=psth_SimonFast_loc;
        RTs_SimonFast{c}(:,:)=RTs_SimonFast_loc;
        
        psth_SimonSlow_loc=psth_cue(slowTrialsSimon,:);
        zpsth_SimonSlow_loc=zpsth_cue(slowTrialsSimon,:);
        RTs_SimonSlow_loc=RT(slowTrialsSimon);
        [RTs_SimonSlow_loc,I]=sort(RTs_SimonSlow_loc);
        psth_SimonSlow_loc=psth_SimonSlow_loc(I,:);
        zpsth_SimonSlow_loc=zpsth_SimonSlow_loc(I,:);
        psth_SimonSlow{c}(:,:)=psth_SimonSlow_loc;
        zpsth_SimonSlow{c}(:,:)=zpsth_SimonSlow_loc;
        RTs_SimonSlow{c}(:,:)=RTs_SimonSlow_loc;
        
        
        psth_ErikSlow_loc=psth_cue(slowTrialsErik,:);
        zpsth_ErikSlow_loc=zpsth_cue(slowTrialsErik,:);
        RTs_ErikSlow_loc=RT(slowTrialsErik);
        [RTs_ErikSlow_loc,I]=sort(RTs_ErikSlow_loc);
        psth_ErikSlow_loc=psth_ErikSlow_loc(I,:);
        zpsth_ErikSlow_loc=zpsth_ErikSlow_loc(I,:);
        psth_ErikSlow{c}(:,:)=psth_ErikSlow_loc;
        zpsth_ErikSlow{c}(:,:)=zpsth_ErikSlow_loc;
        RTs_ErikSlow{c}(:,:)=RTs_ErikSlow_loc;
        %
        %                 if c == 119
        %                     keyboard
        %                 end
        %
        %         psth_ErikSlow{c}(:,:)=psth_cue(slowTrialsErik,epoch_analyze_loc);
        %         psth_SimonFast{c}(:,:)=psth_cue(fastTrialsSimon,epoch_analyze_loc);
        %         psth_SimonSlow{c}(:,:)=psth_cue(slowTrialsSimon,epoch_analyze_loc);
        %
        %Fix the dimensions on this
        RTAllSimon_impute=(RTall.*SimonTrialsFinal + meanRTEriksen*EriksenTrialsFinalExc + meanRTNone*NoneTrialsFinal);
        RTAllEriksen_impute=(RTall.*EriksenTrialsFinal + meanRTSimon*SimonTrialsFinalExc + meanRTNone*NoneTrialsFinal);
        RTAllNone_impute=(RTall.*NoneTrialsFinal + meanRTSimon*SimonTrialsFinalExc + meanRTEriksen*EriksenTrialsFinalExc);
        RTAllEither_impute=(RTall.*(SimonTrialsFinal | EriksenTrialsFinal) + meanRTNone*NoneTrialsFinal);
        %         RTAllCong_impute=(RTall.*NoneTrialsFinal + meanRTSimon*SimonTrialsFinalExc + meanRTEriksen*EriksenTrialsFinalExc);
        %         RTAllEither_impute=(RTall.*(SimonTrialsFinal | EriksenTrialsFinal) + meanRTNone*NoneTrialsFinal);
        
        RTAllEriksenExc_impute=(RTall.*EriksenTrialsFinalExc + meanRTSimon*SimonTrialsFinalExc + meanRTNone*NoneTrialsFinal);
        RTAllSimonExc_impute=(RTall.*SimonTrialsFinalExc + meanRTEriksen*EriksenTrialsFinalExc + meanRTNone*NoneTrialsFinal);
        RTAllBoth_impute=(RTall.*BothTrialsFinal + meanRTEriksen*EriksenTrialsFinalExc + meanRTSimon*SimonTrialsFinalExc + meanRTNone*NoneTrialsFinal);
        %         RTAllSimon_impute=(RTall.*SimonTrialsFinal + (mean(RTall(curEriksenTrialsAll))+mean(RTall(curNoneTrialsAll))).*(~curSimonTrialsAll));
        %         RTAllEriksen_impute=(RTall.*curEriksenTrialsAll + (mean(RTall(curSimonTrialsAll))+mean(RTall(curNoneTrialsAll))).*(~curEriksenTrialsAll));
        %         RTAllNone_impute=(RTall.*curNoneTrialsAll + (mean(RTall(curEriksenTrialsAll))+mean(RTall(curSimonTrialsAll))).*(~curNoneTrialsAll));
        %         RTAllEither_impute=(RTall.*(curSimonTrialsAll | curEriksenTrialsAll) + mean(RTall(curNoneTrialsAll)).*~(curSimonTrialsAll | curEriksenTrialsAll));
        
        %precue_all_ds_PredFR_all_2=table(SpikeCount_precue_all, prevEriksenTrialsAll, prevSimonTrialsAll, prevConflictBinary, lastRTall,RTAllSimon_impute, RTAllEriksen_impute, RTAllNone_impute,RTAllEither_impute);
        % precue_all_ds_PredFR_2=table(SpikeCount_precue_all, prevEriksenTrials, prevSimonTrials, prevConflictBinary, lastRTall,RTAllSimon_impute, RTAllEriksen_impute, RTAllNone_impute,RTAllEither_impute, RTAllEriksenExc_impute, RTAllSimonExc_impute, RTAllBoth_impute, prevCorrect);
        %precue_all_ds_PredFR_2_all{c,ii}=precue_all_ds_PredFR_2;
        precue_all_ds_PredRT_all{c,ii}=precue_all_ds_PredRT_one;
        
        
        if correct_flag==1
            [phat,pci]=gamfit(RTorig,0.05);
            temppdf=gampdf(RTorig,phat(1),phat(2));
            %exclude_loc=[find(~CorrectAll | SpikeCount_precue_all>0)];
            %[B,TF]=rmoutliers(RTall,'ThresholdFactor',3);
            exclude_loc=[find(~CorrectTrials | isnan(lastConflict) | isnan(lastRT) | temppdf<exclude_prob)];
            
            exclude_loc_incor=[find(isnan(lastConflict) | isnan(lastRT))];
            
            exclude_loc_next=[find(~CorrectTrials | isnan(lastConflict) | isnan(lastRT) | (temppdf + 1)<exclude_prob)];
            
            %  exclude_loc=[find(~CorrectTrials |  ExcludeCellsRT  | isnan(lastConflict) | isnan(lastRT))];
            
            clear B
        else
            [phat,pci]=gamfit(RTorig,0.05);
            temppdf=gampdf(RTorig,phat(1),phat(2));
            %exclude_loc=[find(SpikeCount_precue_all>0)];
            %[B,TF]=rmoutliers(RTall,'ThresholdFactor',3);
            % exclude_loc=[find(ExcludeCellsRT  | isnan(lastConflict) | isnan(lastRT) | TF)];
            exclude_loc=[find(isnan(lastConflict) | isnan(lastRT) | temppdf<exclude_prob)];
            exclude_loc_incor=[find(isnan(lastConflict) | isnan(lastRT))];
            exclude_loc_next=[find(isnan(nextRTall))];
            
            % exclude_loc=[find(ExcludeCellsRT  | isnan(lastConflict) | isnan(lastRT))];
            
            %CorrectTrials = categorical(CorrectTrials);
            clear B
        end
        
        exclude_all{c,ii}=exclude_loc;

        
        %% Final Model %%
     
        if cor_uncor_flag == 1
            mdl_three_PredRT=fitglm(precue_all_ds_PredRT_one,'RTall ~ 1 + fr_precue_all_ms + fr_precue_Simon_ms + fr_precue_Eriksen_ms + eriksenTrials + simonTrials + prevEriksenTrials +prevSimonTrials + lastRTall + CorrectTrials','ResponseVar','RTall','PredictorVars',{'CorrectTrials','IncorrectTrials','fr_precue_all_ms','fr_precue_Eriksen_ms','fr_precue_Simon_ms','lastRTall','prevEriksenTrials','prevSimonTrials','simonTrials','eriksenTrials','CorrectTrials'},'Distribution',dist,'Link',linkfunc,'Exclude',exclude_loc);
            
            
        else
             mdl_three_PredRT=fitglm(precue_all_ds_PredRT_one,'RTall ~ 1 + fr_precue_all_ms + fr_precue_Simon_ms + fr_precue_Eriksen_ms + eriksenTrials + simonTrials + prevEriksenTrials +prevSimonTrials + lastRTall','ResponseVar','RTall','PredictorVars',{'CorrectTrials','IncorrectTrials','fr_precue_all_ms','fr_precue_Eriksen_ms','fr_precue_Simon_ms','lastRTall','prevEriksenTrials','prevSimonTrials','simonTrials','eriksenTrials','CorrectTrials'},'Distribution',dist,'Link',linkfunc,'Exclude',exclude_loc);

        end
        
        mdl_three_all{c}=mdl_three_PredRT;
        pvalues_mdl_three_PredRT_Eriksen(c,ii)=table2array(mdl_three_PredRT.Coefficients('fr_precue_Eriksen_ms','pValue'));
        pvalues_mdl_three_PredRT_Simon(c,ii)=table2array(mdl_three_PredRT.Coefficients('fr_precue_Simon_ms','pValue'));
        pvalues_mdl_three_PredRT_all(c,ii)=table2array(mdl_three_PredRT.Coefficients('fr_precue_all_ms','pValue'));
 
     
        AIC_mdl_three_PredRT(c,ii)=mdl_three_PredRT.ModelCriterion.AICc;
        AIC_mdl_three_PredRT_loc = mdl_three_PredRT.ModelCriterion.AICc;
        BIC_mdl_three_PredRT(c,ii)=mdl_three_PredRT.ModelCriterion.BIC;
        BIC_mdl_three_PredRT_loc = mdl_three_PredRT.ModelCriterion.BIC;
        d=devianceTest(mdl_three_PredRT);
        dev_test_mdl_three_PredRT{c,ii}=d;
        dev_pval_mdl_three_PredRT(c,ii)=table2array(d(2,4));
        %coefTest_pval_GLM_all_Cat_NoInt_PredFR(c,ii)=coefTest(mdl_GLM_all_Cat_NoInt_PredFR);
        %  Betas_mdl_three_PredRT_None(c,ii)=table2array(mdl_three_PredRT.Coefficients('fr_precue_None','Estimate'));
        Betas_mdl_three_PredRT_Eriksen_raw(c,ii)=table2array(mdl_three_PredRT.Coefficients('fr_precue_Eriksen_ms','Estimate'));
        Betas_mdl_three_PredRT_Simon_raw(c,ii)=table2array(mdl_three_PredRT.Coefficients('fr_precue_Simon_ms','Estimate'));
        Betas_mdl_three_PredRT_Simon=Betas_mdl_three_PredRT_Simon_raw;
        Betas_mdl_three_PredRT_Eriksen=Betas_mdl_three_PredRT_Eriksen_raw;
       
        sig_mdl_three_PredRT_Eriksen{c,ii}=(table2array(d(2,4))<0.05 & table2array(mdl_three_PredRT.Coefficients('fr_precue_Eriksen_ms','pValue'))<0.05);
        sig_mdl_three_PredRT_all{c,ii}=(table2array(d(2,4))<0.05 & table2array(mdl_three_PredRT.Coefficients('fr_precue_all_ms','pValue'))<0.05);
        sig_mdl_three_PredRT_all_noDev{c,ii}=(table2array(mdl_three_PredRT.Coefficients('fr_precue_all_ms','pValue'))<0.05);
        
        sig_mdl_three_PredRT_Simon{c,ii}=(table2array(d(2,4))<0.05 & table2array(mdl_three_PredRT.Coefficients('fr_precue_Simon_ms','pValue'))<0.05);
        sig_mdl_three_PredRT_Either{c,ii}=(table2array(d(2,4))<0.05 & table2array(mdl_three_PredRT.Coefficients('fr_precue_Eriksen_ms','pValue'))<0.05 | table2array(mdl_three_PredRT.Coefficients('fr_precue_Simon_ms','pValue'))<0.05);
        %sig_mdl_one_PredRT_Eriksen{c,ii}=(table2array(d(2,4))<0.05 & table2array(mdl_one_PredRT.Coefficients('fr_precue_Eriksen','pValue'))<0.05);
        sig_mdl_three_PredRT_Either_noDev{c,ii}=(table2array(mdl_three_PredRT.Coefficients('fr_precue_Eriksen_ms','pValue'))<0.05 | table2array(mdl_three_PredRT.Coefficients('fr_precue_Simon_ms','pValue'))<0.05);
        sig_mdl_three_PredRT_Either_noDev_Bonf{c,ii}=(table2array(mdl_three_PredRT.Coefficients('fr_precue_Eriksen_ms','pValue'))<Bonf_p_value | table2array(mdl_three_PredRT.Coefficients('fr_precue_Simon_ms','pValue'))<Bonf_p_value);
        
        sig_mdl_three_PredRT_Eriksen_noDev{c,ii}=(table2array(mdl_three_PredRT.Coefficients('fr_precue_Eriksen_ms','pValue'))<0.05);
        sig_mdl_three_PredRT_Simon_noDev{c,ii}=(table2array(mdl_three_PredRT.Coefficients('fr_precue_Simon_ms','pValue'))<0.05);
        
        
      

  
        %%post-cue
        if cor_uncor_flag == 1
                        mdl_three_PredRT_post=fitglm(precue_all_ds_PredRT_one,'RTall ~ 1 + fr_postcue_all_ms + fr_postcue_Simon_ms + fr_postcue_Eriksen_ms + eriksenTrials + simonTrials + prevEriksenTrials +prevSimonTrials + lastRTall + CorrectTrials','ResponseVar','RTall','PredictorVars',{'CorrectTrials','IncorrectTrials','fr_postcue_all_ms','fr_postcue_Eriksen_ms','fr_postcue_Simon_ms','lastRTall','prevEriksenTrials','prevSimonTrials','simonTrials','eriksenTrials','CorrectTrials'},'CategoricalVars','CorrectTrials','Distribution',dist,'Link',linkfunc,'Exclude',exclude_loc);



        else
                    mdl_three_PredRT_post=fitglm(precue_all_ds_PredRT_one,'RTall ~ 1 + fr_postcue_all_ms + fr_postcue_Simon_ms + fr_postcue_Eriksen_ms + eriksenTrials + simonTrials + prevEriksenTrials +prevSimonTrials + lastRTall','ResponseVar','RTall','PredictorVars',{'CorrectTrials','IncorrectTrials','fr_postcue_all_ms','fr_postcue_Eriksen_ms','fr_postcue_Simon_ms','lastRTall','prevEriksenTrials','prevSimonTrials','simonTrials','eriksenTrials','CorrectTrials'},'CategoricalVars','CorrectTrials','Distribution',dist,'Link',linkfunc,'Exclude',exclude_loc);

        end
        
        mdl_three_postcue_all{c}=mdl_three_PredRT_post;
  
        pvalues_mdl_three_PredRT_post_Eriksen(c,ii)=table2array(mdl_three_PredRT_post.Coefficients('fr_postcue_Eriksen_ms','pValue'));
        pvalues_mdl_three_PredRT_post_Simon(c,ii)=table2array(mdl_three_PredRT_post.Coefficients('fr_postcue_Simon_ms','pValue'));
        pvalues_mdl_three_PredRT_post_all(c,ii)=table2array(mdl_three_PredRT_post.Coefficients('fr_postcue_all_ms','pValue'));
        
        AIC_mdl_three_PredRT_post(c,ii)=mdl_three_PredRT_post.ModelCriterion.AICc;
        AIC_mdl_three_PredRT_post_loc = mdl_three_PredRT_post.ModelCriterion.AICc;
        BIC_mdl_three_PredRT_post(c,ii)=mdl_three_PredRT_post.ModelCriterion.BIC;
        BIC_mdl_three_PredRT_post_loc = mdl_three_PredRT_post.ModelCriterion.BIC;
        d=devianceTest(mdl_three_PredRT_post);
        dev_test_mdl_three_PredRT_post{c,ii}=d;
        dev_pval_mdl_three_PredRT_post(c,ii)=table2array(d(2,4));
    

  Betas_mdl_three_PredRT_post_Eriksen_raw(c,ii)=table2array(mdl_three_PredRT_post.Coefficients('fr_postcue_Eriksen_ms','Estimate'));
        Betas_mdl_three_PredRT_post_Simon_raw(c,ii)=table2array(mdl_three_PredRT_post.Coefficients('fr_postcue_Simon_ms','Estimate'));

      
        Betas_mdl_three_PredRT_post_Eriksen=Betas_mdl_three_PredRT_post_Eriksen_raw;
        Betas_mdl_three_PredRT_post_Simon=Betas_mdl_three_PredRT_post_Simon_raw;
    
        sig_mdl_three_PredRT_post_Eriksen{c,ii}=(table2array(d(2,4))<0.05 & table2array(mdl_three_PredRT_post.Coefficients('fr_postcue_Eriksen_ms','pValue'))<0.05);
        sig_mdl_three_PredRT_post_all{c,ii}=(table2array(d(2,4))<0.05 & table2array(mdl_three_PredRT_post.Coefficients('fr_postcue_all_ms','pValue'))<0.05);
        sig_mdl_three_PredRT_post_all_noDev{c,ii}=(table2array(mdl_three_PredRT_post.Coefficients('fr_postcue_all_ms','pValue'))<0.05);
        
        sig_mdl_three_PredRT_post_Simon{c,ii}=(table2array(d(2,4))<0.05 & table2array(mdl_three_PredRT_post.Coefficients('fr_postcue_Simon_ms','pValue'))<0.05);
        sig_mdl_three_PredRT_post_Either{c,ii}=(table2array(d(2,4))<0.05 & table2array(mdl_three_PredRT_post.Coefficients('fr_postcue_Eriksen_ms','pValue'))<0.05 | table2array(mdl_three_PredRT_post.Coefficients('fr_postcue_Simon_ms','pValue'))<0.05);
        %sig_mdl_one_PredRT_Eriksen{c,ii}=(table2array(d(2,4))<0.05 & table2array(mdl_one_PredRT.Coefficients('fr_postcue_Eriksen','pValue'))<0.05);
        sig_mdl_three_PredRT_post_Either_noDev{c,ii}=(table2array(mdl_three_PredRT_post.Coefficients('fr_postcue_Eriksen_ms','pValue'))<0.05 | table2array(mdl_three_PredRT_post.Coefficients('fr_postcue_Simon_ms','pValue'))<0.05);
        sig_mdl_three_PredRT_post_Eriksen_noDev{c,ii}=(table2array(mdl_three_PredRT_post.Coefficients('fr_postcue_Eriksen_ms','pValue'))<0.05);
        sig_mdl_three_PredRT_post_Simon_noDev{c,ii}=(table2array(mdl_three_PredRT_post.Coefficients('fr_postcue_Simon_ms','pValue'))<0.05);



    end
end


        %%Now that we have the model for each cell, compute the
        %%correlations and angles between cell regression coefficient
        %%vectors for all the conditions. 


    for ii=1:num_epochs


    
    
    
    %% NoNan  %%
    Inds=(~isnan(Betas_mdl_three_PredRT_Eriksen_raw(:,ii)) & ~isnan(Betas_mdl_three_PredRT_Simon_raw(:,ii)));
    A(:,1)=Betas_mdl_three_PredRT_Eriksen_raw(Inds,ii);
    A(:,2)=Betas_mdl_three_PredRT_Simon_raw(Inds,ii);
    %A=rmoutliers(A);
    x=A(:,1);
    y=A(:,2);
    x_cent=x-mean(x);
    y_cent=y-mean(y);
    Betas_mdl_three_PredRT_Eriksen_raw_NoNan{:,ii}=A(:,1);
    Betas_mdl_three_PredRT_Simon_raw_NoNan{:,ii} =A(:,2);
    
    
    clear A
    A(:,1)=Betas_mdl_three_PredRT_Eriksen(Inds,ii);
    A(:,2)=Betas_mdl_three_PredRT_Simon(Inds,ii);
    Betas_mdl_three_PredRT_Eriksen_NoNan{:,ii}=A(:,1);
    Betas_mdl_three_PredRT_Simon_NoNan{:,ii} =A(:,2);
    [A,ExcludedCells]=rmoutliers(A);
    Betas_mdl_three_PredRT_Eriksen_NoNan_NoOut{:,ii}=A(:,1);
    Betas_mdl_three_PredRT_Simon_NoNan_NoOut{:,ii} =A(:,2);
    IndsIncludedCells=find(~ExcludedCells);
    x=A(:,1);
    y=A(:,2);
    x_cent=x-mean(x);
    y_cent=y-mean(y);
    
    [rvalue_Spear_mdl_three_SimonVEriksen_NoNan(ii), pvalue_Spear_mdl_three_SimonVEriksen_NoNan(ii)] = corr(A(:,2), A(:,1), 'type', 'spearman');
    [rvalue_Pears_mdl_three_SimonVEriksen_NoNan(ii), pvalue_Pears_mdl_three_SimonVEriksen_NoNan(ii)] = corr(A(:,2), A(:,1), 'type', 'pearson');
    [rvalue_Kend_mdl_three_SimonVEriksen_NoNan(ii), pvalue_Kend_mdl_three_SimonVEriksen_NoNan(ii)] = corr(A(:,2), A(:,1), 'type', 'kendall');
    CI_Boot_rval_SimonVEriksen_NoNan_Spear.mdl_three(:,ii)= bootci(numPerms,{@corr_spear,A(:,1),A(:,2)},'alpha',0.05);
    
    
    [rvalue_Unsigned_Spear_mdl_three_SimonVEriksen_NoNan(ii), pvalue_Unsigned_Spear_Smdl_three_SimonVEriksen_NoNan(ii)] = corr(abs(A(:,2)), abs(A(:,1)), 'type', 'spearman');
    [rvalue_Unsigned_Pears_mdl_three_SimonVEriksen_NoNan(ii), pvalue_Unsigned_Pears_mdl_three_SimonVEriksen_NoNan(ii)] = corr(abs(A(:,2)), abs(A(:,1)), 'type', 'pearson');
    [rvalue_Unsigned_Kend_mdl_three_SimonVEriksen_NoNan(ii), pvalue_Unsigned_Kend_SimvErik.mdl_three_PredRT_NoNan(ii)] = corr(abs(A(:,2)), abs(A(:,1)), 'type', 'kendall');
    CI_Boot_Unsigned_rval_SimonVEriksen_NoNan_Spear.mdl_three(:,ii)= bootci(numPerms,{@corr_spear,abs(A(:,1)),abs(A(:,2))},'alpha',0.05);
    
    [rvalue_bootstat_precue_SimvErik_NoNan_spear,~]= bootstrp(numPerms,@corr_spear,x,y);
    [rvalue_bootstat_precue_SimvErik_NoNan_pears,~]= bootstrp(numPerms,@corr_pears,x,y);
    [rvalue_bootstat_precue_SimvErik_NoNan_kend,~]= bootstrp(numPerms,@corr_kend,x,y);
    
    
    A_cent=A-mean(A,1);
    Eriksen_NoNan(ii) = real(acosd(max(min(dot(A(:,1),A(:,2))/(norm(A(:,1))*norm(A(:,2))),1),-1)));
    CosAngCent_mdl_three_SimonVEriksen_NoNan(ii)=dot(A(:,1),A(:,2))/(norm(A(:,1))*norm(A(:,2)));
    CosAng_mdl_three_SimonVEriksen_NoNan(ii)=max(min(dot((A(:,1)-mean(A(:,1))),(A(:,2)-mean(A(:,2))))/(norm(A(:,1)-mean(A(:,1)))*norm(A(:,2)-mean(A(:,2)))),1),-1);
    Ang_Centered_mdl_three_SimonVEriksen_NoNan(ii) = real(acosd(max(min(dot((A(:,1)-mean(A(:,1))),(A(:,2)-mean(A(:,2))))/(norm(A(:,1)-mean(A(:,1)))*norm(A(:,2)-mean(A(:,2)))),1),-1)));
    ang_cent=calc_ang_cent(A);
        CI_Boot_AngCent_SimonVEriksen_NoNan_Solo(:,ii)= bootci(numPerms,{@calc_ang_cent,A},'alpha',0.05);

    
    [bootstat_ang_SimvErik_NoNan(:),~]=bootstrp(numPerms,@compute_Ang2,A);
    [bootstat_CosAng_SimvErik_NoNan(:),~]=bootstrp(numPerms,@compute_cos,A);
    %  [bootstat_CosAng_SimvErik_NoNan(:),~]=bootstrp(numPerms,@compute_cos_ang,x, y);
    % [bootstat_ang_cent_SimvErik_NoNan(:),~]=bootstrp(numPerms,@compute_angle,x_cent, y_cent);
    [bootstat_ang_cent_SimvErik_NoNan(:),~]=bootstrp(numPerms,@compute_Ang2,A_cent);
    
    
    clear A B A_cent B_cent
    
    
   
    Inds=(~isnan(Betas_mdl_three_PredRT_Eriksen(:,ii)) & ~isnan(Betas_mdl_three_PredRT_Simon(:,ii)) & ~isnan(Betas_mdl_three_PredRT_post_Eriksen(:,ii)) & ~isnan(Betas_mdl_three_PredRT_post_Simon(:,ii))) ;
    
    %% Pre Post Simon
    A(:,1)=Betas_mdl_three_PredRT_Eriksen(Inds,ii);
    A(:,2)=Betas_mdl_three_PredRT_Simon(Inds,ii);
    
    B(:,1)=Betas_mdl_three_PredRT_post_Eriksen(Inds,ii);
    B(:,2)=Betas_mdl_three_PredRT_post_Simon(Inds,ii);
    
    [A,I]=rmoutliers(A);
    B(I,:)=[];
    Betas_mdl_three_PredRT_Eriksen_pre_PrePost_NoNan=A(:,1);
    Betas_mdl_three_PredRT_Simon_pre_PrePost_NoNan=A(:,2);
    Betas_mdl_three_PredRT_Eriksen_post_PrePost_NoNan=B(:,1);
    Betas_mdl_three_PredRT_Simon_post_PrePost_NoNan=B(:,2);
    
    
    [rvalue_Spear_mdl_three_SimonPrePost_NoNan(ii), pvalue_Spear_mdl_three_SimonVEriksen_NoNan_post(ii)] = corr(A(:,2), B(:,2), 'type', 'spearman');
    [rvalue_Pears_mdl_three_SimonPrePost_NoNan(ii), pvalue_Pears_mdl_three_SimonVEriksen_NoNan_post(ii)] = corr(A(:,2), B(:,2), 'type', 'pearson');
    [rvalue_Kend_mdl_three_SimonPrePost_NoNan(ii), pvalue_Kend_mdl_three_SimonVEriksen_NoNan_post(ii)] = corr(A(:,2), B(:,2), 'type', 'kendall');
    CI_Boot_rval_SimonPrePost_NoNan_Spear.mdl_three(:,ii)= bootci(numPerms,{@corr_spear,A(:,2),B(:,2)},'alpha',0.05);
    
    [rvalue_Unsigned_Spear_mdl_three_SimonPrePost_NoNan(ii), pvalue_Unsigned_Spear_Sdl_three_SimonVEriksen_NoNan_post(ii)] = corr(abs(A(:,2)), abs(B(:,2)), 'type', 'spearman');
    [rvalue_Unsigned_Pears_mdl_three_SimonPrePost_NoNan(ii), pvalue_Unsigned_Pears_mdl_three_SimonVEriksen_NoNan_post(ii)] = corr(abs(A(:,2)), abs(B(:,2)), 'type', 'pearson');
    [rvalue_Unsigned_Kend_mdl_three_SimonPrePost_NoNan(ii), pvalue_Unsigned_Kend_SimvErik.mdl_three_PredRT_NoNan_post(ii)] = corr(abs(A(:,2)), abs(B(:,2)), 'type', 'kendall');
    CI_Boot_Unsigned_rval_SimonPrePost_NoNan_Spear.mdl_three(:,ii)= bootci(numPerms,{@corr_spear,abs(A(:,2)),abs(B(:,2))},'alpha',0.05);
    
    
    Ang_mdl_three_SimonPrePost_NoNan(ii) = real(acosd(max(min(dot(A(:,2),B(:,2))/(norm(B(:,2))*norm(A(:,2))),1),-1)));
    CosAng_mdl_three_SimonPrePost_NoNan(ii)=dot(A(:,2),B(:,2))/(norm(A(:,2))*norm(B(:,2)));
    Ang_Centered_mdl_three_SimonPrePost_NoNan(ii) = real(acosd(max(min(dot((A(:,2)-mean(A(:,2))),(B(:,2)-mean(B(:,2))))/(norm(A(:,2)-mean(A(:,2)))*norm(B(:,2)-mean(B(:,2)))),1),-1)));
    
    
    
    %% Pre Post Eriksen
    
    
    [rvalue_Spear_mdl_three_EriksenPrePost_NoNan(ii), pvalue_Spear_mdl_three_EriksenVEriksen_NoNan_post(ii)] = corr(A(:,1), B(:,1), 'type', 'spearman');
    [rvalue_Pears_mdl_three_EriksenPrePost_NoNan(ii), pvalue_Pears_mdl_three_EriksenVEriksen_NoNan_post(ii)] = corr(A(:,1), B(:,1), 'type', 'pearson');
    [rvalue_Kend_mdl_three_EriksenPrePost_NoNan(ii), pvalue_Kend_mdl_three_EriksenVEriksen_NoNan_post(ii)] = corr(A(:,1), B(:,1), 'type', 'kendall');
    CI_Boot_rval_EriksenPrePost_NoNan_Spear.mdl_three(:,ii)= bootci(numPerms,{@corr_spear,A(:,1),B(:,1)},'alpha',0.05);
    
    [rvalue_Unsigned_Spear_mdl_three_EriksenPrePost_NoNan(ii), pvalue_Unsigned_Spear_Sdl_three_EriksenVEriksen_NoNan_post(ii)] = corr(abs(A(:,1)), abs(B(:,1)), 'type', 'spearman');
    [rvalue_Unsigned_Pears_mdl_three_EriksenPrePost_NoNan(ii), pvalue_Unsigned_Pears_mdl_three_EriksenVEriksen_NoNan_post(ii)] = corr(abs(A(:,2)), abs(B(:,1)), 'type', 'pearson');
    [rvalue_Unsigned_Kend_mdl_three_EriksenPrePost_NoNan(ii), pvalue_Unsigned_Kend_SimvErik.mdl_three_PredRT_NoNan_post(ii)] = corr(abs(A(:,2)), abs(B(:,1)), 'type', 'kendall');
    CI_Boot_Unsigned_rval_EriksenPrePost_NoNan_Spear.mdl_three(:,ii)= bootci(numPerms,{@corr_spear,abs(A(:,2)),abs(B(:,1))},'alpha',0.05);
    
    
    
    Ang_mdl_three_EriksenPrePost_NoNan(ii) = real(acosd(max(min(dot(A(:,1),B(:,1))/(norm(B(:,1))*norm(A(:,1))),1),-1)));
    CosAng_mdl_three_EriksenPrePost_NoNan(ii)=dot(A(:,1),B(:,1))/(norm(A(:,1))*norm(B(:,1)));
    Ang_Centered_mdl_three_EriksenPrePost_NoNan(ii) = real(acosd(max(min(dot((A(:,1)-mean(A(:,1))),(B(:,1)-mean(B(:,1))))/(norm(A(:,1)-mean(A(:,1)))*norm(B(:,1)-mean(B(:,1)))),1),-1)));
    
    
    
    
   
    
    clear A B
    
    
  


    
    %%%%%%%
    clear A
  inds=find(cell2mat(sig_mdl_three_PredRT_Either_noDev(:,ii))==1);
     Betas_mdl_three_PredRT_Eriksen_SigConf =Betas_mdl_three_PredRT_Eriksen(inds,ii);
    Betas_mdl_three_PredRT_Simon_SigConf = Betas_mdl_three_PredRT_Simon(inds,ii);
    A(:,1) = Betas_mdl_three_PredRT_Eriksen_SigConf;
    A(:,2) = Betas_mdl_three_PredRT_Simon_SigConf;
    A=rmoutliers(A);
    Betas_mdl_three_PredRT_Eriksen_SigConf=A(:,1);
    Betas_mdl_three_PredRT_Simon_SigConf=A(:,2);
    x=A(:,1);
    y=A(:,2);
    x_cent=x-mean(x);
    y_cent=y-mean(y);
    [rvalue_Spear_mdl_three_SimonVEriksen_SigConf(ii), pvalue_Spear_mdl_three_SimonVEriksen_SigConf(ii)] = corr(Betas_mdl_three_PredRT_Eriksen_SigConf, Betas_mdl_three_PredRT_Simon_SigConf, 'type', 'spearman');
    [rvalue_Pears_mdl_three_SimonVEriksen_SigConf(ii), pvalue_Pears_mdl_three_SimonVEriksen_SigConf(ii)] = corr(Betas_mdl_three_PredRT_Eriksen_SigConf, Betas_mdl_three_PredRT_Simon_SigConf, 'type', 'pearson');
    [rvalue_Kend_mdl_three_SimonVEriksen_SigConf(ii), pvalue_Kend_mdl_three_SimonVEriksen_SigConf(ii)] = corr(Betas_mdl_three_PredRT_Eriksen_SigConf, Betas_mdl_three_PredRT_Simon_SigConf, 'type', 'kendall');
    CI_Boot_rval_SimonVEriksen_SigConf_Spear.mdl_three(:,ii)= bootci(numPerms,{@corr_spear,Betas_mdl_three_PredRT_Eriksen_SigConf,Betas_mdl_three_PredRT_Simon_SigConf},'alpha',0.05);
    CI_Boot_rval_SimonVEriksen_SigConf_Kend.mdl_three(:,ii)= bootci(numPerms,{@corr_kend,Betas_mdl_three_PredRT_Eriksen_SigConf,Betas_mdl_three_PredRT_Simon_SigConf},'alpha',0.05);
    [rvalue_bootstat_precue_SimvErik_SigConf_spear,~]= bootstrp(numPerms,@corr_spear,x,y);
    [rvalue_bootstat_precue_SimvErik_SigConf_pears,~]= bootstrp(numPerms,@corr_pears,x,y);
    [rvalue_bootstat_precue_SimvErik_SigConf_kend,~]= bootstrp(numPerms,@corr_kend,x,y);
    
    %%
    
    %A=rmoutliers(A);
    %     x=A(:,1);
    %     y=A(:,2);
    %     x_cent=x-mean(x);
    %     y_cent=y-mean(y);
    Ang_Centered_mdl_three_SimonVEriksen_SigConf(ii) = real(acosd(max(min(dot((x_cent-mean(x_cent)),(y_cent-mean(y_cent)))/(norm(x_cent-mean(x_cent))*norm(y_cent-mean(y_cent))),1),-1)));
    Ang_mdl_three_SimonVEriksen_SigConf(ii) = real(acosd(max(min(dot((A(:,1)),(A(:,2)))/(norm(A(:,1)-mean(A(:,1)))*norm(A(:,2)-mean(A(:,2)))),1),-1)));
    CosAng_Centered_mdl_three_SimonVEriksen_SigConf(ii) = max(min(dot((x_cent-mean(x_cent)),(y_cent-mean(y_cent)))/(norm(x_cent-mean(x_cent))*norm(y_cent-mean(y_cent))),1),-1);
    CosAng_mdl_three_SimonVEriksen_SigConf(ii) = max(min(dot((A(:,1)),(A(:,2)))/(norm(A(:,1)-mean(A(:,1)))*norm(A(:,2)-mean(A(:,2)))),1),-1);
    CI_Boot_rval_SimonVEriksen_SigConf_Spear_NoOut.mdl_three(:,ii)= bootci(numPerms,{@corr_spear,x,y},'alpha',0.05);
    CI_Boot_rval_SimonVEriksen_SigConf_Kend_NoOut.mdl_three(:,ii)= bootci(numPerms,{@corr_kend,x,y},'alpha',0.05);
            CI_Boot_AngCent_SimonVEriksen_SigConf_Solo(:,ii)= bootci(numPerms,{@calc_ang_cent,A},'alpha',0.05);

    
  %  CI_Boot_Ang_SimonVEriksen_SigConf.mdl_three(:,ii)= bootci(numPerms,{@compute_Ang2,A},'alpha',0.05);
    A_cent=A-mean(A,1);
    %CI_Boot_AngCent_SimonVEriksen_SigConf.mdl_three(:,ii)= bootci(numPerms,{@compute_Ang2,A_cent},'alpha',0.05);
    [bootstat_ang_SimvErik_SigConf(:),~]=bootstrp(numPerms,@compute_Ang2,A);
    [bootstat_ang_cent_SimvErik_SigConf(:),~]=bootstrp(numPerms,@compute_Ang2,A_cent);
    
    
    CI_Boot_CosAng_SimonVEriksen_SigConf.mdl_three(:,ii)= bootci(numPerms,{@compute_cos,A},'alpha',0.05);
    A_cent=A-mean(A,1);
    CI_Boot_CosAngCent_SimonVEriksen_SigConf.mdl_three(:,ii)= bootci(numPerms,{@compute_cos,A_cent},'alpha',0.05);
    [bootstat_cosang_SimvErik_SigConf(:),~]=bootstrp(numPerms,@compute_cos,A);
    [bootstat_cosang_cent_SimvErik_SigConf(:),~]=bootstrp(numPerms,@compute_cos,A_cent);
 
  
  
  %%%%%%%%
   
        clear A B
    %% Pre Post Simon
        Inds=(~isnan(Betas_mdl_three_PredRT_Eriksen(:,ii)) & ~isnan(Betas_mdl_three_PredRT_Simon(:,ii)) & ~isnan(Betas_mdl_three_PredRT_post_Eriksen(:,ii)) & ~isnan(Betas_mdl_three_PredRT_post_Simon(:,ii))) ;

    A(:,1)=Betas_mdl_three_PredRT_Eriksen(Inds,ii);
    A(:,2)=Betas_mdl_three_PredRT_Simon(Inds,ii);
    
    B(:,1)=Betas_mdl_three_PredRT_post_Eriksen(Inds,ii);
    B(:,2)=Betas_mdl_three_PredRT_post_Simon(Inds,ii);
    [A,I]=rmoutliers(A);
    B(I,:)=[];
    
    
    
 %%%%%% No Nan %%%%%%%       
 clear A B
        A(:,1)=cell2mat(Betas_mdl_three_PredRT_Eriksen_NoNan);
        A(:,2)=cell2mat(Betas_mdl_three_PredRT_Simon_NoNan);
        %  ExcInds=(abs(A(:,1)-median(A(:,1))) > MAD_Factor) | (abs(A(:,2)-median(A(:,2))) > MAD_Factor);
        % A(ExcInds,:)=[];
        A=rmoutliers(A);
        %  Thesh1=MAD_Factor*(-1/(sqrt(2)*erfcinv(3/2))*median(abs(A(:,1)-median(A(:,1)))));
        % Thesh2=MAD_Factor*(-1/(sqrt(2)*erfcinv(3/2))*median(abs(A(:,2)-median(A(:,2)))));
        %Inds=A(:,1)<Thesh1 & A(:,2)<Thesh2;
        
        x=A(:,1);
        y=A(:,2);
        % inds=find((abs(x)>(Beta_thresh*abs(std(x))+abs(mean(x))) | abs(y)>(Beta_thresh*abs(std(y))+abs(mean(y)))));
        
        %x(inds)=[];
        % y(inds)=[];
        x_cent=x-mean(x);
        y_cent=y-mean(y);

       Betas_NoNan_NoOut=A;
%        A(:,1)=B(sig_cells_eriksen_and_simon,1);
       % A(:,2)=B(sig_cells_eriksen_and_simon,2);
        %ExcInds=(abs(A(:,1)-median(A(:,1))) > MAD_Factor) | (abs(A(:,2)-median(A(:,2))) > MAD_Factor);
        % A(ExcInds,:)=[];
      %  A=rmoutliers(A);
       % x_sig=A(:,1);
       % y_sig=A(:,2);
        %inds=find((abs(x_sig)>(Beta_thresh*abs(std(x_sig))+abs(mean(x_sig))) | abs(y_sig)>(Beta_thresh*abs(std(y_sig))+abs(mean(y_sig)))));
        % x_sig(inds)=[];
        % y_sig(inds)=[];
       % x_sig_cent=x_sig-mean(x_sig);
        %y_sig_cent=y_sig-mean(y_sig);
        
        %clear A
       % A=horzcat(x,y);
     
        [bootstat_ang_SimvErik_NoNan(:),~]=bootstrp(numPerms,@compute_angle,x, y);
        [bootstat_CosAng_SimvErik_NoNan(:),~]=bootstrp(numPerms,@compute_cos_ang,x, y);
        [bootstat_ang_cent_SimvErik_NoNan(:),~]=bootstrp(numPerms,@compute_angle,x_cent, y_cent);
        [bootstat_CosAng_cent_SimvErik_NoNan(:),~]=bootstrp(numPerms,@compute_cos_ang,x_cent, y_cent);
        [rvalue_bootstat_precue_SimvErik_NoNan_spear_noout,~]= bootstrp(numPerms,@corr_spear,x,y);
        [rvalue_bootstat_precue_SimvErik_NoNan_pears_noout,~]= bootstrp(numPerms,@corr_pears,x,y);
        [rvalue_bootstat_precue_SimvErik_NoNan_kend_noout,~]= bootstrp(numPerms,@corr_kend,x,y);
%         CI_Boot_Ang_SimVErik_FDR_noout_cent(:)= bootci(numPerms,{@compute_angle,x_cent,y_cent},'alpha',0.05);
%         CI_Boot_Ang_SimVErik_FDR_noout(:)= bootci(numPerms,{@compute_angle,x,y},'alpha',0.05);
%         CI_Boot_CosAng_no_center_SimVErik_FDR_Spear_noout(:)= bootci(numPerms,{@compute_cos_ang,x,y},'alpha',0.05);
%         CI_Boot_CosAng_centered_SimVErik_FDR_Spear_noout(:)= bootci(numPerms,{@compute_cos_ang,x_cent,y_cent},'alpha',0.05);
%         CI_Boot_rval_SimVErik_FDR_Spear_noout(:)= bootci(numPerms,{@corr_spear,x,y},'alpha',0.05);
%         CI_Boot_rval_SimVErik_FDR_pears_noout(:)= bootci(numPerms,{@corr_pears,x,y},'alpha',0.05);
%         CI_Boot_rval_SimVErik_FDR_kend_noout(:)= bootci(numPerms,{@corr_kend,x,y},'alpha',0.05);
%         [rvalue_Spear_SimvErik.mdl_three_PredRT_FDR_NoOut, pvalue_Spear_SimvErik.mdl_three_PredRT_FDR_NoOut] = corr(A(:,2), A(:,1), 'type', 'spearman');
%         [rvalue_Pears_SimvErik.mdl_three_PredRT_FDR_NoOut, pvalue_Pears_SimvErik.mdl_three_PredRT_FDR_NoOut] = corr(A(:,2), A(:,1), 'type', 'pearson');
%         [rvalue_Kend_SimvErik.mdl_three_PredRT_FDR_NoOut, pvalue_Kend_SimvErik.mdl_three_PredRT_FDR_NoOut] = corr(A(:,2), A(:,1), 'type', 'kendall');
%         


        
        clear A B x y
   
        
        [rvalue_bootstat_precue_SimvErik_NoNan_spear,~]= bootstrp(numPerms,@corr_spear,cell2mat(Betas_mdl_three_PredRT_Simon_NoNan'),cell2mat(Betas_mdl_three_PredRT_Eriksen_NoNan'));
        [rvalue_bootstat_precue_SimvErik_NoNan_pears,~]= bootstrp(numPerms,@corr_pears,cell2mat(Betas_mdl_three_PredRT_Simon_NoNan'),cell2mat(Betas_mdl_three_PredRT_Eriksen_NoNan'));
        [rvalue_bootstat_precue_SimvErik_NoNan_kend,~]= bootstrp(numPerms,@corr_kend,cell2mat(Betas_mdl_three_PredRT_Simon_NoNan'),cell2mat(Betas_mdl_three_PredRT_Eriksen_NoNan'));
      
%%Set up vectors for angle and correlation permutation tests

        a=Betas_mdl_three_PredRT_Eriksen;
        a1=Betas_mdl_three_PredRT_post_Eriksen;
        
        
        b=Betas_mdl_three_PredRT_Simon;
        b1=Betas_mdl_three_PredRT_post_Simon;
        
        a_cent=a-mean(a);
        a1_cent=a1-mean(a1);
        b_cent=b-mean(b);
        b1_cent=b1-mean(b1);
        
        o=Betas_NoNan_NoOut(:,1);
        q=Betas_NoNan_NoOut(:,2);
        o_cent=o-mean(o);
        q_cent=q-mean(q);
        

        clear A B
      
        A(:,1) = Betas_mdl_three_PredRT_Eriksen_pre_PrePost_NoNan;
        A(:,2) = Betas_mdl_three_PredRT_Simon_pre_PrePost_NoNan;
        
        
        B(:,1) = Betas_mdl_three_PredRT_Eriksen_post_PrePost_NoNan;
        B(:,2) = Betas_mdl_three_PredRT_Simon_post_PrePost_NoNan;
        
        v5 = A(:,1);
        w5 = A(:,2);
        
        v6=B(:,1);
        w6=B(:,2);
     clear A B
        B(:,1) = Betas_mdl_three_PredRT_post_Eriksen;
        B(:,2) = Betas_mdl_three_PredRT_post_Simon;
        B=rmoutliers(B);
        v_post=B(:,1);
        w_post=B(:,1);
        clear A B
        
    
        
        B(:,1) = Betas_mdl_three_PredRT_post_Eriksen;
        B(:,2) = Betas_mdl_three_PredRT_post_Simon;
        B=rmoutliers(B);
        v_post=B(:,1);
        w_post=B(:,1);
        clear A B
        
        %%
        
        for p = 1:numPerms
            

            
            
           
            
            shuffBetas=o(randperm(length(o)));
            [rvals, ~]=corr(shuffBetas, q, 'type', 'spearman');
            [rvals_pears, ~]=corr(shuffBetas, q, 'type', 'pearson');
            [rvals_kend, ~]=corr(shuffBetas, q, 'type', 'kend');
            shuffBetas=o_cent(randperm(length(o_cent)));
            [Angs]=real(acosd(max(min(dot(shuffBetas,q)/(norm(shuffBetas)*norm(q)),1),-1)));
            [CosAngs]=max(min(dot(shuffBetas,q)/(norm(shuffBetas)*norm(q)),1),-1);
            
            [Angs_cent]=real(acosd(max(min(dot(shuffBetas,q_cent)/(norm(shuffBetas)*norm(q_cent)),1),-1)));
            [CosAngs_cent]=max(min(dot(shuffBetas,q_cent)/(norm(shuffBetas)*norm(q_cent)),1),-1);
            
            %size(rvalue_precue_perm_null_eriksen_simon)
            
            
            
            %size(rvalue_precue_perm_null_eriksen_simon)
            [rvalue_null_eriksen_simon_spears_mdl_three_PredRT_NoNan_NoOut(p)] = rvals;
            [rvalue_null_eriksen_simon_pears_mdl_three_PredRT_NoNan_NoOut(p)] = rvals_pears;
            [rvalue_null_eriksen_simon_kend_mdl_three_PredRT_NoNan_NoOut(p)] = rvals_kend;
            [Angs_cent_null_eriksen_simon_three_PredRT_NoNan_NoOut(p)] = Angs_cent;
            [Angs_null_eriksen_simon_three_PredRT_NoNan_NoOut(p)] = Angs;
            [CosAngs_cent_null_eriksen_simon_three_PredRT_NoNan_NoOut(p)] = CosAngs_cent;
            [CosAngs_null_eriksen_simon_three_PredRT_NoNan_NoOut(p)] = CosAngs;
        
            
            
     
            
            %%% Post All Cells %%
            shuffBetas=v_post(randperm(length(v_post)));
            %shuffBetas2=w_post(randperm(length(w_post)));
            [rvals, ~]=corr(shuffBetas, w_post, 'type', 'spearman');
            [rvals_pears, ~]=corr(shuffBetas, w_post, 'type', 'pearson');
            [rvals_kend, ~]=corr(shuffBetas, w_post, 'type', 'kend');
            %dotprod=dot(shuffBetas, I)/(norm(shuffBetas)*norm(I));
            %  eucdist=euc_dist(shuffBetas, I);
            
            % size(rvalue_precue_perm_null_eriksen_simon)
            [rvalue_null_eriksen_simon_spears_mdl_three_PredRT_AllCells_post(p)] = rvals;
            [rvalue_null_eriksen_simon_pears_mdl_three_PredRT_AllCells_post(p)] = rvals_pears;
            [rvalue_null_eriksen_simon_kend_mdl_three_PredRT_AllCells_post(p)] = rvals_kend;
            % shuffBetas=v2(randperm(length(v2)));
            %shuffBetas2=w_post2(randperm(length(w_post)));
            
            cosang=max(min(dot((w_post),(shuffBetas))/(norm(w_post)*norm(shuffBetas)),1),-1);
            cosang_cent=max(min(dot((w_post-mean(w_post)),(shuffBetas-mean(shuffBetas)))/(norm(w_post-mean(shuffBetas))*norm(shuffBetas-mean(shuffBetas))),1),-1);
            
            % cosang__null_eriksen_simon_mdl_three_PredRT_Dev(p)=dotprod;
            Ang_null_mdl_three_SimonVEriksen_AllCells_post(p) = real(acosd(max(min(dot((w_post),(shuffBetas))/(norm(w_post)*norm(shuffBetas)),1),-1)));
            AngCent_null_mdl_three_SimonVEriksen_AllCells_post(p) = real(acosd(max(min(dot((w_post-mean(w_post)),(shuffBetas-mean(shuffBetas)))/(norm(w_post-mean(w_post))*norm(shuffBetas-mean(shuffBetas))),1),-1)));
            CosAng_null_mdl_three_SimonVEriksen_AllCells_post(p)=cosang;
            CosAngCent_null_mdl_three_SimonVEriksen_AllCells_post(p)=cosang_cent;
            
            
            
            
            
         
            %% Simon Simon Post All Cells
            shuffBetas1=w5(randperm(length(w5)));
            % w6=w6(randperm(length(w6)));
            
            % w6=w(randperm(length(w)));
            [rvals, ~]=corr(shuffBetas1, w6, 'type', 'spearman');
            [rvals_pears, ~]=corr(shuffBetas1, w6, 'type', 'pearson');
            [rvals_kend, ~]=corr(shuffBetas1, w6, 'type', 'kend');
            %dotprod=dot(shuffBetas1, I)/(norm(shuffBetas1)*norm(I));
            %  eucdist=euc_dist(shuffBetas1, I);
            
            % size(rvalue_precue_perm_null_eriksen_simon)
            [rvalue_null_simon_simon_post_spears_mdl_three_PredRT_NoNan(p)] = rvals;
            [rvalue_null_simon_simon__postpears_mdl_three_PredRT_NoNan(p)] = rvals_pears;
            [rvalue_null_simon_simon_post_kend_mdl_three_PredRT_NoNan(p)] = rvals_kend;
            % shuffBetas1=v2(randperm(length(v2)));
            %shuffBetas12=w6(randperm(length(w6)));
            
            cosang=max(min(dot((w6),(shuffBetas1))/(norm(w6)*norm(shuffBetas1)),1),-1);
            cosang_cent=max(min(dot((w6-mean(w6)),(shuffBetas1-mean(shuffBetas1)))/(norm(w6-mean(w6))*norm(shuffBetas1-mean(shuffBetas1))),1),-1);
            
            % cosang__null_eriksen_simon_mdl_three_PredRT_Dev(p)=dotprod;
            Ang_null_mdl_three_SimonVSimon_NoNan(p) = real(acosd(max(min(dot((w6),(shuffBetas1))/(norm(w6)*norm(shuffBetas1)),1),-1)));
            AngCent_null_mdl_three_SimonVSimon_NoNan(p) = real(acosd(max(min(dot((w6-mean(w6)),(shuffBetas1-mean(shuffBetas1)))/(norm(w6-mean(w6))*norm(shuffBetas1-mean(shuffBetas1))),1),-1)));
            CosAng_null_mdl_three_SimonVSimon_NoNan(p)=cosang;
            CosAngCent_null_mdl_three_SimonVSimon_NoNan(p)=cosang_cent;
            
            
            
            
            
            
            %%%
            
            %% Eriksen Eriksen Post All Cells
            shuffBetas1=v5(randperm(length(v5)));
            %shuffBetas2=v6(randperm(length(v6)));
            % shuffBetas2=w(randperm(length(w)));
            [rvals, ~]=corr(shuffBetas1, v6, 'type', 'spearman');
            [rvals_pears, ~]=corr(shuffBetas1, v6, 'type', 'pearson');
            [rvals_kend, ~]=corr(shuffBetas1, v6, 'type', 'kend');
            %dotprod=dot(v6, I)/(norm(v6)*norm(I));
            %  eucdist=euc_dist(v6, I);
            
            % size(rvalue_precue_perm_null_eriksen_simon)
            [rvalue_null_simon_simon_post_spears_mdl_three_PredRT_NoNan(p)] = rvals;
            [rvalue_null_simon_simon__postpears_mdl_three_PredRT_NoNan(p)] = rvals_pears;
            [rvalue_null_simon_simon_post_kend_mdl_three_PredRT_NoNan(p)] = rvals_kend;
            % v6=v6(randperm(length(v6)));
            %v62=v6(randperm(length(v6)));
            
            cosang=max(min(dot((v6),(shuffBetas1))/(norm(v6)*norm(shuffBetas1)),1),-1);
            cosang_cent=max(min(dot((v6-mean(v6)),(shuffBetas1-mean(shuffBetas1)))/(norm(v6-mean(v6))*norm(shuffBetas1-mean(shuffBetas1))),1),-1);
            
            % cosang__null_eriksen_simon_mdl_three_PredRT_Dev(p)=dotprod;
            Ang_null_mdl_three_EriksenVEriksen_NoNan(p) = real(acosd(max(min(dot((v6),(shuffBetas1))/(norm(v6)*norm(shuffBetas1)),1),-1)));
            AngCent_null_mdl_three_EriksenVEriksen_NoNan(p) = real(acosd(max(min(dot((v6-mean(v6)),(shuffBetas1-mean(shuffBetas1)))/(norm(v6-mean(v6))*norm(shuffBetas1-mean(shuffBetas1))),1),-1)));
            CosAng_null_mdl_three_EriksenVEriksen_NoNan(p)=cosang;
            CosAngCent_null_mdl_three_EriksenVEriksen_NoNan(p)=cosang_cent;
            %%%
            
            
        end
   
        
        
    end
    
%% More permutation tests

    K=cell2mat(Betas_mdl_three_PredRT_Eriksen_NoNan(:,ii));
    L=cell2mat(Betas_mdl_three_PredRT_Simon_NoNan(:,ii));
   
 
    K=cell2mat(Betas_mdl_three_PredRT_Eriksen_NoNan(:,ii));
    L=cell2mat(Betas_mdl_three_PredRT_Simon_NoNan(:,ii));
    
    
    for p = 1:numPerms
        
        
        
        shuffBetas=K(randperm(length(K)));
        [rvals, ~]=corr(shuffBetas, L, 'type', 'spearman');
        [rvals_pears, ~]=corr(shuffBetas, L, 'type', 'pearson');
        [rvals_kend, ~]=corr(shuffBetas, L, 'type', 'kend');
        dotprod=dot(shuffBetas, L)/(norm(shuffBetas)*norm(L));
        %eucdist=euc_dist(shuffBetas, L);
        
        % size(rvalue_precue_perm_null_eriksen_simon)
        [rvalue_null_eriksen_simon_spears_mdl_three_PredRT_NoNan(p,ii)] = rvals;
        [rvalue_null_eriksen_simon_pears_mdl_three_PredRT_NoNan(p,ii)] = rvals_pears;
        [rvalue_null_eriksen_simon_kend_mdl_three_PredRT_NoNan(p,ii)] = rvals_kend;
        cosang_null_eriksen_simon_NoInt_mdl_three_PredRT_NoNan(p)=dotprod;
        
        


    end

    
    
    [pvalue_null_eriksen_simon_spears_mdl_three_PredRT_NoNan(ii)] = 1-numel(find(abs(rvalue_null_eriksen_simon_spears_mdl_three_PredRT_NoNan(:,ii))<abs(rvalue_Spear_mdl_three_SimonVEriksen_NoNan(ii))))/numPerms;
    [pvalue_null_eriksen_simon_pears_mdl_three_PredRT_NoNan(ii)] = 1-numel(find(abs(rvalue_null_eriksen_simon_pears_mdl_three_PredRT_NoNan(:,ii))<abs(rvalue_Pears_mdl_three_SimonVEriksen_NoNan(ii))))/numPerms;
    [pvalue_null_eriksen_simon_kend_mdl_three_PredRT_NoNan(ii)] = 1-numel(find(abs(rvalue_null_eriksen_simon_kend_mdl_three_PredRT_NoNan(:,ii))<abs(rvalue_Kend_mdl_three_SimonVEriksen_NoNan(ii))))/numPerms;
    %  [pvalue_null_eriksen_simon_cosang_mdl_three_PredRT_NoNan(ii)] = 1-numel(find(abs(cosang_null_eriksen_simon_NoInt_mdl_three_PredRT_NoNan(:,ii))<abs(Cos_ang_SimVErik.mdl_three_PredRT_NoNan(ii))))/numPerms;

keyboard
end



