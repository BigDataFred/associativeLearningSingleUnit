function [lfp,erp,delIx,selIx,n,chanLabLFP] = preprocLFP( lfpDat, trlPool, trlIdx )

%% extract the labels of Behnke-Fried electrodes
BFlab = cell(length(lfpDat.chanLab),1);
for curChan = 1:length( lfpDat.chanLab )% loop over microwires
    tmp = lfpDat.chanLab{curChan};
    ix = regexp(tmp,'\d{1}');
    BFlab(curChan) = {tmp(1:ix-1)};% only keep the BF label, ignore MW index
end;
BFid = unique(BFlab);% extract unique labels

BFix = cell(length(BFid),1);
for curBF = 1:length(BFid)
    BFix{curBF} = find(strcmp(BFid{curBF},BFlab));% map micro-wires to BF-label (ie extract indexes)
end;

%% calculate root mean square energy for micro-wires
for curBF = 1:length(BFid)
    rms = zeros(length(BFix{curBF}),1);
    for curMW = 1:length(BFix{curBF})
        rms(curMW) = mean(sqrt(mean(lfpDat.LFPseg{BFix{curBF}(curMW)}(lfpDat.dsTrlTime>=-0.5 & lfpDat.dsTrlTime <=5,:).^2,1)));
    end;
    rms = (rms-mean(rms))./std(rms);
    %delete Micro-wires with RMS above threshold
    lfpDat.LFPseg(BFix{curBF}(abs(rms)>3)) = [];
    lfpDat.chanLab(BFix{curBF}(abs(rms)>3)) = [];
end;

%%
[ lfp ]   = cell(1,length( lfpDat.LFPseg));% lfp data
[ erp ]   = cell(1,length( lfpDat.LFPseg));% erp data
[ delIx ] = cell(1,length( lfpDat.LFPseg));% indexes with delete lfp trials (important for spk data!)
[ n ]     = zeros(1,length( lfpDat.LFPseg));

%% pre-allocate and initialize
parfor curChan = 1:length( lfpDat.LFPseg )%loop over channels
        
    [ tmpLFP ] = lfpDat.LFPseg{curChan}(:,ismember(lfpDat.trlSel, trlIdx )); % extract trials corresponding to encoding
    
    %% create average signal of all channels expect current MW channel
    tmpAVGsig = 0;
    
    [ selIx ] = setdiff( 1:length( lfpDat.LFPseg ), curChan );% include indexes of all channels, except current MW
    for curChan2 = 1:length( selIx )
        tmpAVGsig = tmpAVGsig + lfpDat.LFPseg{selIx(curChan2)};
    end;
    tmpAVGsig = tmpAVGsig/length(selIx);
    
    [ tmpAVGsig] = tmpAVGsig( : ,ismember( lfpDat.trlSel, trlIdx ) );
    
    %% detect trials with artefacts
    % compute RMS of LFP over time
    rms = sqrt(mean(tmpLFP(lfpDat.dsTrlTime>=-0.5 & lfpDat.dsTrlTime <=5,:).^2,1));
    rms = ( rms-mean(rms) )./std(rms);% standardize RMS
    
    % compute z-score of LFP
    m  = ones(size(tmpLFP,1),1)*mean(tmpLFP,1);
    sd = ones(size(tmpLFP,1),1)*std(tmpLFP,0,1);
    z = ( tmpLFP-m )./ sd;
    z = max(abs(z(lfpDat.dsTrlTime>=-0.5&lfpDat.dsTrlTime <=5,:)),[],1);
    
    delIx{curChan} = unique([find(abs(rms) >= 4) find(z >= 4)]);% keep indexes of trials where LFP-activity is out of range (outliers)
    
    %% delete trials with artefacts
    tmpLFP(:,delIx{curChan}) = []; %delete trials where LFP-activity is out of range (outliers)
    tmpAVGsig(:,delIx{curChan}) = [];
    
    %% remove activity from LFP that is phase-locked to average activity, trial by trial
    [ tmpLFP ] = orthogonalizeTimeDomain( tmpAVGsig , tmpLFP );
    
    %%
    % LFP-data must have at least 25 trials in total (hits + misses) per session to perform spectral analysis
    tmptrlPool = trlPool;
    tmptrlPool(delIx{curChan}) = [];
    n(curChan) = length(find(ismember(tmptrlPool,trlIdx)));
    if (~isempty( tmpLFP) ) && (n(curChan)>25)
        
        m  = ones(size(tmpLFP,1),1)*mean(tmpLFP,1); % mean of LFP over time
        %sd = ones(size(tmpLFP,1),1)*std(tmpLFP,0,1); % SD of LFP over time
        
        lfp{curChan} = (tmpLFP-m);%./sd;% mean-centered single trial lfp data
        erp{curChan} = lfp{curChan}; % evoked poential
    end;
end;
clear LFPseg;

%%
selIx = [];
for curChan = 1:length( lfp )
    if ~isempty(lfp{curChan})
        selIx = [selIx;curChan];
    end;
end;
chanLabLFP = lfpDat.chanLab(selIx);

lfp = lfp(selIx);
erp = erp(selIx);
delIx = delIx(selIx);
n     = n(selIx);


