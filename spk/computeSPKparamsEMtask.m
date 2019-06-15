%%
addpath('~/prj/Bham/code/mcode/utils/');
addpath('~/tbx/CircStat2012a/');
addpath(genpath('/home/rouxf/tbx/eeglab14_1_1b/functions/sigprocfunc/'));
addpath(genpath('~/tbx/chronux_2_11/spectral_analysis/'));
addpath('~/tbx/fieldtrip-20170618/');
ft_defaults;

%%
pId = {'P02','P04','P05','P07','P08','P09','P03ERL','P22AMS','P23AMS'};%
expMode = {'fVSpEM','cnEM'};
savePath = '~/resultsSpikeFieldOct18II/';

%% start parallel ppol
if isempty(gcp('nocreate'))
    parpool(36,'SpmdEnabled', false);
end;

%%
[ spkMode ] = {'Sorting','noSorting'};% can be either noSorting or Sorting

%%
[ dt2 ] = -300:300;

%%
for curMode = 1:length(spkMode)
    for curPat = 1:length(pId) % loop over patients
        for curExp = 1:length(expMode) % loop over exp-modes (either fVSpEM or cnEM)
            
            [ tmp ] = dir(['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pId{curPat},'/',expMode{curExp},'/']);% check for existing data
            
            if ~isempty(tmp) % get the timestamp labels for each session
                
                cnt = 0;
                sesh = cell(1,length(tmp));
                for curSesh = 1:length(sesh)
                    if ~isempty(regexp(tmp(curSesh).name,'\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}'))
                        cnt = cnt+1;
                        ix = regexp(tmp(curSesh).name,'\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}');
                        ix2 = ix+18;
                        sesh(cnt) = {tmp(curSesh).name(ix:ix2)};
                    end;
                end;
                sesh(cnt+1:end) = [];
                
                %%
                pId{curPat}
                expMode{curExp}
                sesh'
                
                %%
                for curSesh = 1:length(sesh) % loop over the recording-sesssions
                    
                    %% get the spk data corresponding to the current session
                    [ p2d ] = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pId{curPat},'/',expMode{curExp},'/',sesh{curSesh},'/']; %path for preproc data
                    
                    [ fN ]     = dir([p2d,pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_lfpDataStimLockedSegmenteddownsampled.mat']);% filename of preproc data
                    [ lfpDat ] = load([p2d,fN.name])% load data
                    [ dt ]     = min(lfpDat.dsTrlTime).*1e3-1:1:max(lfpDat.dsTrlTime).*1e3+1;    
                    [ trlENC ] = lfpDat.trlENC;
                    
                    [ fN ]     = dir([p2d,pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spkDataStimLockedSegmented.mat']) %check for existing data
                    [ spkDat ] = load([p2d,fN.name])% load data
                    
                    %% pre-allocate and initialize
                    [ SPKseg ] = {}; [ spkTs ]  = {};                                                    
                    [ frM ]    = []; [ frSD ]   = [];
                    [ xc ]     = [];
                    
                    %%
                    uCnt = 0;
                    [ chanLabSPK ] = [];
                    for  curMW = 1:length(spkDat.sortedSpikesSEG ) % loop over micro-wires
                        
                        fprintf([num2str(curMW),'/',num2str(length(spkDat.sortedSpikesSEG ))]);
                        if strcmp(spkMode,'noSorting')
                            [ cID ] = 1;
                        else
                            [ cID ] = unique(spkDat.sortedSpikesSEG{curMW}.assignedClusterSeg);
                        end;
                        cID = cID(cID ~=0);
                        
                        for curUnit = 1:length(cID)
                            
                            uIx = [];
                            if strcmp(spkMode{curMode},'noSorting')
                                uIx = 1:length(spkDat.sortedSpikesSEG{curMW}.assignedClusterSeg);%takes all spike-times
                            elseif strcmp(spkMode{curMode},'Sorting')
                                uIx = find(spkDat.sortedSpikesSEG{curMW}.assignedClusterSeg == cID(curUnit));% only spike-times of cluster
                            else
                                error('spkMode must match either ''noSorting'' or ''Sorting'' ');
                            end;
                            
                            % extract time stamps for encoding trials
                            if ~isempty(uIx)
                                ts = spkDat.sortedSpikesSEG{curMW}.SpikeTimesSeg(uIx).*1e3;% convert from s to ms
                                trl = spkDat.sortedSpikesSEG{curMW}.trl(uIx);% trial indexes
                                uCnt = uCnt+1;% increment unit-counter
                                chanLabSPK = [chanLabSPK;spkDat.chanLab(curMW)];%unit-channel
                                
                                % init some vars
                                dum1 = cell(1,length( trlENC) );% init
                                dum2 = zeros(length(lfpDat.dsTrlTime),length( trlENC ));% init
                                fr = zeros(length( trlENC ),length(dt));% init
                                xcTrl= zeros(length( trlENC ),length(dt2));%init
                                parfor curTrl = 1:length( trlENC )% loop over encoding trials
                                    ts2 = ts;
                                    tsSel = ts2( trl == trlENC( curTrl ));% filter spike times
                                    dum1{curTrl} = tsSel;
                                    [n,~] = hist(tsSel,dt);% transform spike timestamps into binary data (0=no spike,1=spike)
                                    
                                    dum2(:,curTrl) = n(2:end-1);% get rid of bin-edge artifacts
                                    fr(curTrl,:) = conv(n,gausswin(251),'same')./0.251;% compute smoothed firing rate using gaussian kernel
                                    xcTrl(curTrl,:) = xcorr(n(2:end-1),300); % compute auto-correlation of spike-train
                                end;
                                SPKseg{uCnt}.trial = dum1;
                                 spkTs{uCnt}       = dum2;
                                frM(uCnt,:)        = mean(fr,1); % compute mean firing rage
                                frSD(uCnt,:)       = std(fr,0,1)./sqrt(size(fr,1)-1);% compute SD of firing rate
                                xc(uCnt,:)         = mean(xcTrl,1);% compute mean auto-corrlation over trials
                                
                            end;
                            
                        end;
                        
                        fprintf('\n');
                    end;
                    xc(:,dt2 ==0) = NaN;
                    
                    %% find units that are active during encoding
                    [ spkSelIx ] = [];
                    [ dsTrlTime ] = lfpDat.dsTrlTime;
                    parfor curMW = 1:length( spkTs )
                        fprintf([num2str(curMW),'/',num2str(length( spkTs ))]);
                        frM2 = frM;
                        x = spkTs{curMW}(dsTrlTime >=0 & dsTrlTime<=4,:);
                        if ( sum(x(:)) > 100 ) && ( mean(frM2(curMW,dsTrlTime >=0. & dsTrlTime<=4)) > 1 )
                            spkSelIx = [spkSelIx curMW]; 
                        end;
                        fprintf('\n');
                    end;
                    
                    %% select units that are active during encoding
                    [xc]           = xc(spkSelIx,:);
                    [ spkTs ]      = spkTs(spkSelIx);
                    [ frM ]        = frM(spkSelIx,:);
                    [ frSD ]       = frSD(spkSelIx,:);
                    [ chanLabSPK ] = chanLabSPK(spkSelIx);
                    
                    %% save spk-params to file
                    saveName = [pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spkParams_',spkMode{curMode},'.mat'];
                    save([savePath,saveName],'spkTs','frM','frSD','xc','chanLabSPK','dt','dt2');
                    
                end;
            end;
        end;
    end;
end;