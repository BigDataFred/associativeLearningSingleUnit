function [SPKseg,spkTs,frM,frSD,xc,chanLabSPK] = computeSPKParams(spkDat,spkMode,trlIx,dt)

%% pre-allocate and initialize
[ SPKseg ] = {}; [ spkTs ]  = {};
[ frM ]    = []; [ frSD ]   = [];
[ xc ]     = []; [ chanLabSPK ] = [];

%%
uCnt = 0;
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
            dum1 = cell(1,length( trlIx) );% init
            dum2 = zeros(length(lfpDat.dsTrlTime),length( trlIx ));% init
            fr = zeros(length( trlIx ),length(dt));% init
            xcTrl= zeros(length( trlIx ),length(dt));%init
            parfor curTrl = 1:length( trlIx )% loop over encoding trials
                ts2 = ts;
                tsSel = ts2( trl == trlIx( curTrl ));% filter spike times
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
xc(:,dt ==0) = NaN;