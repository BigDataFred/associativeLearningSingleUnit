function computeSPKparamsEMtask( pId, expMode, spkMode, savePath, rdsPath  )

%% start parallel ppol
if isempty(gcp('nocreate'))
    parpool(36,'SpmdEnabled', false);
end;

%%
for curMode = 1:length(spkMode)
    for curPat = 1:length(pId) % loop over patients
        for curExp = 1:length(expMode) % loop over exp-modes (either fVSpEM or cnEM)
            
            [ tmp ] = dir([rdsPath,pId{curPat},'/',expMode{curExp},'/']);% check for existing data
            
            if ~isempty(tmp) % get the timestamp labels for each session
                
                [sesh] = extractSeshLabels(tmp);clear tmp;
                
                %%
                pId{curPat}
                expMode{curExp}
                sesh'
                
                %%
                for curSesh = 1:length(sesh) % loop over the recording-sesssions
                    
                    %% get the spk data corresponding to the current session
                    [ p2d ] = [rdsPath,pId{curPat},'/',expMode{curExp},'/',sesh{curSesh},'/']; %path for preproc data
                    
                    [ fN ]     = dir([p2d,pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_lfpDataStimLockedSegmenteddownsampled.mat']);% filename of preproc data
                    [ lfpDat ] = load([p2d,fN.name])% load data
                    [ dt ]     = min(lfpDat.dsTrlTime).*1e3-1:1:max(lfpDat.dsTrlTime).*1e3+1;    
                    [ trlENC ] = lfpDat.trlENC;
                    
                    [ fN ]     = dir([p2d,pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spkDataStimLockedSegmented.mat']) %check for existing data
                    [ spkDat ] = load([p2d,fN.name])% load data
                    
                    %%
                    [ SPKseg, spkTs, frM, frSD, xc, chanLabSPK ] = computeSPKParams( spkDat, spkMode{curMode}, trlENC, -300:300, dt );
                    
                    %% find units that are active during encoding
                    [ spkSelIx ] = filterUnits4SPKrate( spkTs, frM, [0. 4.], 100, lfpDat.dsTrlTime );
                    
                    %% select units that are active during encoding
                    [SPKseg]       = SPKseg(spkSelIx,:);
                    [xc]           = xc(spkSelIx,:);
                    [ spkTs ]      = spkTs(spkSelIx);
                    [ frM ]        = frM(spkSelIx,:);
                    [ frSD ]       = frSD(spkSelIx,:);
                    [ chanLabSPK ] = chanLabSPK(spkSelIx);
                    
                    %% save spk-params to file
                    saveName = [pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spkParams_',spkMode{curMode},'.mat'];
                    save([savePath,saveName],'SPKseg','spkTs','frM','frSD','xc','chanLabSPK','dt','dt2');
                    
                end;
            end;
        end;
    end;
end;