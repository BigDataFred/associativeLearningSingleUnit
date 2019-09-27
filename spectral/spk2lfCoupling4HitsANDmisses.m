function spk2lfCoupling4HitsANDmisses( pId, expMode, spkMode, spk2LFPmode, savePath, alpha, nRand )

%%
addpath(genpath('~/tbx/chronux_2_11/spectral_analysis/'));
addpath('~/tbx/fieldtrip-20170618/');
ft_defaults;

%%
[ p2SPKd ] = savePath;

%%
if isempty( gcp('Nocreate') )
    parpool(36,'SpmdEnabled',false);
end;

%%
for curPat = 1:length(pId)
    for curExp = 1:length(expMode)
        
        [ tmp ] = dir( [ rdsPath, pId{curPat},'/', expMode{curExp}, '/' ] );
        
        if ~isempty(tmp)
            
            [ sesh ] = extractSeshLabels( tmp );
            clear tmp;
            
            %%
            for curSesh = 1:length(sesh)
                
                p2d = [ rdsPath,pId{curPat}, '/' ,expMode{curExp}, '/',sesh{curSesh}, '/' ];
                fN = dir( [ p2d,pId{curPat}, '_',expMode{curExp}, '_',sesh{curSesh}, '_lfpDataStimLockedSegmenteddownsampled.mat' ] );
                
                [ lfpDat ] = load( [ p2d, fN.name ] )% load the LFP data
                
                %%
                for curSpk2LFPmode = 1:length(spk2LFPmode)
                    
                    for curSpkSortingMode = 1:length(spkMode)
                        
                        %% load the spk2LFP coupling data
                        p2spk2LFPplvDat = savePath;
                        fNspk2LFP =[ pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spk2LFPCoupling_',spk2LFPmode{curSpk2LFPmode},'_',spkMode{curSpkSortingMode},'_alpha',num2str(alpha),'_nRand',num2str(nRand),'_lowFreq.mat' ];
                        fNspk2LFP = dir([p2spk2LFPplvDat,fNspk2LFP]);
                        
                        [ spk2LFPplvDat ] = load([p2spk2LFPplvDat,fNspk2LFP.name]);
                        
                        %% Readout sig LFP channels
                        [ lfpDat.LFPseg ] = lfpDat.LFPseg(spk2LFPplvDat.sigIxLFP)
                        
                        %%
                        [ trlPool, hitIdx, missIdx, trlENC ] = organizeTrlIdxEM( lfpDat );
                        
                        [ lfp, ~, delIx, ~, ~, chanLabLFP ] = preprocLFP( lfpDat, trlPool, trlENC );% preprocess LFP
                        
                        %% extract phase from LFP
                        [ phi, phiTime, phiFreq ] = computeLFPphase( lfp, lfpDat.dsFs, [min(lfpDat.dsTrlTime) max(lfpDat.dsTrlTime)] );
                        
                        %% load the spike data
                        fN = dir( [ p2SPKd, pId{curPat}, '_', expMode{curExp}, '_', sesh{curSesh}, '_spkParams_', spkMode{curSpkSortingMode},'.mat' ] )
                        [ spkDat ] = load( [ p2SPKd, fN.name ] )
                        
                        %% Readout sig SPK channels
                        [ spkTs ] = spkDat.spkTs( spk2LFPplvDat.sigIxSPK ) ;
                        
                        %% safety checks
                        chck = zeros(1,length(phi));
                        for curChan = 1:length(phi)
                            chck(curChan) = ~isempty(phi{curChan});
                        end;
                        
                        chck2 = zeros(1,length(spkTs));
                        for curChan = 1:length(spkTs)
                            chck2(curChan) = ~isempty(spkTs{curChan});
                        end;
                        clear spkDat;
                        clear lfpDat;
                        
                        %%
                        if any(chck~=0) && any(chck2~=0)
                            
                            %%
                            for curMW = 1:length(phi)
                                if ~isempty(phi{curMW})
                                    tIx = find( phiTime >=-.5 & phiTime < 5 );
                                    fIx = find( phiFreq >=0 & phiFreq <=30 );
                                    spk2LFPfreqAx = phiFreq;
                                    break;
                                end;
                            end;
                                                        
                            %% compute spk2LFP coupling for Hits & Misses
                            [ spk2LFPCouplingH ] = NaN( length( sigIxLFP ), length( fIx ) );
                            [ spk2LFPCouplingM ] = NaN( length( sigIxLFP ), length( fIx ) );
                            if ~isempty( sigIxLFP ) && ~isempty( sigIxSPK )
                                parfor curMW = 1:length( sigIxLFP )
                                    fprintf( [ num2str(curMW), '/', num2str(length(sigIxLFP)) ] );
                                    if ~isempty( phi{curMW} )
                                        [ tmp ] = squeeze( phi{curMW}(:,:,fIx,tIx) );
                                        
                                        ts = spkTs{curMW}(tIx,:);
                                        ts(:,delIx{curMW}) = [];
                                        ix = find(ts==1);
                                        ts = [];
                                        
                                        [ spk2LFPCouplingH(curMW,:) ] = computeSPK2LFPcoupling( tmp(hitIdx,:,:), ix, spk2LFPmode{curSpk2LFPmode} );
                                        [ spk2LFPCouplingM(curMW,:) ] = computeSPK2LFPcoupling( tmp(missIdx,:,:), ix, spk2LFPmode{curSpk2LFPmode} );
                                    end;
                                    fprintf('\n');
                                end;
                            end;
                            
                            %%
                            chanLabLFPsig = [];
                            chanLabSPKsig = [];
                            if ~isempty (sigIxLFP) && ~isempty(sigIxSPK)
                                [ chanLabLFPsig ] = chanLabLFP;
                                [ chanLabSPKsig ] = chanLabSPK;
                            end;
                            
                            %%
                            saveName =[ pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spk2LFPCouplingHitsANDmisses_',spk2LFPmode{curSpk2LFPmode},'_',spkMode{curSpkSortingMode},'_alpha',num2str(alpha),'_nRand',num2str(nRand),'_lowFreq.mat' ];
                            save([savePath,saveName],'spk2LFPfreqAx','spk2LFPCouplingH','spk2LFPCouplingM','chanLabLFPsig','chanLabSPKsig','sigIxLFP','sigIxSPK','-v7.3');
                            
                            clear spk2LFPfreqAx spk2LFPCoupling* chanLabLFP chanLabSPK sigIxLFP sigIxSPK;
                            
                        end;
                    end;
                end;
            end;
        end;
    end;
end;



