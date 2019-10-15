function spk2lfCoupling4HitsANDmisses( pId, expMode, spkMode, spk2LFPmode, rdsPath, savePath, alpha, nRand, stratMode )

%%
addpath(genpath('~/tbx/chronux_2_11/spectral_analysis/'));
addpath('~/tbx/fieldtrip-20170618/');
ft_defaults;

%%
[ p2SPKd ] = savePath;
[ p2spk2LFPplvDat ] = savePath;

%%

if isempty( gcp('Nocreate') )
    parpool(36,'SpmdEnabled',false);
end;

%%
for curPat = 1:length( pId )
    for curExp = 1:length( expMode )
        
        [ tmp ] = dir( [ rdsPath, pId{curPat},'/', expMode{curExp}, '/' ] );
        
        if ~isempty(tmp)
            
            [ sesh ] = extractSeshLabels( tmp );
            clear tmp;
            
            %%
            for curSesh = 1:length( sesh )
                
                [ p2d ] = [ rdsPath, pId{curPat}, '/' ,expMode{curExp}, '/',sesh{curSesh}, '/' ];
                [ lfpDatfN ] = dir( [ p2d, pId{curPat}, '_',expMode{curExp}, '_',sesh{curSesh}, '_lfpDataStimLockedSegmenteddownsampled.mat' ] );
                
                if ~isempty( lfpDatfN )
                    
                    %%
                    for curSpk2LFPmode = 1:length( spk2LFPmode )
                        
                        %%
                        for curSpkSortingMode = 1:length( spkMode )
                            
                            %% load the spk2LFP coupling data
                            [ fNspk2LFP ] =[ pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spk2LFPCoupling_plv_',spkMode{curSpkSortingMode},'_alpha',num2str(alpha),'_nRand',num2str(nRand),'_lowFreq.mat' ];
                            [ fNspk2LFP ] = dir( [ p2spk2LFPplvDat, fNspk2LFP ] );
                            
                            if ~isempty( fNspk2LFP )
                                
                                [ spk2LFPplvDat ] = load([ p2spk2LFPplvDat,fNspk2LFP.name ]);
                                [ chanLabLFPsig ] = spk2LFPplvDat.chanLabLFPsig;
                                
                                %%
                                if ~isempty( chanLabLFPsig )
                                    
                                    %%
                                    [ lfpDat ] = load( [ p2d, lfpDatfN.name ] ) % load the LFP data
                                    [ trlPool, hitIdx, missIdx, trlENC ] = organizeTrlIdxEM( lfpDat );
                                    
                                    %% preproc LFP channels
                                    [ lfp, ~, delIx, ~, ~, chanLabLFP ] = preprocLFP( lfpDat, trlPool, trlENC );% preprocess LFP
                                    
                                    %%
                                    [ sigIxLFP ] = zeros( length( spk2LFPplvDat.chanLabLFPsig ), 1 );
                                    for curMW = 1:length( spk2LFPplvDat.chanLabLFPsig )
                                        [ sigIxLFP(curMW) ] = find(strcmp( chanLabLFP, chanLabLFP{ spk2LFPplvDat.sigIxLFP(curMW) } ) );
                                    end;
                                    
                                    if ~isequal( sigIxLFP, spk2LFPplvDat.sigIxLFP )
                                        error('Index mismatch detected')
                                    end;
                                    
                                    if ~( strcmp(chanLabLFP(sigIxLFP),chanLabLFPsig) )
                                        error('Index mismatch detected')
                                    end;
                                    
                                    [ delIxTmp ] = delIx( sigIxLFP );
                                    [ lfpTmp ]   = lfp( sigIxLFP );
                                    
                                    %% extract phase from LFP
                                    [ phi, phiTime, phiFreq ] = computeLFPphase( lfpTmp, lfpDat.dsFs, [min(lfpDat.dsTrlTime) max(lfpDat.dsTrlTime)] );
                                    clear lfpTmp;
                                    
                                    %% load the spike data
                                    [spkDatfN ] = dir( [ p2SPKd, pId{curPat}, '_', expMode{curExp}, '_', sesh{curSesh}, '_spkParams_', spkMode{curSpkSortingMode},'.mat' ] )
                                    [ spkDat ] = load( [ p2SPKd, spkDatfN.name ] )
                                    
                                    %% Readout sig SPK channels
                                    [ sigIxSPK ]      = spk2LFPplvDat.sigIxSPK;
                                    [ spkTs ]         = spkDat.spkTs(sigIxSPK) ;
                                    [ chanLabSPKsig ] = spkDat.chanLabSPK(sigIxSPK);
                                    
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
                                        
                                        %%
                                        for curStratMode = 1:length( stratMode )
                                            
                                            %% compute spk2LFP coupling for Hits & Misses
                                            [ spk2LFPCouplingH ] = NaN( length( sigIxLFP ), length( fIx ) );
                                            [ spk2LFPCouplingM ] = NaN( length( sigIxLFP ), length( fIx ) );
                                            if ~isempty( sigIxLFP ) && ~isempty( sigIxSPK )
                                                for curMW = 1:length( sigIxLFP )
                                                    fprintf( [ num2str(curMW), '/', num2str(length(sigIxLFP)) ] );
                                                    if ~isempty( phi{curMW} )
                                                        
                                                        [ tmp ] = squeeze( phi{curMW}(:,:,fIx,tIx) );
                                                        
                                                        [ ts ]  = spkTs{curMW}(tIx,:);
                                                        ts(:,delIxTmp{curMW}) = [];
                                                        
                                                        [ trlPoolTMP ] = trlPool;
                                                        trlPoolTMP(delIxTmp{curMW}) = [];
                                                        hIx = find(ismember(trlPoolTMP,hitIdx));
                                                        mIx = find(ismember(trlPoolTMP,missIdx));
                                                        
                                                        if (  strcmp( stratMode{curStratMode},'On') || strcmp( stratMode{curStratMode},'on') )% stratifsy number of trials
                                                            
                                                            if ( length( hIx ) > length( mIx ) )
                                                                rIx = randperm( length( hIx ) );% random permute trial indexes
                                                                cnt = 0;
                                                                nsamp = fix( length( hIx )/ length( mIx ) );
                                                                hIx2 = [];
                                                                for curSamp = 1:nsamp
                                                                    if ~ismepty( rIx ) && ( length( rIx) >= length( mIx ) )
                                                                        cnt = cnt+1;
                                                                        hIx2(cnt,:) = hIx( rIx( 1:length( mIx ) ) );% equate number of hits and misses
                                                                        rIx(1:length( mIx )) = [];
                                                                    end;
                                                                end;
                                                                hIx = hIx2; clear hIx2;
                                                            elseif ( length( hIx ) < length( mIx ) )
                                                                rIx = randperm( length( mIx ) );% random permute trial indexes
                                                                cnt = 0;
                                                                nsamp = fix( length( mIx )/ length( hIx ) );
                                                                mIx2 = [];
                                                                for curSamp = 1:nsamp
                                                                    if ~ismepty( rIx ) && ( length( rIx) >= length( hIx ) )
                                                                        cnt = cnt+1;
                                                                        mIx2(cnt,:) = mIx( rIx( 1:length( hIx ) ) );% equate number of hits and misses
                                                                        rIx(1:length( hIx )) = [];
                                                                    end;
                                                                end;
                                                                mIx = mIx2;clear mIx2;
                                                            end;
                                                            
                                                            x = [];
                                                            for curSamp = 1:size( hIx,1) 
                                                                ix1 = find( ts(:,hIx(curSamp,:) ) ==1 );
                                                                [ x(curSamp,:) ] = computeSPK2LFPcoupling( tmp(hIx(curSamp,:),:,:), ix1, spk2LFPmode{curSpk2LFPmode} );
                                                            end;
                                                            [ spk2LFPCouplingH(curMW,:) ] = mean(x,1);
                                                            
                                                            x = [];
                                                            for curSamp = 1:size( mIx,1) 
                                                                ix2 = find(ts(:,mIx(curSamp,:))==1);
                                                                [ x(curSamp,:) ] = computeSPK2LFPcoupling( tmp(mIx(curSamp,:),:,:), ix2, spk2LFPmode{curSpk2LFPmode} );
                                                            end;
                                                            [ spk2LFPCouplingM(curMW,:) ] = mean(x,1);
                                                            
                                                            
                                                            ts = [];
                                                                                                                                                                                   
                                                        else
                                                            
                                                            ix1 = find(ts(:,hIx)==1);
                                                            ix2 = find(ts(:,mIx)==1);
                                                            ts = [];
                                                            
                                                            [ spk2LFPCouplingH(curMW,:) ] = computeSPK2LFPcoupling( tmp(hIx,:,:), ix1, spk2LFPmode{curSpk2LFPmode} );
                                                            [ spk2LFPCouplingM(curMW,:) ] = computeSPK2LFPcoupling( tmp(mIx,:,:), ix2, spk2LFPmode{curSpk2LFPmode} );
                                                        end;
                                                    end;
                                                    fprintf('\n');
                                                end;
                                            end;
                                            
                                            %% save results to disk
                                            saveName =[ pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spk2LFPCouplingHitsANDmisses_',spk2LFPmode{curSpk2LFPmode},'_',spkMode{curSpkSortingMode},'_alpha',num2str(alpha),'_nRand',num2str(nRand),'_stratMode',stratMode{curStratMode},'_lowFreq.mat' ];
                                            save([savePath,saveName],'spk2LFPfreqAx','spk2LFPCouplingH','spk2LFPCouplingM','chanLabLFPsig','chanLabSPKsig','sigIxLFP','sigIxSPK','-v7.3');
                                        end;
                                        clear delIxTmp;
                                        clear spk2LFPfreqAx spk2LFPCoupling* chanLabLFPsig chanLabSPKsig sigLFPix sigIxSPK;
                                    end;
                                end;
                            end;
                        end;
                    end;
                end;
            end
        end;
    end;
end;
