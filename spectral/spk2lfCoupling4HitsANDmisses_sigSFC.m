function [Npairs]=spk2lfCoupling4HitsANDmisses_sigSFC( pId, expMode, timewindows,spk2LFPmode, rdsPath, savePath, Notch)

%%
addpath(genpath('~/tbx/chronux_2_11/spectral_analysis/'));
addpath('/media/hanslmas/rds-share/fred4simon/tbx/fieldtrip-master');
ft_defaults;

%%
[ p2SPKd ] = savePath;
[ p2spk2LFPplvDat ] = savePath;

%%

if isempty( gcp('Nocreate') )
    parpool(36,'SpmdEnabled',false);
end;

%%
counter=0;
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
                        %% load the spike data
                        [spkDatfN ] = dir( [ p2SPKd, pId{curPat}, '_', expMode{curExp}, '_', sesh{curSesh}, '_spkParams_Sorting_scr.mat' ] );
                        if ~isempty(spkDatfN)
                            load( [ p2SPKd, spkDatfN.name ] )

                            [ lfpDat ] = load( [ p2d, lfpDatfN.name ] ); % load the LFP data
                            [ trlPool, hitIdx, missIdx, trlENC ] = organizeTrlIdxEM( lfpDat );

                            %% preproc LFP channels
                            [ lfp, ~, delIx, ~, ~, chanLabLFP ] = preprocLFP( lfpDat, trlPool, trlENC, Notch );% preprocess LFP

                            %% extract phase from LFP
                            [ phi, phiTime, phiFreq ] = computeLFPphase( lfp, lfpDat.dsFs, [min(lfpDat.dsTrlTime) max(lfpDat.dsTrlTime)],'LF');
                            clear lfp;
                            counter=counter+1;
                            %% Get time periods of interest    
                            for tw=1:size(timewindows,1)
                                tmptw=timewindows(tw,:);    
                                ptIx = find( phiTime >=tmptw(1) & phiTime < tmptw(2) );
                                stIx = find( spkDat.dt2(2:end-1) >=tmptw(1)*1000 & spkDat.dt2(2:end-1) < (tmptw(2)*1000));
                                fIx = find( phiFreq >=2 & phiFreq <=40);
                                spk2LFPfreqAx = phiFreq(fIx);
                            
                                %% compute spk2LFP coupling for Hits & Misses
                                for k=1:length(phi)
                                    phitmp{k}=phi{k}(:,:,:,ptIx);
                                end
                                nSPKchans=length(spkDat.chanLabSPK);
                                nLFPchans=length(chanLabLFP);
                                nCombs=nSPKchans*nLFPchans;
                                paircnt=0;spaircnt=0;
                                for curSpk = 1:nSPKchans
                                    fprintf( [ num2str(curSpk), '/', num2str(nSPKchans) ] );
                                    spklab=spkDat.chanLabSPK(curSpk);
                                    clusID=spkDat.clusID(curSpk);
                                    for curMW = 1:nLFPchans
                                        if ~isempty( phitmp{curMW} )
                                            lfplab=chanLabLFP(curMW);
                                            [ tmp ] = squeeze( phitmp{curMW}(:,:,fIx,:) );
                                            [ ts ]  = spkDat.spkTs{curSpk}(stIx,:);
                                            ts(:,delIx{curMW}) = [];
                                            %ts=ts(:);
                                            [ trlPoolTMP ] = trlPool;
                                            trlPoolTMP(delIx{curMW}) = [];
                                            hIx = find(ismember(trlPoolTMP,hitIdx));
                                            mIx = find(ismember(trlPoolTMP,missIdx));
                                            ix1 = find(ts(:,hIx)==1);
                                            ix2 = find(ts(:,mIx)==1);
                                            minspkchk=min([length(ix1) length(ix2)]);
                                            if minspkchk>=30
                                                paircnt=paircnt+1;
                                                tmp1 = permute(tmp(hIx,:,:),[3 1 2]);
                                                tmp2 = permute(tmp(mIx,:,:),[3 1 2]);
                                                tmp1 = reshape(tmp1,[size(tmp1,1)*size(tmp1,2) length(fIx)]);
                                                tmp2 = reshape(tmp2,[size(tmp2,1)*size(tmp2,2) length(fIx)]);
                                                tmpall=angle([tmp1(ix1,:);tmp2(ix2,:)]);
                                                for freq=1:size(tmpall,2)
                                                    [pval(freq), ~] = circ_rtest(tmpall(:,freq));
                                                end
                                                [~,pcor,~] = fdr(pval);
                                                sigchk=find(pcor<0.05);
                                                if ~isempty(sigchk)
                                                    spaircnt=spaircnt+1;
                                                    chan_labels(spaircnt,1)= spklab;% first entry is spike label
                                                    chan_labels(spaircnt,2)= {num2str(clusID)}; % second entry is cluster ID
                                                    chan_labels(spaircnt,3)=lfplab;% third entry is lfp label
                                                    [ spk2LFPCouplingHM(spaircnt,:) ] = computeSPK2LFPcoupling_S( [tmp1(ix1,:);tmp2(ix2,:)], spk2LFPmode{curSpk2LFPmode} );
                                                    [ spk2LFPCouplingH(spaircnt,:) ] = computeSPK2LFPcoupling_S( tmp1(ix1,:), spk2LFPmode{curSpk2LFPmode} );
                                                    [ spk2LFPCouplingM(spaircnt,:) ] = computeSPK2LFPcoupling_S( tmp2(ix2,:), spk2LFPmode{curSpk2LFPmode} );
                                                    sigfs{spaircnt,1}=sigchk;
                                                    nspksHM(spaircnt,:)=[length(ix1) length(ix2)];
                                                end
                                            end
                                        end
                                    end
                                end
                                fprintf('\n');
                                % save Number of possible pairs for
                                   % later Binom test
                                   
                                Npairs(counter,tw)=paircnt;
                               if spaircnt > 1
                                   %% save results to disk
                                   timebit=[num2str(tmptw(1)), 'to', num2str(tmptw(2)), 'Secs'];
                                   saveName =[ pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spk2LFPCouplingHitsANDmisses_',...
                                       spk2LFPmode{curSpk2LFPmode},'_Sorting_2to40Hz_FDRfs_',timebit,'_NotchOff.mat' ];
                                   save([savePath,saveName],'spk2LFPfreqAx','spk2LFPCouplingH','spk2LFPCouplingM','spk2LFPCouplingHM','chan_labels', 'nspksHM', 'sigfs', 'paircnt');
                                   clear spk2LFPCoupling* chan_labels nspksHM paircnt spaircnt;
                               end
                            end
                        clear delIx;
                        end
                    
                    end
                
                end
            
            end
        
        end
    
    end
end
