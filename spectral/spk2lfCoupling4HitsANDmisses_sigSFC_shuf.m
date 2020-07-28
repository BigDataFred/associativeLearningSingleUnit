function spk2lfCoupling4HitsANDmisses_sigSFC_shuf( pId, expMode, timewindows,spk2LFPmode, rdsPath, savePath, Notch,Nrand,locdist)

%%
addpath(genpath('~/tbx/chronux_2_11/spectral_analysis/'));
addpath('/media/hanslmas/rds-share/fred4simon/tbx/fieldtrip-master');
ft_defaults;

%%
[ p2SPKd ] = savePath;
[ p2spk2LFPplvDat ] = savePath;

%%

% if isempty( gcp('Nocreate') )
%     parpool(36,'SpmdEnabled',false);
% end;

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
                            [ phi, phiTime, phiFreq ] = computeLFPphase( lfp, lfpDat.dsFs, [min(lfpDat.dsTrlTime) max(lfpDat.dsTrlTime)],'HF');
                            clear lfp;
                            counter=counter+1;
                            %% Get time periods of interest    
                            for tw=1:size(timewindows,1)
                                tmptw=timewindows(tw,:);    
                                ptIx = find( phiTime >=tmptw(1) & phiTime < tmptw(2) );
                                stIx = find( spkDat.dt2(2:end-1) >=tmptw(1)*1000 & spkDat.dt2(2:end-1) < (tmptw(2)*1000));
                                fIx = [find( phiFreq >=40) & find(phiFreq <=80)];
                                spk2LFPfreqAx = phiFreq(fIx);
                            
                                %% compute spk2LFP coupling for Hits & Misses
                                for k=1:length(phi)
                                    phitmp{k}=phi{k}(:,:,:,ptIx);
                                end
                                nSPKchans=length(spkDat.chanLabSPK);
                                nLFPchans=length(chanLabLFP);
                                nCombs=nSPKchans*nLFPchans;
                                paircnt=0;
                                for curSpk = 1:nSPKchans
                                    spklab=spkDat.chanLabSPK(curSpk);
                                    clusID=spkDat.clusID(curSpk);
                                    for curMW = 1:nLFPchans
                                        if ~isempty( phitmp{curMW} )
                                            % Exclude distal or local
                                            % couplings based on input
                                            % argument locdist;
                                            % locdist = 1; --> only local couplings will be considered
                                            % locdist = 0; --> only distal couplings will be considered 
                                            %
                                            % (i.e.where spike and lfp are on same Behnke-Fried)
                                            lfplab=chanLabLFP(curMW);
                                            chk=strcmp(spklab{1,1}(1:end-1),lfplab{1,1}(1:end-1));
                                            if chk == locdist
                                                rncnt=0;
                                                nrnd=0;
                                                spk2LFPCouplingHM=cell(1,Nrand);
                                                [ tmp ] = squeeze( phitmp{curMW}(:,:,fIx,:) );
                                                [ ts_orig ]  = spkDat.spkTs{curSpk}(stIx,:);
                                                for nrnd=1:Nrand
                                                    shufvec = randperm(numel(trlPool));
                                                    ts_shuf = ts_orig(:,shufvec);
                                                    ts_shuf(:,delIx{curMW}) = [];
                                                    [ trlPoolTMP ] = trlPool;
                                                    trlPoolTMP(delIx{curMW}) = [];
                                                    hIx = find(ismember(trlPoolTMP,hitIdx));
                                                    mIx = find(ismember(trlPoolTMP,missIdx));
                                                    ix1 = find(ts_shuf(:,hIx)==1);
                                                    ix2 = find(ts_shuf(:,mIx)==1);
                                                    minspkchk(nrnd,1)=min([length(ix1) length(ix2)]);
                                                    if nrnd==1
                                                       paircnt=paircnt+1;
                                                       fprintf( [ num2str(paircnt), '/', num2str(nCombs) ] );
                                                    end
                                                    tmp1 = permute(tmp(hIx,:,:),[3 1 2]);
                                                    tmp2 = permute(tmp(mIx,:,:),[3 1 2]);
                                                    tmp1 = reshape(tmp1,[size(tmp1,1)*size(tmp1,2) length(fIx)]);
                                                    tmp2 = reshape(tmp2,[size(tmp2,1)*size(tmp2,2) length(fIx)]);
                                                    tmpall=angle([tmp1(ix1,:);tmp2(ix2,:)]);
                                                    for freq=1:size(tmpall,2)
                                                        [pval{1,nrnd}(freq), ~] = circ_rtest(tmpall(:,freq));
                                                    end
                                                    [~,pcor,~] = fdr(pval{1,nrnd});
                                                    sigchk=find(pcor<0.05);
                                                    if ~isempty(sigchk)
                                                        rncnt=rncnt+1;
                                                        %chan_labels(nrnd,1)= spklab;% first entry is spike label
                                                        %chan_labels(rncnt,2)= {num2str(clusID)}; % second entry is cluster ID
                                                        %chan_labels(rncnt,3)=lfplab;% third entry is lfp label
                                                        [ spk2LFPCouplingHM{1,nrnd} ] = computeSPK2LFPcoupling_S( [tmp1(ix1,:);tmp2(ix2,:)], spk2LFPmode{curSpk2LFPmode} );
                                                        %                                                             [ spk2LFPCouplingH{1,nrnd} ] = computeSPK2LFPcoupling_S( tmp1(ix1,:), spk2LFPmode{curSpk2LFPmode} );
                                                        %                                                             [ spk2LFPCouplingM{1,nrnd} ] = computeSPK2LFPcoupling_S( tmp2(ix2,:), spk2LFPmode{curSpk2LFPmode} );
                                                        %                                                             sigfs{1,nrnd}=sigchk;
                                                    end
                                                end
                                            spk2LFPshuf(paircnt,1).nsig=rncnt;
                                            spk2LFPshuf(paircnt,1).labels={spklab{1,1}, num2str(clusID), lfplab{1,1}};
                                            spk2LFPshuf(paircnt,1).HM=spk2LFPCouplingHM;
                                            spk2LFPshuf(paircnt,1).minspk=minspkchk;
                                            %spk2LFPshuf(paircnt,1).H=spk2LFPCouplingH;
                                            %spk2LFPshuf(paircnt,1).M=spk2LFPCouplingM;
                                            %spk2LFPshuf(paircnt,1).sigfs=sigfs;
                                            fprintf('\n');
                                            end
                                        end
                                        
                                    end
                                end
                                
                                if ~isempty(spkDat.chanLabSPK)
                                    %% save results to disk
                                    timebit=[num2str(tmptw(1)), 'to', num2str(tmptw(2)), 'Secs'];
                                    saveName =[ pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spk2LFPCouplingHitsANDmisses_',...
                                        spk2LFPmode{curSpk2LFPmode},'_Sorting_40to80Hz_FDRfs_',timebit,'_shuff_dist.mat' ];
                                    save([savePath,saveName],'spk2LFPfreqAx','spk2LFPshuf');
                                    clear spk2LFPshuf paircnt spaircnt;
                                end
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
