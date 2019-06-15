%%
addpath('~/prj/Bham/code/mcode/utils/');
addpath('~/tbx/CircStat2012a/');
addpath(genpath('~/tbx/chronux_2_11/spectral_analysis/'));
addpath(genpath('~/associativeLearningSingleUnit/'));
addpath('~/tbx/fieldtrip-20170618/');
ft_defaults;

%%
pId = {'P02','P04','P05','P07','P08','P09','P03ERL','P22AMS','P23AMS'};%{'P05'};%,'P09'{'P09'};%
expMode = {'fVSpEM','cnEM'};%
savePath = '~/resultsSpikeFieldOct18II/';

%%
if isempty( gcp('Nocreate') ) 
    parpool(36,'SpmdEnabled',false);
end;

%%
for curPat = 1:length(pId)
    for curExp = 1:length(expMode)
        
        %tmp = dir(['~/MICRO/pool/',pId{curPat},'_',expMode{curExp},'_*_lfpDataStimLockedSegmenteddownsampled.mat']);
        [ tmp ] = dir(['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pId{curPat},'/',expMode{curExp},'/']);
        
        if ~isempty(tmp)
            
            [sesh] = extractSeshLabels(tmp);
            clear tmp;
            
            %%
            for curSesh = 1:length(sesh)
                
                %p2d = ['~/MICRO/pool/'];
                p2d = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pId{curPat},'/',expMode{curExp},'/',sesh{curSesh},'/'];
                fN = dir([p2d,pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_lfpDataStimLockedSegmenteddownsampled.mat']);
                
                [ lfpDat ] = load([p2d,fN.name])                
                
                %%                
                [trlPool,hitIdx,missIdx,trlENC] = organizeTrlIdxEM(lfpDat);
                
                [lfp,erp,delIx,selIx,n,chanLabLFP] = preprocLFP( lfpDat, trlPool, hitIdx );
                
                %%
                [phi,phiTime,phiFreq] = computeLFPphase( lfp, lfpDat.dsFs, [min(lfpDat.dsTrlTime) max(lfpDat.dsTrlTime)] );
                
                %%
                fN = dir([p2d,pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spkParams_',spkMode{curMode},'.mat'])
                [ spkDat ] = load([p2d,fN.name])
              
                %%
                [ spkTs ] = spkDat.spkTs;
                
                [ spkSelIx ] = [];
                for curMW = 1:length( spkTs )
                    
                    x = spkTs{curMW}(lfpDat.dsTrlTime >=-0.5 & lfpDat.dsTrlTime<=5,:);
                    mSpkCnt = median(sum(spkTs{curMW}(lfpDat.dsTrlTime >0 & lfpDat.dsTrlTime <=4,:),1));
                    
                    if ( sum(x(:)) >= 50 ) && (mSpkCnt(curMW)>2) && ( mean(spkDat.frM(curMW,lfpDat.dsTrlTime >=-0.5 & lfpDat.dsTrlTime<=5)) > 1 )
                        spkSelIx = [spkSelIx curMW];
                    end;
                end; 
                
                [ spkTs ]      = spkTs(spkSelIx);
                [ chIx ]       = chIx(spkSelIx);                   
                [ chanLabSPK ] = lfpDat.chanLab(chIx);clear chIx;
                
                %%
                [ spkSelIx2 ] = [];
                for curMW = 1:length( spkTs )
                    x = spkTs{curMW}(lfpDat.dsTrlTime >=-0.5 & lfpDat.dsTrlTime<=5,:);
                    x(:,delIx{curMW}) = [];
                    if ( sum(x(:)) >= 50 )
                        spkSelIx2 = [spkSelIx2 curMW];
                    end;
                end;  clear x;              
                [ spkTs ]      = spkTs(spkSelIx2);               
                [ chanLabSPK ] = chanLabSPK(spkSelIx2);                                
                
                %%
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
                nRand = 200;
                if any(chck~=0) && any(chck2~=0)
                    
                    %%
                    spk2LFPmode = 'plv';
                    
                    [ spk2lfp ] = zeros(length(phi),length(spkTs));
                    [ spk2lfpRnd ] = zeros(length(phi),length(spkTs),nRand);
                    
                    for curMW = 1:length(phi)
                        if ~isempty(phi{curMW})
                            tIx = find( phi{curMW}.time >=-.5 & phi{curMW}.time < 5);
                            fIx = find( phi{curMW}.freq >=4 & phi{curMW}.freq <=18 );  
                            spk2LFPfreqAx = phi{curMW}.freq;
                            break;
                        end;
                    end;
                    
                   for curMW = 1:length(phi)
                        fprintf([num2str(curMW),'/',num2str(length(phi))]);
                                                                         
                        if ~isempty(phi{curMW})                                                        
                            [ tmp ] = squeeze(phi{curMW}.fourierspctrm(:,:,fIx,tIx));
                            [ phiTrl ] = size(tmp,1);
                            
                            [ mxRnd ] = zeros(length(spkTs),nRand);
                            [ mX ] = zeros(1,length(spkTs));       
                            
                            for curMW2 = 1:length( spkTs )
                                ts = spkTs{curMW2}(tIx,:);
                                ts(:,delIx{curMW}) = [];
                                ix = find(ts==1);
                                ts = [];
                                
                                switch spk2LFPmode
                                    case 'plv'
                                        [ plv ] = NaN( 1,length( fIx) );
                                        parfor curFreq = 1:length( fIx)
                                            ph = squeeze(tmp(:,curFreq,:))';
                                            ph = ph(:);
                                            ph = ph(ix);
                                            ph(isnan(ph)) = [];
                                            ph = ph./abs(ph);
                                            plv(curFreq) = abs(nansum(ph))/length(ph);
                                        end;
                                        mX(curMW2) = max( plv );           
                                        
                                        [ plvRand ] = NaN( length( fIx), nRand );
                                        parfor randIter = 1:nRand
                                            dum = NaN(1,length(fIx));
                                            rIx = randperm(phiTrl);
                                            for curFreq = 1:length( fIx)
                                                ph = squeeze(tmp(:,curFreq,:))';
                                                ph = ph(:,rIx);
                                                ph = ph(:);
                                                ph = ph(ix);
                                                ph(isnan(ph)) = [];
                                                ph = ph./abs(ph);
                                                dum(curFreq) = abs(nansum(ph))/length(ph);
                                            end;
                                            plvRand(:,randIter) = dum;
                                        end;
                                        [mxRnd(curMW2,:)] = max( plvRand,[],1 );
                                        
                                    case 'ppc'
                                        N = length(ix);
                                        p = zeros(N*(N-1)/2,2);
                                        cnt = 0;
                                        for curTrl = 1:N-1
                                            for curTrl2 = curTrl+1:N
                                                cnt= cnt+1;
                                                p(cnt,:) = [curTrl curTrl2];
                                            end;
                                        end;
                                        
                                        [ ppc ] = NaN( 1,length(fIx) );
                                        parfor curFreq = 1:length( fIx )
                                            ph = squeeze(tmp(:,curFreq,:))';
                                            ph = ph(:);
                                            ph = ph(ix);
                                            d = diff(angle(ph(p)),[],2);
                                            ppc(curFreq) = nansum(cos(d))/size(p,1);
                                        end;
                                        mX(curMW2) = max( ppc );
                                        
                                        [ ppcRand ] = NaN( length( fIx), nRand );
                                        parfor randIter = 1:nRand
                                            dum = NaN(1,length(fIx));
                                            rIx = randperm(phiTrl);
                                            for curFreq = 1:length( fIx)
                                                ph = squeeze(tmp(:,curFreq,:))';
                                                ph = ph(:,rIx);
                                                ph = ph(:);
                                                ph = ph(ix);
                                                d = diff(angle(ph(p)),[],2);
                                                dum(curFreq) = nansum(cos(d))/size(p,1);
                                            end;
                                            ppcRand(:,randIter) = dum;
                                        end;
                                        [mxRnd(curMW2,:)] = max( ppcRand,[],1 );                                        
                                end
                                                              
                            end;                            
                            [spk2lfp(curMW,:)] = mX;            mX      = [];
                            [spk2lfpRnd(curMW,:,:)] = mxRnd;    mxRnd   = [];
                            
                        end;
                        fprintf('\n');
                    end;
                    
                    %%                    
                    [ Pval ] = zeros(size(spk2lfp));
                    for curMW = 1:size(spk2lfp,1)
                        for curMW2 = 1:size(spk2lfp,2)
                            Pval(curMW,curMW2) = length(find(spk2lfpRnd(curMW,curMW2,:) >= spk2lfp(curMW,curMW2)))/nRand;
                        end;
                    end;
                    clear spk2lfpRnd spk2lfp;
                    [sigIxLFP,sigIxSPK] = find( Pval < (0.05/length(Pval(:))) )
                    clear Pval;
                    
                    %%                    
                    %[ spk2LFPCoupling ] = NaN( length(sigIxLFP),  length(sigIxSPK), length( fIx) );
                    [ spk2LFPCoupling ] = NaN( length(sigIxLFP), length( fIx) );
                    for curMW = 1:length(sigIxLFP)
                        fprintf([num2str(curMW),'/',num2str(length(sigIxLFP))]);
                        
                        if ~isempty(phi{sigIxLFP(curMW)})
                            [ tmp ] = squeeze(phi{sigIxLFP(curMW)}.fourierspctrm(:,:,:,tIx));
                            
                            curMW2 = curMW;
                            %for curMW2 = 1:length( sigIxSPK )
                                ts = spkTs{sigIxSPK(curMW2)}(tIx,:);
                                ts(:,delIx{sigIxLFP(curMW)}) = [];
                                ix = find(ts==1);
                                
                                switch spk2LFPmode
                                    case 'plv'
                                        parfor curFreq = 1:size( tmp,2 )
                                            ph = squeeze(tmp(:,curFreq,:))';
                                            ph = ph(:);
                                            ph = ph(ix);
                                            ph(isnan(ph)) = [];
                                            ph = ph./abs(ph);
                                            %spk2LFPCoupling(curMW,curMW2,curFreq) = abs(nansum(ph))/length(ph);
                                            spk2LFPCoupling(curMW,curFreq) = abs(nansum(ph))/length(ph);
                                        end;
                                    case 'ppc'
                                        N = length(ix);
                                        p = zeros(N*(N-1)/2,2);
                                        cnt = 0;
                                        for curTrl = 1:N-1
                                            for curTrl2 = curTrl+1:N
                                                cnt= cnt+1;
                                                p(cnt,:) = [curTrl curTrl2];
                                            end;
                                        end;
                                        
                                        [ ppc ] = NaN( 1,length(fIx) );
                                        parfor curFreq = 1:length( fIx )
                                            ph = squeeze(tmp(:,curFreq,:))';
                                            ph = ph(:);
                                            ph = ph(ix);
                                            d = diff(angle(ph(p)),[],2);
                                            %spk2LFPCoupling(curMW,curMW2,curFreq) = nansum(cos(d))/size(p,1);
                                            spk2LFPCoupling(curMW,curFreq) = nansum(cos(d))/size(p,1);
                                        end;
                                end;
                                
                            %end;
                        end;
                        fprintf('\n');
                    end;
                    clear phi tmp ph d ix cnt p N ts spkTs;
                    
                    [ chanLabLFP ] = chanLabLFP(sigIxLFP);
                    [ chanLabSPK ] = chanLabSPK(sigIxSPK);                    
                    %[ cluLab ]     = cluLab(sigIxSPK); 
                    
                    %%
                    saveName =[pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spk2LFPCouplingII_',spk2LFPmode,'_lowFreq.mat'];
                    %saveName =[pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spk2LFPCouplingSpikeSorting_',spk2LFPmode,'_lowFreq.mat'];
                    save([savePath,saveName],'spk2LFPfreqAx','spk2LFPCoupling','chanLabLFP','chanLabSPK','sigIxLFP','sigIxSPK','-v7.3');
                    clear spk2LFPfreqAx spk2LFPCoupling chanLabLFP chanLabSPK sigIxLFP sigIxSPK;
                    
                end;
            end;
        end;
    end;
end;



