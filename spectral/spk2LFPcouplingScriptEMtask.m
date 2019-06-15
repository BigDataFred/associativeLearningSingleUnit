%%
addpath('~/prj/Bham/code/mcode/utils/');
addpath('~/tbx/CircStat2012a/');
addpath(genpath('~/tbx/chronux_2_11/spectral_analysis/'));
addpath(genpath('~/associativeLearningSingleUnit/'));
addpath('~/tbx/fieldtrip-20170618/');
ft_defaults;

%%
[ pId ] = {'P02','P04','P05','P07','P08','P09','P03ERL','P22AMS','P23AMS'};%{'P05'};%,'P09'{'P09'};%
[ expMode ] = {'fVSpEM','cnEM'};
[ spkMode ] = {'noSorting', 'Sorting'};
[ spk2LFPmode ] = {'plv','ppc'};

%%
[ nRand ] = 200;
[ alpha ] = 0.05;

curSpkSortingMode = 2;
curSpk2LFPmode = 1;

%%
[p2SPKd]       = '~/resultsSpikeFieldOct18II/';
[ savePath ] = '~/resultsSpikeFieldOct18II/';

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
                fN = dir([p2SPKd,pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spkParams_',spkMode{curSpkSortingMode},'.mat'])
                [ spkDat ] = load([p2SPKd,fN.name])
                
                %%
                [ spkTs ] = spkDat.spkTs;
                [ chIx ] = 1:length(spkTs);
                
                [ spkSelIx ] = [];
                for curMW = 1:length( spkTs )
                    
                    x = spkTs{curMW}(lfpDat.dsTrlTime >=-0.5 & lfpDat.dsTrlTime<=5,:);
                    mSpkCnt = median(sum(spkTs{curMW}(lfpDat.dsTrlTime >0 & lfpDat.dsTrlTime <=4,:),1));
                    
                    if ( sum(x(:)) >= 50 ) && ( mSpkCnt>2 ) && ( mean(spkDat.frM(curMW,lfpDat.dsTrlTime >=-0.5 & lfpDat.dsTrlTime<=5)) > 1 )
                        spkSelIx = [spkSelIx curMW];
                    end;
                end;
                
                [ spkTs ]      = spkTs(spkSelIx);
                [ chIx ]       = chIx(spkSelIx);
                [ chanLabSPK ] = lfpDat.chanLab(chIx);
                clear chIx;
                
                %%
                [ spkSelIx2 ] = [];
                for curMW = 1:length( spkTs )
                    x = spkTs{curMW}(lfpDat.dsTrlTime >=-0.5 & lfpDat.dsTrlTime<=5,:);
                    x(:,delIx{curMW}) = [];
                    if ( sum(x(:)) >= 50 )
                        spkSelIx2 = [spkSelIx2 curMW];
                    end;
                end;
                clear x;
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
                if any(chck~=0) && any(chck2~=0)
                    
                    %%
                    for curMW = 1:length(phi)
                        if ~isempty(phi{curMW})
                            tIx = find( phiTime >=-.5 & phiTime < 5 );
                            fIx = find( phiFreq >=4 & phiFreq <=20 );
                            spk2LFPfreqAx = phiFreq;
                            break;
                        end;
                    end;
                    
                    %%
                    [sigIxLFP,sigIxSPK] = computeSpk2LFPCouplingSig( phi, spkTs, nRand, fIx, tIx, delIx,  spk2LFPmode{curSpk2LFPmode}, alpha );
                    
                    %%
                    [ spk2LFPCoupling ] = NaN( length(sigIxLFP), length( fIx) );
                    for curMW = 1:length(sigIxLFP)
                        fprintf([num2str(curMW),'/',num2str(length(sigIxLFP))]);
                        if ~isempty(phi{sigIxLFP(curMW)})
                            [ tmp ] = squeeze(phi{sigIxLFP(curMW)}(:,:,:,tIx));
                            
                            curMW2 = curMW;
                            ts = spkTs{sigIxSPK(curMW2)}(tIx,:);
                            ts(:,delIx{sigIxLFP(curMW)}) = [];
                            ix = find(ts==1);
                            
                            [spk2LFPCoupling(curMW,:)] = computeSPK2LFPcoupling( phi, spkIx, fIx, spk2LFPmode{curSpk2LFPmode} );                                                                                    
                        end;
                        fprintf('\n');
                    end;
                    
                    %%
                    [ chanLabLFP ] = chanLabLFP(sigIxLFP);
                    [ chanLabSPK ] = chanLabSPK(sigIxSPK);
                    
                    %%
                    saveName =[pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spk2LFPCouplingII_',spk2LFPmode{curSpk2LFPmode},'_',spkMode{curSpkSortingMode},'_lowFreq.mat'];
                    save([savePath,saveName],'spk2LFPfreqAx','spk2LFPCoupling','chanLabLFP','chanLabSPK','sigIxLFP','sigIxSPK','-v7.3');
                    clear spk2LFPfreqAx spk2LFPCoupling chanLabLFP chanLabSPK sigIxLFP sigIxSPK;
                    
                end;
            end;
        end;
    end;
end;



