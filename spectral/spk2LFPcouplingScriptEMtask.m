%%
addpath('~/prj/Bham/code/mcode/utils/');
addpath('~/tbx/CircStat2012a/');
addpath(genpath('/home/rouxf/tbx/eeglab14_1_1b/functions/sigprocfunc/'));
addpath(genpath('~/tbx/chronux_2_11/spectral_analysis/'));
addpath('~/tbx/fieldtrip-20170618/');
ft_defaults;

% %%
% paramsPowS                  = [];
% paramsPowS.Fs               = 1e3; 
% paramsPowS.pad              = 0;
% paramsPowS.fpass                = [3 12];
% paramsPowS.tapers               = [2 3];
% paramsPowS.trialave             = 1;

%%
pId = {'P02','P04','P05','P07','P08','P09','P03ERL'};%,'P22AMS','P23AMS'};%{'P05'};%,'P09'{'P09'};%
expMode = {'fVSpEM'};%{'fVSpEM','cnEM'};%
savePath = '~/resultsSpikeFieldOct18/';

%%
if isempty( gcp('Nocreate') ) 
    parpool(36,'SpmdEnabled',false);
end;

%%
for curPat = 1:length(pId)
    for curExp = 1:length(expMode)
        
        tmp = dir(['~/MICRO/pool/',pId{curPat},'_',expMode{curExp},'_*_lfpDataStimLockedSegmenteddownsampled.mat']);
        %tmp = dir(['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pId{curPat},'/',expMode{curExp},'/']);
        
        if ~isempty(tmp)
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
            for curSesh = 1:length(sesh)
                
                p2d = ['~/MICRO/pool/'];
                %p2d = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pId{curPat},'/',expMode{curExp},'/',sesh{curSesh},'/'];
                fN = dir([p2d,pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_lfpDataStimLockedSegmenteddownsampled.mat']);
                
                lfpDat = load([p2d,fN.name])
                
%                 %%
%                 movingwin = [.5 .01];
%                 
%                 params                  = [];
%                 params.Fs               = dsFs;
%                 params.pad              = 4;
%                 params.fpass            = [0 30];
%                 params.tapers           = [1 1];
%                 params.trialave         = 0;
%                 
                cfgtf                     = [];
                cfgtf.method              = 'wavelet';
                cfgtf.output              = 'pow';
                cfgtf.pad                 = 'nextpow2';
                cfgtf.foi                 = 0:30;
                cfgtf.width               = 4;
                cfgtf.toi                 = lfpDat.dsTrlTime(1):1/lfpDat.dsFs:lfpDat.dsTrlTime(end);
                
                cfgtf2                     = cfgtf;
                cfgtf2.output              = 'fourier';
                
                %%
                BFlab = cell(length(lfpDat.chanLab),1);
                for curChan = 1:length( lfpDat.chanLab )
                    tmp = lfpDat.chanLab{curChan};
                    ix = regexp(tmp,'\d{1}');
                    BFlab(curChan) = {tmp(1:ix-1)};
                end;
                BFid = unique(BFlab);
                
                BFix = cell(length(BFid),1);
                %refDat = cell(1,length(BFlab));
                for curBF = 1:length(BFid)
                    %selBF = find(strcmp(BFid,BFlab{curBF}));
                    %BFix{curBF} = min(find(strcmp(BFid{curBF},BFlab)));
                    %refDat{curBF} = LFPseg{BFix{curBF}};
                    BFix{curBF} = find(strcmp(BFid{curBF},BFlab));
                end;                                               
                
                %%
                for curBF = 1:length(BFid)
                    rms = zeros(length(BFix{curBF}),1);
                    for curMW = 1:length(BFix{curBF})
                        rms(curMW) = mean(sqrt(mean(lfpDat.LFPseg{BFix{curBF}(curMW)}.^2,2)));
                    end;
                    rms = (rms-mean(rms))./std(rms);
                    lfpDat.LFPseg(abs(rms)>8) = [];
                    lfpDat.chanLab(abs(rms)>8) = [];
                end;
                
                %%
%                 Sxx   = cell(1,length( LFPseg));
%                 Ser   = cell(1,length( LFPseg));
%                 pow   = cell(1,length( LFPseg));
%                 erp   = cell(1,length( LFPseg));
                phi   = cell(1,length( lfpDat.LFPseg));
%                 itc   = cell(1,length( LFPseg));
                delIx = cell(1,length( lfpDat.LFPseg));
%                 refIx = cell(1,length( LFPseg));
                for curChan = 1:length( lfpDat.LFPseg )
                    
                    lfp = lfpDat.LFPseg{curChan}(:,ismember(lfpDat.trlSel,lfpDat.trlENC));%sum([ismember(trlSel,trlENC(hitIdx)) ismember(trlSel,trlRET(hitIdx))],2)==1                                        
                    %lfp = lfp-refDat{curChan}(:,ismember(trlSel,trlENC));
                    
                    rms = sqrt(mean(lfp(lfpDat.dsTrlTime>=-0.5 & lfpDat.dsTrlTime <=5,:).^2,1));
                    rms = ( rms-mean(rms) )./std(rms);
                    
                    m  = ones(size(lfp,1),1)*mean(lfp,1);
                    sd = ones(size(lfp,1),1)*std(lfp,0,1);
                    z = ( lfp-m )./ sd;
                    z = max(abs(z(lfpDat.dsTrlTime>=-0.5 & lfpDat.dsTrlTime <=5,:)),[],1);
                    
                    delIx{curChan} = unique([find(rms > 4) find(z > 4)]);                    
                    lfp(:,delIx{curChan}) = [];
                                                            
                    if (~isempty(lfp)) && (size(lfp,2)>25)
                        
%                         rIx = find(strcmp(BFlab{curChan},BFid));
%                         powS = zeros(length(rIx),147);
%                         for refChan = 1:length( BFix{rIx} )
%                             refDat = LFPseg{BFix{rIx}(refChan)}(:,ismember(trlSel,trlENC));                            
%                             refDat(:,delIx{curChan}) = [];
%                             lfp2 = lfp-refDat;
%                             m  = ones(size(lfp,1),1)*mean(lfp2,1);
%                             sd = ones(size(lfp,1),1)*std(lfp2,0,1);                                                                                    
%                             lfp2 = (lfp2-m)./sd;
%                             [powS(refChan,:),fx] = mtspectrumc( lfp2, paramsPowS );                            
%                         end;
%                         [~,mIx] = max(mean(powS,2));
%                         
%                         refIx{curChan} = mIx;
                        
%                         refDat = LFPseg{BFix{rIx}(mIx)}(:,ismember(trlSel,trlENC));
%                         refDat(:,delIx{curChan}) = [];
%                         lfp = lfp-refDat;
                        m  = ones(size(lfp,1),1)*mean(lfp,1);
                        sd = ones(size(lfp,1),1)*std(lfp,0,1);
                        
%                         erp{curChan} = lfp-m;
                        
                        lfp = (lfp-m);%./sd;
                        
%                         [Sxx{curChan},~,~] = mtspecgramc( gradient(lfp')', movingwin, params);
%                         [Ser{curChan},~,~] = mtspecgramc( gradient(mean(lfp,2)')', movingwin, params);
                        
                        dum                     = [];
                        dum.fsample             = lfpDat.dsFs;
                        dum.label               = {'dumChan1'};
                        dum.trial               = cell(1,size(lfp,2));
                        dum.time                = cell(1,size(lfp,2));
                        for curTrl = 1:size(lfp,2)
                            dum.trial{curTrl}            = gradient( lfp(:,curTrl)' );
                            dum.time{curTrl}             = lfpDat.dsTrlTime;
                        end;
                        
                        %[pow{curChan}] = ft_freqanalysis( cfgtf, dum );
                        
                        [phi{curChan}] = ft_freqanalysis( cfgtf2, dum );
                        
%                         tmp = phi{curChan}.fourierspctrm./abs(phi{curChan}.fourierspctrm);
%                         [ itc{curChan} ] = 1/size(tmp,1)*abs(nansum(tmp,1));
                    end;
                    
                end;
                
                for curChan = 1:length(phi)
                    if ~isempty(phi{curChan})
%                         [~,tx,fx] = mtspecgramc( ones(size(lfpDat.LFPseg{curChan},1),1) , movingwin, params);
                        itcTime = phi{curChan}.time;
                        itcFreq = phi{curChan}.freq;
                        break;
                    end;
                end;
                
                %%
                clear BFid BFix BFlab;
                
                %%
                fN = dir([p2d,pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spkDataStimLockedSegmented.mat'])
                spkDat = load([p2d,fN.name])
                
                %%
                %[ SPKseg ] = cell(1,length(spkDat.sortedSpikesSEG));
                [ spkTs ]  = cell(1,length(spkDat.sortedSpikesSEG));
                
                dt = min(lfpDat.dsTrlTime).*1e3-1:1:max(lfpDat.dsTrlTime).*1e3+1;
                %dt2 = -300:300;
                [ frM ] = [];%zeros(length(sortedSpikesSEG),length(dt));
                %[ frSD ] = [];%zeros(length(sortedSpikesSEG),length(dt));
                %[ xc ] = [];%zeros(length(sortedSpikesSEG),length(dt2));
                [ mSpkCnt ] = [];%zeros(length(sortedSpikesSEG),1);
                cluCnt = 0;
                chIx = [];
                cluLab = [];
                for  curMW = 1:length(spkDat.sortedSpikesSEG )
                    fprintf([num2str(curMW),'/',num2str(length(spkDat.sortedSpikesSEG ))]);
                    cID = unique(spkDat.sortedSpikesSEG{curMW}.assignedClusterSeg);
                    for curClu = 1:length(cID)
                        cluCnt = cluCnt+1;
                        %SPKseg{cluCnt}.trial = cell(1,length( spkDat.trlENC) );
                        spkTs{cluCnt} = zeros(length(lfpDat.dsTrlTime),length(spkDat.trlENC));                        
                        fr = zeros(length(spkDat.trlENC),length(dt));
                        %xcTrl= zeros(length(spkDat.trlENC),length(dt2));
                        cluIx = find( spkDat.sortedSpikesSEG{curMW}.assignedClusterSeg==cID(curClu) );
                        for curTrl = 1:length( spkDat.trlENC )
                            ts = spkDat.sortedSpikesSEG{curMW}.SpikeTimesSeg(spkDat.sortedSpikesSEG{curMW}.trl(cluIx) == spkDat.trlENC( curTrl )).*1e3;
                            %SPKseg{cluCnt}.trial{curTrl} = ts;
                            [n,~] = hist(ts,dt);
                            
                            spkTs{cluCnt}(:,curTrl) = n(2:end-1);
                            fr(curTrl,:) = conv(n,gausswin(251),'same')./0.251;
                            %xcTrl(curTrl,:) = xcorr(n(2:end-1),300);
                        end;
                        mSpkCnt(cluCnt) = median(sum(spkTs{cluCnt}(lfpDat.dsTrlTime >0 & lfpDat.dsTrlTime <=4,:),1));
                        frM(cluCnt,:) = mean(fr,1); clear fr;
                        %frSD(cluCnt,:) = std(fr,0,1)./sqrt(size(fr,1)-1);
                        %xc(cluCnt,:) = mean(xcTrl,1);
                        chIx(cluCnt) = curMW;
                        cluLab = [cluLab;cID(curClu)];
                    end;
                    fprintf('\n');
                end;
                %xc(:,dt2 ==0) = NaN;
                
                %%
                [ spkSelIx ] = [];
                for curMW = 1:length( spkTs )
                    x = spkTs{curMW}(lfpDat.dsTrlTime >=-0.5 & lfpDat.dsTrlTime<=5,:);
                    if ( sum(x(:)) >= 50 ) && (mSpkCnt(curMW)>2) && ( mean(frM(curMW,lfpDat.dsTrlTime >=-0.5 & lfpDat.dsTrlTime<=5)) > 1 )
                        spkSelIx = [spkSelIx curMW];
                    end;
                end; clear frM mSpkCnt x;
                %[xc]           = xc(spkSelIx,:);
                [ spkTs ]      = spkTs(spkSelIx);
                %[ frM ]        = frM(spkSelIx,:);
                %[ frSD ]       = frSD(spkSelIx,:);
                [ chIx ]       = chIx(spkSelIx);   
                [ cluLab ]       = cluLab(spkSelIx);   
                
                [ chanLabSPK ] = lfpDat.chanLab(chIx);clear chIx;
                [ chanLabLFP ] = lfpDat.chanLab;
                
%                 %%
%                 saveName = [pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spkParams.mat'];
%                 save([savePath,saveName],'spkTs','frM','frSD','xc','chanLabSPK','dt','dt2');
                
                %%
                [ spkSelIx2 ] = [];
                for curMW = 1:length( spkTs )
                    x = spkTs{curMW}(lfpDat.dsTrlTime >=-0.5 & lfpDat.dsTrlTime<=5,:);
                    x(:,delIx{curMW}) = [];
                    if ( sum(x(:)) >= 50 )
                        spkSelIx2 = [spkSelIx2 curMW];
                    end;
                end;  clear x;              
                %[xc]           = xc(spkSelIx2,:);
                [ spkTs ]      = spkTs(spkSelIx2);
                %[ frM ]        = frM(spkSelIx2,:);
                %[ frSD ]       = frSD(spkSelIx2,:);
                [ chanLabSPK ] = chanLabSPK(spkSelIx2);                               
                [ cluLab ]       = cluLab(spkSelIx2); 
                
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
                    
%                     [ spk2lfpITC ] = zeros(length(phi),length(spkTs));
%                     %[ spk2lfpPPC ] = zeros(length(phi),length(spkTs));
%                     [ spk2lfpITCrand ] = zeros(length(phi),length(spkTs),nRand);
%                     
%                     for curMW = 1:length(phi)
%                         if ~isempty(phi{curMW})
%                             tIx = find( phi{curMW}.time >=0 & phi{curMW}.time <2 );
%                             fIx = find( phi{curMW}.freq >=4 & phi{curMW}.freq <=11 );
%                             break;
%                         end;
%                     end;
%                     
%                     for curMW = 1:length(phi)
%                         fprintf([num2str(curMW),'/',num2str(length(phi))]);
%                         
%                         if ~isempty(phi{curMW})
%                             [ tmp ] = squeeze(phi{curMW}.fourierspctrm);
%                             [ tmp ] = squeeze(tmp(:,fIx,tIx));
%                             
%                             for curMW2 = 1:length( spkTs )
%                                 ts = spkTs{curMW2}(tIx,:);
%                                 ts(:,delIx{curMW}) = [];
%                                 ix = find(ts==1);
%                                 
%                                 [ itcTMP ] = NaN( 1,length( fIx) );
%                                 parfor curFreq = 1:length( fIx)
%                                     ph = squeeze(tmp(:,curFreq,:))';
%                                     ph = ph(:);
%                                     ph = ph(ix);
%                                     ph(isnan(ph)) = [];
%                                     ph = ph./abs(ph);
%                                     itcTMP(curFreq) = abs(nansum(ph))/length(ph);
%                                 end;
%                                 [spk2lfpITC(curMW,curMW2)] = max( itcTMP );
%                                 
%                                 [ itcTMPrand ] = NaN( length( fIx), nRand );
%                                 parfor randIter = 1:nRand
%                                     dum = NaN(1,length(fIx));
%                                     for curFreq = 1:length( fIx)
%                                         ph = squeeze(tmp(:,curFreq,:))';
%                                         ph = ph(:);                                        
%                                         ph = ph(randperm(length(ph)));
%                                         ph = ph(ix);
%                                         ph(isnan(ph)) = [];
%                                         ph = ph./abs(ph);
%                                         dum(curFreq) = abs(nansum(ph))/length(ph);
%                                     end;
%                                     itcTMPrand(:,randIter) = dum;
%                                 end;
%                                 [spk2lfpITCrand(curMW,curMW2,:)] = max( itcTMPrand,[],1 );
%                                 
%                                 %                             N = length(ix);
%                                 %                             p = zeros(N*(N-1)/2,2);
%                                 %                             cnt = 0;
%                                 %                             for curTrl = 1:N-1
%                                 %                                 for curTrl2 = curTrl+1:N
%                                 %                                     cnt= cnt+1;
%                                 %                                     p(cnt,:) = [curTrl curTrl2];
%                                 %                                 end;
%                                 %                             end;
%                                 %
%                                 %                             [ ppc ] = NaN( 1,length(fIx) );
%                                 %                             parfor curFreq = 1:length( fIx )
%                                 %                                 ph = squeeze(tmp(:,curFreq,:))';
%                                 %                                 ph = ph(:);
%                                 %                                 ph = ph(ix);
%                                 %                                 d = diff(angle(ph(p)),[],2);
%                                 %                                 ppc(curFreq) = nansum(cos(d))/size(p,1);
%                                 %                             end;
%                                 %                             spk2lfpPPC(curMW,curMW2) = max( ppc );
%                                 
%                             end;
%                         end;
%                         fprintf('\n');
%                     end;
%                     
%                     itcPval = zeros(size(spk2lfpITC));
%                     for curMW = 1:size(spk2lfpITC,1)
%                         for curMW2 = 1:size(spk2lfpITC,2)
%                             itcPval(curMW,curMW2) = length(find(spk2lfpITCrand(curMW,curMW2,:) >= spk2lfpITC(curMW,curMW2)))/nRand;
%                         end;
%                     end;
%                     
%                     [plvSigIx1,plvSigIx2] = find( itcPval < (0.05/length(itcPval(:))) );
%                     %[plvSigIx1,plvSigIx2] = find( spk2lfpITC >= 0.35 )
%                     %                 [plvSigIx1,plvSigIx2] = find( spk2lfpPPC >= 0.1 )
%                     
%                     %%
%                     [chanLabLFP] = chanLab(unique(plvSigIx1));
%                     
%                     Sxx2 = Sxx(unique(plvSigIx1));
%                     saveName =[pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spectrogramTrialSegCue_lowFreq.mat'];
%                     save([savePath,saveName],'Sxx2','tx','fx','chanLabLFP','-v7.3');
%                     clear Sxx2;
%                     
%                     Ser2 = Ser(unique(plvSigIx1));
%                     saveName =[pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spectrogramTrialSegERPCue_lowFreq.mat'];
%                     save([savePath,saveName],'Ser2','tx','fx','chanLabLFP');
%                     clear Ser2;
%                     
%                     itc2 = itc(unique(plvSigIx1));
%                     saveName =[pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spectrogramTrialSegITCCue_lowFreq.mat'];
%                     save([savePath,saveName],'itc2','itcTime','itcFreq','chanLabLFP');
%                     clear itc2;
%                     
%                     %%
%                     [xc2]           = xc(unique(plvSigIx2),:);
%                     [ spkTs2 ]      = spkTs(unique(plvSigIx2));
%                     [ frM2 ]        = frM(unique(plvSigIx2),:);
%                     [ frSD2 ]       = frSD(unique(plvSigIx2),:);
%                     [ chanLabSPK2 ] = chanLabSPK(unique(plvSigIx2));
%                 
%                     %%
%                     saveName = [pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spkParamsCueSFC.mat'];
%                     save([savePath,saveName],'spkTs2','frM2','frSD2','xc2','chanLabSPK2','dt','dt2');
%                     
%                     %%                    
%                     paramsC                     = [];
%                     paramsC.pad                 = 4; 
%                     paramsC.Fs                  = 1e3;
%                     paramsC.fpass               = [0 30];
%                     paramsC.tapers              = [3 5]; 
%                     paramsC.trialave            = 1;
%                     
%                     n = length(plvSigIx1);
%                     np = n*(n-1)/2;
%                     p = zeros(np,2);
%                     
%                     cnt = 0;
%                     for curChan = 1:length( plvSigIx1)-1
%                         for curChan2 = curChan+1:length( plvSigIx1)
%                             cnt = cnt+1;
%                             p(cnt,:) = [plvSigIx1(curChan) plvSigIx1(curChan2)];
%                         end;
%                     end;
%                     
%                     [ Cxy ] = cell(1,size(p,1));
%                     for curMW = 1:size(p,1)                        
%                         fprintf([num2str(curMW),'/',num2str(size(p,1))]);
%                         dIx = unique([delIx{p(curMW,1)} delIx{p(curMW,2)}]);
%                         [ lfp1 ] = LFPseg{p(curMW,1)}(dsTrlTime >= 0 & dsTrlTime <2,ismember(trlSel,trlENC));
%                         lfp1(:,dIx) = [];
%                         lfp1 = lfp1-ones(size(lfp1,1),1)*mean(lfp1,1);
%                         [ lfp2 ] = LFPseg{p(curMW,2)}(dsTrlTime >= 0 & dsTrlTime <2,ismember(trlSel,trlENC));                        
%                         lfp2(:,dIx) = [];
%                         lfp2 = lfp2-ones(size(lfp2,1),1)*mean(lfp2,1);
%                         [Cxy{curMW},~,~,~,~,cxyFx] = coherencyc(lfp1,lfp2,paramsC);
%                         fprintf('\n');
%                     end;
%                     
%                     %%
%                     saveName = [pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_LFP2LFPCohCue.mat'];
%                     save([savePath,saveName],'Cxy','cxyFx');
%                     
%                     %%
%                     for curMW = 1:length(phi)
%                         if ~isempty(phi{curMW})
%                             tIx = find( phi{curMW}.time >=0 & phi{curMW}.time <2 );
%                             fIx = find( phi{curMW}.freq >=4 & phi{curMW}.freq <=11 );
%                             freqAx = phi{curMW}.freq;
%                             break;
%                         end;
%                     end;
%                     
%                     movingwinSTA = [0.5 .025];
%                     
%                     paramsSTA                  = [];
%                     paramsSTA.pad              = 6;
%                     paramsSTA.Fs               = dsFs;
%                     paramsSTA.fpass            = [0 30];
%                     paramsSTA.tapers           = [1 1];
%                     paramsSTA.trialave         = 1;
%                     
%                     params                  = [];
%                     params.pad              = 4;
%                     params.Fs               = dsFs;
%                     params.fpass            = [0 30];
%                     params.tapers           = [2 3];
%                     params.trialave         = 1;
%                     
%                     [ zval ] = NaN( length(plvSigIx1),length( freqAx ));
%                     [ sta ] = cell(1,length(plvSigIx1));
%                     [ sta2 ] = cell(1,length(plvSigIx1));
%                     [ staSxx ] = cell(1,length( plvSigIx1));
%                     [ staSxx2 ] = cell(1,length( plvSigIx1));
%                     [ nPA ] = cell(1,length( plvSigIx1));
%                     [ staStf ] = cell(1,length( plvSigIx1));
%                     
%                     for curMW = 1:length(plvSigIx1)
%                         fprintf([num2str(curMW),'/',num2str(length(plvSigIx1))]);
%                         
%                         [ tmp ] = squeeze(phi{plvSigIx1(curMW)}.fourierspctrm(:,:,:,tIx));
%                         [ tmp ] = angle(tmp);
%                         
%                         ts = spkTs{plvSigIx2(curMW)}(tIx,:);
%                         ts(:,delIx{plvSigIx1(curMW)}) = [];
%                         
%                         if any(size(ts) ~= [size(tmp,3) size(tmp,1)])
%                             error(':(');
%                         end;
%                         
%                         [ ix ] = find(ts==1);
%                         
%                         for curFreq = 1:length( freqAx )
%                             ph = squeeze(tmp(:,curFreq,:))';
%                             ph = ph(:);
%                             ph = ph(ix);
%                             ph(isnan(ph)) = [];
%                             if ~isempty(ph)
%                                 [~,zval(curMW,curFreq)] = circ_rtest( ph );
%                             end;
%                         end;
%                         
%                         [~,mIx] = max(zval(curMW,fIx));
%                         ph = squeeze(tmp(:,fIx(mIx),:))';
%                         ph = ph(:);
%                         ph = ph(ix);
%                         ph(isnan(ph)) = [];
%                         nPA{curMW} = ph;
%                         
%                         [ st ] = spkTs{plvSigIx2(curMW)}(tIx,:);
%                         st(:,delIx{plvSigIx1(curMW)}) = [];
%                         [ lfp ] = LFPseg{plvSigIx1(curMW)}(:,ismember(trlSel,trlENC));
% %                         lfp = lfp- LFPseg{refIx{plvSigIx1(curMW)}}(:,ismember(trlSel,trlENC));
%                         lfp(:,delIx{plvSigIx1(curMW)}) = [];
%                         
% %                         m  = ones(size(lfp,1),1)*mean(lfp,1);
% %                         sd = ones(size(lfp,1),1)*std(lfp,0,1);                                               
% %                         lfp = (lfp-m)./sd;
%                         
%                         if any(size(st) ~= size(lfp(tIx,:)))
%                             error(':(');
%                         end;
%                         
%                         cnt = 0;
%                         cnt2 = 0;
%                         for curTrl = 1:size(lfp,2)
%                             spkIx = find(st(:,curTrl)==1);
%                             for curSPK = 1:length( spkIx )
%                                 if ( tIx(spkIx(curSPK))-480>0 ) && ( tIx(spkIx(curSPK))+480<=size(lfp,1) ) 
%                                     cnt = cnt+1;
%                                     dum = lfp(tIx(spkIx(curSPK))-480:tIx(spkIx(curSPK))+480,curTrl);
%                                     sta{curMW}(cnt,:) = dum-mean(dum);
%                                 end;
%                                 if ( tIx(spkIx(curSPK))-1e3>0 ) && ( tIx(spkIx(curSPK))+1e3<=size(lfp,1) ) 
%                                     cnt2 = cnt2+1;
%                                     dum = lfp(tIx(spkIx(curSPK))-1e3:tIx(spkIx(curSPK))+1e3,curTrl);
%                                     sta2{curMW}(cnt2,:) = dum-mean(dum);
%                                 end;
%                             end;
%                         end;
%                                                 
%                         if (isempty(sta{curMW})) || (isempty(sta2{curMW}))
%                             error(':(');
%                         end;
%                         
%                         [staSxx{curMW}] = mtspectrumc( gradient(sta{curMW})',params );
%                         [staSxx2{curMW}] = mtspectrumc( gradient(mean(sta{curMW},1))',params );
%                         
%                         [staStf{curMW}] = mtspecgramc(gradient(sta2{curMW})',movingwinSTA, paramsSTA);
%                         fprintf('\n');
%                     end;
%                     
%                     [rayleighFreq] = freqAx;
%                     [staTime] = -480:480;
%                     [staTime2] = -1e3:1e3;
%                     
%                     [~,staFx] = mtspectrumc( ones(length(staTime),1),params );
%                     [~,txSTA,fxSTA] = mtspecgramc(gradient(ones(length(staTime2),1))',movingwinSTA, paramsSTA);
%                     
%                     txSTA = txSTA - 1;
%                     
%                     %%
%                     [ chanLabLFP ] = chanLab(plvSigIx1);
%                     
%                     saveName =[pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_STACue.mat'];
%                     save([savePath,saveName],'sta','staTime','chanLabLFP','chanLabSPK2','plvSigIx1','plvSigIx2');
%                     saveName =[pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_rayleighZvalCue_lowFreq.mat'];
%                     save([savePath,saveName],'zval','rayleighFreq','nPA','chanLabLFP','chanLabSPK2','plvSigIx1','plvSigIx2');
%                     saveName =[pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_STApowCue_lowFreq.mat'];
%                     save([savePath,saveName],'staSxx','staFx','chanLabLFP','chanLabSPK2','plvSigIx1','plvSigIx2');
%                     saveName =[pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_STApowCueAVG_lowFreq.mat'];
%                     save([savePath,saveName],'staSxx2','staFx','chanLabLFP','chanLabSPK2','plvSigIx1','plvSigIx2');
%                     saveName =[pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_STASpectrogramCue_lowFreq.mat'];
%                     save([savePath,saveName],'staStf','txSTA','fxSTA','chanLabLFP','chanLabSPK2','plvSigIx1','plvSigIx2');
                    
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
                    [ spk2LFPCoupling ] = NaN( length(sigIxLFP),  length( fIx) );
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
                    [ cluLab ]     = cluLab(sigIxSPK); 
                    
                    %%
                    %saveName =[pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spk2LFPCoupling_',spk2LFPmode,'_lowFreq.mat'];
                    saveName =[pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spk2LFPCouplingSpikeSortingII_',spk2LFPmode,'_lowFreq.mat'];
                    save([savePath,saveName],'spk2LFPfreqAx','spk2LFPCoupling','chanLabLFP','chanLabSPK','cluLab','sigIxLFP','sigIxSPK','-v7.3');
                    clear spk2LFPfreqAx spk2LFPCoupling chanLabLFP chanLabSPK sigIxLFP sigIxSPK;

                    %[plvSigIx3,plvSigIx4] = find( spk2lfpITC >= 0.35 )
                    %                 [plvSigIx3,plvSigIx4] = find( spk2lfpPPC >= 0.1 )
                    
%                     %%                   
%                     movingwin = [.5 0.01];
%                     
%                     params                  = [];
%                     params.Fs               = 1e3;
%                     params.pad              = 2;
%                     params.fpass            = [0 30];
%                     params.tapers           = [2 3];
%                     params.trialave         = 1;
%                     
%                     [ sFC ] = cell(1,length(plvSigIx3));
%                     for curChan = 1:length( plvSigIx3 )
%                         [ lfpDat ] = LFPseg{plvSigIx3(curChan)}(:,ismember(trlSel,trlENC));
%                         [ spkDat ] = spkTs{plvSigIx4(curChan)};
%                         
%                         lfpDat(:,delIx{plvSigIx3(curChan)}) = [];
%                         lfpDat = lfpDat - ones(length(dsTrlTime),1)*mean(lfpDat,1);
%                         
%                         spkDat(:,delIx{plvSigIx3(curChan)}) = [];
%                         
%                         [sFC{curChan}] = cohgramcpb(lfpDat,spkDat,movingwin,params);
%                     end;
%                     [~,~,~,~,~,txSFC,fxSFC] = cohgramcpb(lfpDat,spkDat,movingwin,params);
%                     
%                     [chanLabLFP] = chanLab(unique(plvSigIx3));
%                     
%                     saveName =[pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_SFCcohgramTrialSegEnc_lowFreq.mat'];
%                     save([savePath,saveName],'sFC','txSFC','fxSFC','chanLabLFP','-v7.3');
                    
%                     %%
%                     [chanLabLFP] = chanLab(unique(plvSigIx3));
%                     
%                     Sxx2 = Sxx(unique(plvSigIx3));
%                     saveName =[pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spectrogramTrialSegEnc_lowFreq.mat'];
%                     save([savePath,saveName],'Sxx2','tx','fx','chanLabLFP','-v7.3');
%                     clear Sxx2;
%                     
%                     Ser2 = Ser(unique(plvSigIx3));
%                     saveName =[pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spectrogramTrialSegERPEnc_lowFreq.mat'];
%                     save([savePath,saveName],'Ser2','tx','fx','chanLabLFP');
%                     clear Ser2;
%                     
%                     itc2 = itc(unique(plvSigIx3));
%                     saveName =[pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spectrogramTrialSegITCEnc_lowFreq.mat'];
%                     save([savePath,saveName],'itc2','itcTime','itcFreq','chanLabLFP');
%                     clear itc2;
%                     
%                      %%
%                     [xc2]           = xc(unique(plvSigIx4),:);
%                     [ spkTs2 ]      = spkTs(unique(plvSigIx4));
%                     [ frM2 ]        = frM(unique(plvSigIx4),:);
%                     [ frSD2 ]       = frSD(unique(plvSigIx4),:);
%                     [ chanLabSPK2 ] = chanLabSPK(unique(plvSigIx4));
%                 
%                     %%
%                     saveName = [pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spkParamsEncSFC.mat'];
%                     save([savePath,saveName],'spkTs2','frM2','frSD2','xc2','chanLabSPK2','dt','dt2');
%                     
%                     %%                    
%                     paramsC                     = [];
%                     paramsC.pad                 = 2; 
%                     paramsC.Fs                  = 1e3;
%                     paramsC.fpass               = [0 30];
%                     paramsC.tapers              = [5 9]; 
%                     paramsC.trialave            = 1;
%                     
%                     n = length(plvSigIx3);
%                     np = n*(n-1)/2;
%                     p = zeros(np,2);
%                     
%                     cnt = 0;
%                     for curChan = 1:length( plvSigIx3)-1
%                         for curChan2 = curChan+1:length( plvSigIx3)
%                             cnt = cnt+1;
%                             p(cnt,:) = [plvSigIx3(curChan) plvSigIx3(curChan2)];
%                         end;
%                     end;
%                     
%                     [ Cxy ] = cell(1,size(p,1));
%                     for curMW = 1:size(p,1)                        
%                         fprintf([num2str(curMW),'/',num2str(size(p,1))]);
%                         dIx = unique([delIx{p(curMW,1)} delIx{p(curMW,2)}]);
%                         [ lfp1 ] = LFPseg{p(curMW,1)}(dsTrlTime >= 0 & dsTrlTime <=5,ismember(trlSel,trlENC));
%                         lfp1(:,dIx) = [];
%                         lfp1 = lfp1-ones(size(lfp1,1),1)*mean(lfp1,1);
%                         [ lfp2 ] = LFPseg{p(curMW,2)}(dsTrlTime >= 0 & dsTrlTime <=5,ismember(trlSel,trlENC));                        
%                         lfp2(:,dIx) = [];
%                         lfp2 = lfp2-ones(size(lfp2,1),1)*mean(lfp2,1);
%                         [Cxy{curMW},~,~,~,~,cxyFx] = coherencyc(lfp1,lfp2,paramsC);
%                         fprintf('\n');
%                     end;
%                     
%                     %%
%                     saveName = [pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_LFP2LFPCohEnc.mat'];
%                     save([savePath,saveName],'Cxy','cxyFx');
%                     
%                     %%
%                     for curMW = 1:length(phi)
%                         if ~isempty(phi{curMW})
%                             tIx = find( phi{curMW}.time >=0 & phi{curMW}.time <=5 );
%                             fIx = find( phi{curMW}.freq >=4 & phi{curMW}.freq <=18 );
%                             freqAx = phi{curMW}.freq;
%                             break;
%                         end;
%                     end;
%                     
%                     movingwinSTA = [0.5 .025];
%                     
%                     paramsSTA                  = [];
%                     paramsSTA.pad              = 6;
%                     paramsSTA.Fs               = dsFs;
%                     paramsSTA.fpass            = [0 30];
%                     paramsSTA.tapers           = [1 1];
%                     paramsSTA.trialave         = 1;
%                     
%                     params                  = [];
%                     params.pad              = 4;
%                     params.Fs               = dsFs;
%                     params.fpass            = [0 30];
%                     params.tapers           = [2 3];
%                     params.trialave         = 1;
%                     
%                     [ zval ] = NaN( length(plvSigIx3),length( freqAx ));
%                     [ sta ] = cell(1,length(plvSigIx3));
%                     [ sta2 ] = cell(1,length(plvSigIx3));
%                     [ staSxx ] = cell(1,length( plvSigIx3));
%                     [ staSxx2 ] = cell(1,length( plvSigIx3));
%                     [ nPA ] = cell(1,length( plvSigIx3));
%                     [ staStf ] = cell(1,length( plvSigIx3));
%                     
%                     for curMW = 1:length(plvSigIx3)
%                         fprintf([num2str(curMW),'/',num2str(length(plvSigIx3))]);
%                         
%                         [ tmp ] = squeeze(phi{plvSigIx3(curMW)}.fourierspctrm(:,:,:,tIx));
%                         [ tmp ] = angle(tmp);
%                         
%                         ts = spkTs{plvSigIx4(curMW)}(tIx,:);
%                         ts(:,delIx{plvSigIx3(curMW)}) = [];
%                         
%                         if any(size(ts) ~= [size(tmp,3) size(tmp,1)])
%                             error(':(');
%                         end;
%                         
%                         [ ix ] = find(ts==1);
%                         
%                         for curFreq = 1:length( freqAx )
%                             ph = squeeze(tmp(:,curFreq,:))';
%                             ph = ph(:);
%                             ph = ph(ix);
%                             ph(isnan(ph)) = [];
%                             if ~isempty(ph)
%                                 [~,zval(curMW,curFreq)] = circ_rtest( ph );
%                             end;
%                         end;
%                         
%                         [~,mIx] = max(zval(curMW,fIx));
%                         ph = squeeze(tmp(:,fIx(mIx),:))';
%                         ph = ph(:);
%                         ph = ph(ix);
%                         ph(isnan(ph)) = [];
%                         nPA{curMW} = ph;
%                         
%                         [ st ] = spkTs{plvSigIx4(curMW)}(tIx,:);
%                         st(:,delIx{plvSigIx3(curMW)}) = [];
%                         [ lfp ] = LFPseg{plvSigIx3(curMW)}(:,ismember(trlSel,trlENC));
% %                         lfp = lfp- LFPseg{refIx{plvSigIx3(curMW)}}(:,ismember(trlSel,trlENC));
%                         %lfp = lfp- refDat{plvSigIx3(curMW)}(:,ismember(trlSel,trlENC));
%                         lfp(:,delIx{plvSigIx3(curMW)}) = [];
%                         
% %                         m  = ones(size(lfp,1),1)*mean(lfp,1);
% %                         sd = ones(size(lfp,1),1)*std(lfp,0,1);                                               
% %                         lfp = (lfp-m)./sd;
%                         
%                         if any(size(st) ~= size(lfp(tIx,:)))
%                             error(':(');
%                         end;
%                         
%                         cnt = 0;
%                         cnt2 = 0;
%                         for curTrl = 1:size(lfp,2)
%                             spkIx = find(st(:,curTrl)==1);
%                             for curSPK = 1:length( spkIx )
%                                 if ( tIx(spkIx(curSPK))-480>0 ) && ( tIx(spkIx(curSPK))+480<=size(lfp,1) ) 
%                                     cnt = cnt+1;
%                                     dum = lfp(tIx(spkIx(curSPK))-480:tIx(spkIx(curSPK))+480,curTrl);
%                                     sta{curMW}(cnt,:) = dum-mean(dum);
%                                 end;
%                                 if ( tIx(spkIx(curSPK))-1e3>0 ) && ( tIx(spkIx(curSPK))+1e3<=size(lfp,1) ) 
%                                     cnt2 = cnt2+1;
%                                     dum = lfp(tIx(spkIx(curSPK))-1e3:tIx(spkIx(curSPK))+1e3,curTrl);
%                                     sta2{curMW}(cnt2,:) = dum-mean(dum);
%                                 end;
%                             end;
%                         end;
%                         
%                         if ~(isempty(sta{curMW})) || ~(isempty(sta2{curMW}))
%                             [staSxx{curMW}] = mtspectrumc( gradient(sta{curMW})',params );
%                             [staSxx2{curMW}] = mtspectrumc( gradient(mean(sta{curMW},1))',params );
%                             
%                             [staStf{curMW}] = mtspecgramc(gradient(sta2{curMW})',movingwinSTA, paramsSTA);
%                         end;
%                         
% 
%                         fprintf('\n');
%                     end;
%                     
%                     [rayleighFreq] = freqAx;
%                     [staTime] = -480:480;
%                     [staTime2] = -1e3:1e3;
%                     
%                     [~,staFx] = mtspectrumc( ones(length(staTime),1),params );
%                     [~,txSTA,fxSTA] = mtspecgramc(gradient(ones(length(staTime2),1))',movingwinSTA, paramsSTA);
%                     
%                     txSTA = txSTA - 1;
%                     
%                     %%
%                     [ chanLabLFP ] = chanLab(plvSigIx3);
%                     
%                     saveName =[pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_STAEnc.mat'];
%                     save([savePath,saveName],'sta','staTime','chanLabLFP','chanLabSPK','plvSigIx3','plvSigIx4');
%                     saveName =[pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_rayleighZvalEnc_lowFreq.mat'];
%                     save([savePath,saveName],'zval','rayleighFreq','nPA','chanLabLFP','chanLabSPK','plvSigIx3','plvSigIx4');
%                     saveName =[pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_STApowEnc_lowFreq.mat'];
%                     save([savePath,saveName],'staSxx','staFx','chanLabLFP','chanLabSPK','plvSigIx3','plvSigIx4');
%                     saveName =[pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_STApowEncAVG_lowFreq.mat'];
%                     save([savePath,saveName],'staSxx2','staFx','chanLabLFP','chanLabSPK','plvSigIx3','plvSigIx4');
%                     saveName =[pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_STASpectrogramEnc_lowFreq.mat'];
%                     save([savePath,saveName],'staStf','txSTA','fxSTA','chanLabLFP','chanLabSPK','plvSigIx3','plvSigIx4');
                    
                end;
            end;
        end;
    end;
end;



