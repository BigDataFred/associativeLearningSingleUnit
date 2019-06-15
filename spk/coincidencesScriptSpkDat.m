%%
addpath('~/prj/Bham/code/mcode/utils/');
addpath('~/tbx/CircStat2012a/');
addpath(genpath('/home/rouxf/tbx/eeglab14_1_1b/functions/sigprocfunc/'));
addpath(genpath('~/tbx/chronux_2_11/spectral_analysis/'));
addpath('~/tbx/fieldtrip-20170618/');
ft_defaults;

%%
pId = {'P02','P04','P05','P07','P08','P09','P03ERL'};
expMode = {'fVSpEM'};%{'fVSpEM','cnEM'};%
savePath = '~/resultsSpikeFieldOct18/';

p2d = '/home/rouxf/MICRO/pool/';
p2d2 = '/home/rouxf/resultsSpikeFieldOct18/';

%%
if isempty( gcp('Nocreate') )
    parpool(36);
end;

%%
coFI1 = [];
coFI2 = [];
coFIH = [];
coFIM = [];
poolCnt1 = 0;
poolCnt2 = 0;

nPLVunits   = 0;
nNoPLVunits = 0;

%%
for curPat = 1:length(pId)
    for curExp = 1:length(expMode)
        
        %%
        tmp = dir(['/home/rouxf/MICRO/pool/',pId{curPat},'_',expMode{curExp},'_*_lfpDataStimLockedSegmenteddownsampled.mat']);
        
        %%
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
                
                %%
                fN = dir([p2d,pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_lfpDataStimLockedSegmenteddownsampled.mat']);
                [ lfpDat ] = load([p2d,fN.name])
                
                %%
                fN = dir([p2d,pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spkDataStimLockedSegmented.mat'])
                [ spkDat ] = load([p2d,fN.name])
                
                %%
                fNspk2LFP = [pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spk2LFPCouplingSpikeSorting_plv_lowFreq.mat'];
                fNspk2LFP = dir([p2d2,fNspk2LFP]);
                [ plvSigDat ] = load([p2d2,fNspk2LFP.name]);
                
                %                                 %%
                %                 BFlab = cell(length(lfpDat.chanLab),1);
                %                 for curChan = 1:length( lfpDat.chanLab )
                %                     tmp = lfpDat.chanLab{curChan};
                %                     ix = regexp(tmp,'\d{1}');
                %                     BFlab(curChan) = {tmp(1:ix-1)};
                %                 end;
                %                 BFid = unique(BFlab);
                %
                %                 BFix = cell(length(BFid),1);
                %                 for curBF = 1:length(BFid)
                %                     BFix{curBF} = find(strcmp(BFid{curBF},BFlab));
                %                 end;
                %
                %                 %%
                %                 for curBF = 1:length(BFid)
                %                     rms = zeros(length(BFix{curBF}),1);
                %                     for curMW = 1:length(BFix{curBF})
                %                         rms(curMW) = mean(sqrt(mean(lfpDat.LFPseg{BFix{curBF}(curMW)}.^2,2)));
                %                     end;
                %                     rms = (rms-mean(rms))./std(rms);
                %                     lfpDat.LFPseg(abs(rms)>8) = [];
                %                     lfpDat.chanLab(abs(rms)>8) = [];
                %                 end;
                
                %                 %%
                %                 cfgtf                     = [];
                %                 cfgtf.method              = 'wavelet';
                %                 cfgtf.output              = 'pow';
                %                 cfgtf.pad                 = 'nextpow2';
                %                 cfgtf.foi                 = 0:30;
                %                 cfgtf.width               = 4;
                %                 cfgtf.toi                 = lfpDat.dsTrlTime(1):1/lfpDat.dsFs:lfpDat.dsTrlTime(end);
                %
                %                 cfgtf2                     = cfgtf;
                %                 cfgtf2.output              = 'fourier';
                
                %%
                %                 phi   = cell(1,length( lfpDat.LFPseg));
                %                 delIx = cell(1,length( lfpDat.LFPseg));
                %                 for curChan = 1:length( lfpDat.LFPseg )
                %
                %                     lfp = lfpDat.LFPseg{curChan}(:,ismember(lfpDat.trlSel,lfpDat.trlENC));%sum([ismember(trlSel,trlENC(hitIdx)) ismember(trlSel,trlRET(hitIdx))],2)==1
                %
                %                     rms = sqrt(mean(lfp(lfpDat.dsTrlTime>=-0.5 & lfpDat.dsTrlTime <=5,:).^2,1));
                %                     rms = ( rms-mean(rms) )./std(rms);
                %
                %                     m  = ones(size(lfp,1),1)*mean(lfp,1);
                %                     sd = ones(size(lfp,1),1)*std(lfp,0,1);
                %                     z = ( lfp-m )./ sd;
                %                     z = max(abs(z(lfpDat.dsTrlTime>=-0.5 & lfpDat.dsTrlTime <=5,:)),[],1);
                %
                %                     delIx{curChan} = unique([find(rms > 4) find(z > 4)]);
                %                     lfp(:,delIx{curChan}) = [];
                %
                %                     if (~isempty(lfp)) && (size(lfp,2)>25)
                %
                %                         m  = ones(size(lfp,1),1)*mean(lfp,1);
                %                         sd = ones(size(lfp,1),1)*std(lfp,0,1);
                %
                %                         lfp = (lfp-m);%./sd;
                %
                %                         dum                     = [];
                %                         dum.fsample             = lfpDat.dsFs;
                %                         dum.label               = {'dumChan1'};
                %                         dum.trial               = cell(1,size(lfp,2));
                %                         dum.time                = cell(1,size(lfp,2));
                %                         for curTrl = 1:size(lfp,2)
                %                             dum.trial{curTrl}            = gradient( lfp(:,curTrl)' );
                %                             dum.time{curTrl}             = lfpDat.dsTrlTime;
                %                         end;
                %
                %                         [phi{curChan}] = ft_freqanalysis( cfgtf2, dum );
                %                     end;
                %                 end;
                
                %%
                %                 clear BFid BFix BFlab;
                
                %%
                [ spkTs ]  = cell(1,length(spkDat.sortedSpikesSEG));
                
                dt = min(lfpDat.dsTrlTime).*1e3-1:1:max(lfpDat.dsTrlTime).*1e3+1;
                dt2 = min(lfpDat.dsTrlTime).*1e3-1:1:max(lfpDat.dsTrlTime).*1e3+1;
                [ frM ] = [];%zeros(length(sortedSpikesSEG),length(dt));
                [ mSpkCnt ] = [];%zeros(length(sortedSpikesSEG),1);
                cluCnt = 0;
                chIx = [];
                cluLab = [];
                for  curMW = 1:length(spkDat.sortedSpikesSEG )
                    fprintf([num2str(curMW),'/',num2str(length(spkDat.sortedSpikesSEG ))]);
                    cID = unique(spkDat.sortedSpikesSEG{curMW}.assignedClusterSeg);
                    for curClu = 1:length(cID)
                        cluCnt = cluCnt+1;
                        spkTs{cluCnt} = zeros(length(dt(2:end-1)),length(spkDat.trlENC));
                        fr = zeros(length(spkDat.trlENC),length(dt2(2:end-1)));
                        cluIx = find( spkDat.sortedSpikesSEG{curMW}.assignedClusterSeg==cID(curClu) );
                        for curTrl = 1:length( spkDat.trlENC )
                            ts = spkDat.sortedSpikesSEG{curMW}.SpikeTimesSeg(spkDat.sortedSpikesSEG{curMW}.trl(cluIx) == spkDat.trlENC( curTrl )).*1e3;
                            [n,~] = hist(ts,dt);
                            [n2,~] = hist(ts,dt2);
                            spkTs{cluCnt}(:,curTrl) = n(2:end-1);
                            fr(curTrl,:) = conv(n2(2:end-1),gausswin(251),'same')./0.251;
                        end;
                        mSpkCnt(cluCnt) = median(sum(spkTs{cluCnt}(dt(2:end-1) >0 & dt(2:end-1) <=4e3,:),1));
                        frM(cluCnt,:) = mean(fr,1); clear fr;
                        chIx(cluCnt) = curMW;
                        cluLab = [cluLab;cID(curClu)];
                    end;
                    fprintf('\n');
                end;
                
                %%
                [ spkSelIx ] = [];
                for curMW = 1:length( spkTs )
                    x = spkTs{curMW}(dt(2:end-1) >=-500 & dt(2:end-1)<=5e3,:);
                    if ( sum(x(:)) > 50 ) && (mSpkCnt(curMW)>2) && ( mean(frM(curMW,dt2(2:end-1) >=-500 & dt2(2:end-1)<=5e3)) > 1 )
                        spkSelIx = [spkSelIx curMW];
                    end;
                end; clear frM mSpkCnt x;
                [ spkTs ]      = spkTs(spkSelIx);
                [ chIx ]       = chIx(spkSelIx);
                [ cluLab ]       = cluLab(spkSelIx);
                
                [ chanLabSPK ] = lfpDat.chanLab(chIx);clear chIx;
                [ chanLabLFP ] = lfpDat.chanLab;
                
                %%
                spkSelIx = zeros(length(chanLabSPK),1);
                for curMW = 1:length(chanLabSPK)
                    spkSelIx(curMW) = find(strcmp(chanLabSPK(curMW),lfpDat.chanLab));
                end;
                
                tmp = [spkSelIx cluLab];
                
                plvSelIx = zeros(length(plvSigDat.chanLabSPK),1);
                for curMW = 1:length(plvSigDat.chanLabSPK)
                    plvSelIx(curMW) = find(strcmp(plvSigDat.chanLabSPK(curMW),lfpDat.chanLab));
                end;
                
                tmp2 = [plvSelIx plvSigDat.cluLab];
                
                [~,iC,iA] = unique(tmp2,'rows');
                tmp2 = tmp2(iC,:);
                
                [ plvIx ] = [];
                if ~isempty(tmp2)
                    for curMW = 1:size(tmp,1)
                        if any(sum(ones(size(tmp2,1),1)*tmp(curMW,:)==tmp2,2)==2)
                            plvIx  = [plvIx curMW];
                        end;
                    end;
                    
                    if isempty(plvIx)
                        error('totoNtata');
                    end;
                    
                    chck = zeros(length( plvIx ),1);
                    for curMW = 1:length( plvIx )
                        chck(curMW) = any(strcmp(chanLabSPK(plvIx(curMW)),plvSigDat.chanLabSPK));
                    end;
                    if sum(chck) ~= length( plvIx )
                        error('tataNtoto');
                    end;
                end;
                
                [ noPlvIx ] = setdiff(1:length(spkTs),plvIx);
                
                nPLVunits   = nPLVunits + length(plvIx);
                nNoPLVunits = nNoPLVunits + length(noPlvIx);
                
                %%
                if ~isempty(plvIx)
                    
                    %%
                    spkTsPLV   = spkTs(plvIx);
                    
                    %%
                    %                     for curMW = 1:length( spkTsPLV )
                    %                         x = spkTsPLV{curMW}(dt(2:end-1)>0 & dt(2:end-1) <=2e3,:);
                    %                         tmp = zeros(1001,size(x,2));
                    %                         for curTrl = 1:size(x,2)
                    %                             [tmp(:,curTrl)] = xcorr(x(:,curTrl),500);
                    %                         end;
                    %                         figure;
                    %                         bar(linspace(-500,500,size(tmp,1)),mean(tmp,2));
                    %                     end;
                    
                                        %%
                                        trlPool = sort([lfpDat.hitIdx;lfpDat.missIdx]);
                                        hitIdx = lfpDat.hitIdx;
                                        missIdx = lfpDat.missIdx;
                    
                                        clear lfpDat;
                    %
                    %%
                    np = length(spkTsPLV)*(length(spkTsPLV)-1)/2;
                    cnt = 0;
                    pIx = zeros(np,2);
                    for curMW = 1:length( spkTsPLV )-1
                        for curMW2 = curMW + 1:length( spkTsPLV )
                            cnt = cnt+1;
                            pIx(cnt,:) = [curMW curMW2];
                        end;
                    end;
                    
                    %%
                    for curP = 1:size(pIx,1)
                        poolCnt1 = poolCnt1+1;
                        x1 = spkTsPLV{pIx(curP,1)}(dt(2:end-1)>0 & dt(2:end-1) <=2e3,:);
                        x2 = spkTsPLV{pIx(curP,2)}(dt(2:end-1)>0 & dt(2:end-1) <=2e3,:);
                        tmp = zeros(101,size(x1,2));
                        for curTrl = 1:size(x1,2)
                            [tmp(:,curTrl)] = xcorr(x1(:,curTrl),x2(:,curTrl),50);
                        end;
                        coFI1(poolCnt1) = max(sum(tmp,2))/sum([x1(:);x2(:)]);
                        coFIH(poolCnt1) = max(sum(tmp(:,hitIdx),2))/sum([sum(x1(:,hitIdx),2);sum(x2(:,hitIdx),2)]);
                        coFIM(poolCnt1) = max(sum(tmp(:,missIdx),2))/sum([sum(x1(:,missIdx),2);sum(x2(:,missIdx),2)]);
                    end;
                    
                end;
                
                %%
                if length(noPlvIx) >1
                    
                    spkTsnoPLV = spkTs(noPlvIx);
                    
                    %%
                    np = length(spkTsnoPLV)*(length(spkTsnoPLV)-1)/2;
                    cnt = 0;
                    pIx = zeros(np,2);
                    for curMW = 1:length( spkTsnoPLV )-1
                        for curMW2 = curMW + 1:length( spkTsnoPLV )
                            cnt = cnt+1;
                            pIx(cnt,:) = [curMW curMW2];
                        end;
                    end;
                    
                    for curP = 1:size(pIx,1)
                        poolCnt2 = poolCnt2+1;
                        x1 = spkTsnoPLV{pIx(curP,1)}(dt(2:end-1)>2e3 & dt(2:end-1) <=3e3,:);
                        x2 = spkTsnoPLV{pIx(curP,2)}(dt(2:end-1)>2e3 & dt(2:end-1) <=3e3,:);
                        tmp = zeros(101,size(x1,2));
                        for curTrl = 1:size(x1,2)
                            [tmp(:,curTrl)] = xcorr(x1(:,curTrl),x2(:,curTrl),50);
                        end;
                        coFI2(poolCnt2) = max(sum(tmp,2))/sum([x1(:);x2(:)]);
                    end;
                    
                end;
                
                %%
            end;
            
        end;
    end;
end;

%%
figure;
subplot(121);
errorbar([1 2],[mean(coFI1) mean(coFI2)],[std(coFI1)./sqrt(length(coFI1)-1) std(coFI2)/sqrt(length(coFI2)-1)],'bs-');
set(gca,'XTick',[1 2]);
set(gca,'XTickLabel',{'PLV','no-PLV'});
ylabel('Co-Firing Index [a.u.]');

[hH,xH] = hist(coFI1,linspace( min([coFI1 coFI2]),max([coFI1 coFI2]),10 ));
[hM,xM] = hist(coFI2,linspace( min([coFI1 coFI2]),max([coFI1 coFI2]),10 ));
subplot(122);
hold on;
h = [];
h(1) = plot(xH,(hH./sum(hH)),'bo-');
h(2) = plot(xM,(hM./sum(hM)),'rx-');
legend(h,'PLV','no-PLV');
xlabel('Co-Figiring Index [a.u.]');
ylabel('Cumulative distribution [%]');

%%
figure;
subplot(121);
errorbar([1 2],[mean(coFIH) mean(coFIM)],[std(coFIH)./sqrt(length(coFIH)-1) std(coFIM)/sqrt(length(coFIM)-1)],'bs-');
set(gca,'XTick',[1 2]);
set(gca,'XTickLabel',{'Hits','Misses'});
ylabel('Co-Firing Index [a.u.]');

[hH,xH] = hist(coFIH,linspace( min([coFIH coFIM]),max([coFIH coFIM]),10 ));
[hM,xM] = hist(coFIM,linspace( min([coFIH coFIM]),max([coFIH coFIM]),10 ));
subplot(122);
hold on;
h = [];
h(1) = plot(xH,(hH./sum(hH)),'bo-');
h(2) = plot(xM,(hM./sum(hM)),'rx-');
legend(h,'Hits','Misses');
xlabel('Co-Figiring Index [a.u.]');
ylabel('Cumulative distribution [%]');