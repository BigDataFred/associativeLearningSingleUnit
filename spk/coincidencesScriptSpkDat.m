function coincidencesScriptSpkDat( pId, expMode, spkMode, spk2LFPmode, savePath, rdsPath  )

%%
if nargin ==0
    pId = {'P02'};
    expMode = {'fVSpEM'};
    spkMode = {'Sorting'};
    spk2LFPmode = {'plv'};
    savePath = '/home/rouxf/resultsSpikeFieldJun2019/';
    rdsPath = '/media/rouxf/rds-share/iEEG_DATA/MICRO/';
end;

%%
% if isempty( gcp('Nocreate') )
%     parpool(36);
% end;

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
        for curMode = 1:length(spkMode)
            
            %%
            [ tmp ] = dir([rdsPath,pId{curPat},'/',expMode{curExp},'/']);
            
            %%
            if ~isempty(tmp)
                
                [sesh] = extractSeshLabels(tmp);clear tmp;
                
                %%
                for curSesh = 1:length(sesh)
                    
                    %%
                    [ p2d ] = [rdsPath,pId{curPat},'/',expMode{curExp},'/',sesh{curSesh},'/'];
                    
                    %%
                    fN = dir([p2d,pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_lfpDataStimLockedSegmenteddownsampled.mat']);
                    [ lfpDat ] = load([p2d,fN.name])
                    
%                     %%
%                     fN = dir([p2d,pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spkDataStimLockedSegmented.mat']);
%                     [ spkDat ] = load([p2d,fN.name])
                    
                    for curSPK2LFP = 1:length(spk2LFPmode)
                        %%
                        [ fNspk2LFP ] = [pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spk2LFPCoupling_',spk2LFPmode{curSPK2LFP},'_',spkMode{curMode},'_alpha0.05_nRand200_lowFreq.mat'];
                        [ fNspk2LFP ] = dir([savePath,fNspk2LFP]);
                        [ plvSigDat ] = load([savePath,fNspk2LFP.name]);
                        
                        %%
                        [ spkParams ] = load([savePath,pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spkParams_',spkMode{curMode},'.mat']);
                        
                        %% find units that are active during encoding
                        [spkSelIx] = filterUnits4SPKrate(spkParams.spkTs,spkParams.frM, [0. 4.], 50, lfpDat.dsTrlTime);
                        [ chIx ] = 1:length(spkParams.spkTs);
                        
                        %%
                        [ spkTs ]      = spkParams.spkTs( spkSelIx );
                        [ chIx ]       = chIx( spkSelIx );                        
                        [ chanLabSPK ] = lfpDat.chanLab( chIx );clear chIx;
                        
                        %%
                        [ spkTs ]      = spkTs( plvSigDat.sigIxSPK );
                        [ chanLabSPK ] = chanLabSPK( plvSigDat.sigIxSPK );
                        
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
                            hitIdx = lfpDat.hitIdx;
                            missIdx = lfpDat.missIdx;
                            
                            clear lfpDat;
                            
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
                                [ x1 ] = spkTsPLV{pIx(curP,1)}(dt(2:length(dt)-1)>0 & dt(2:length(dt)-1) <=2e3,:);
                                [ x2 ] = spkTsPLV{pIx(curP,2)}(dt(2:length(dt)-1)>0 & dt(2:length(dt)-1) <=2e3,:);
                                [ tmp ] = zeros(101,size(x1,2));
                                for curTrl = 1:size(x1,2)
                                    [ tmp(:,curTrl) ] = xcorr(x1(:,curTrl),x2(:,curTrl),50);
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