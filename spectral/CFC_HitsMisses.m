function [PAC]=CFC_HitsMisses(fls_HF, rdsPath, savePath);

cnt=0;
for curfile=1:size(fls_HF,1)
    tmpfl=fls_HF(curfile,1).name;
    spk2lfpHF=load([savePath, '/', tmpfl]);
    % extract patient, session, etc. info for file import
    ix = regexp(tmpfl,'_');
    pId = tmpfl(1:ix(1)-1);
    em = tmpfl(ix(1)+1:ix(2)-1);
    sesh = tmpfl(ix(2)+1:ix(4)-1);
    tw(1) = str2num(tmpfl(ix(9)+1));
    tw(2) = str2num(tmpfl(ix(9)+4));
    %% select locally coupled gamma channels
    chanCnt=0;
    HFchan={};
    HFreg={};
    for curpair=1:size(spk2lfpHF.chan_labels,1)
        tmpspkchan = spk2lfpHF.chan_labels{curpair,1};
        tmplfpchan = spk2lfpHF.chan_labels{curpair,3};
        loc_dist=strcmp(tmpspkchan(1:end-1),tmplfpchan(1:end-1));
        if loc_dist == 1
            chanCnt = chanCnt+1;
            plvall{1,1}(chanCnt,:) = spk2lfpHF.spk2LFPCouplingHM(curpair,:);
            plvH{1,1}(chanCnt,:) = spk2lfpHF.spk2LFPCouplingH(curpair,:);
            plvM{1,1}(chanCnt,:) = spk2lfpHF.spk2LFPCouplingM(curpair,:);
            HFchan{chanCnt,1}=tmplfpchan;
            HFreg{chanCnt,1}=tmplfpchan(1:end-1);
        end
    end
    
    %% select distally coupled theta channels
    tmpfl2=dir([savePath, pId, '_', em, '_', sesh, '_',...
    'spk2LFPCouplingHitsANDmisses_ppc_Sorting_2to40Hz_FDRfs_'...   
    num2str(tw(1)),'to', num2str(tw(2)), 'Secs_NotchOff.mat']);
    if ~isempty(tmpfl2)
        spk2lfpLF=load([savePath, tmpfl2.name]);
        chanCnt=0;
        LFchan={};
        LFreg={};
        for curpair=1:size(spk2lfpLF.chan_labels,1)
            tmpspkchan = spk2lfpLF.chan_labels{curpair,1};
            tmplfpchan = spk2lfpLF.chan_labels{curpair,3};
            loc_dist=strcmp(tmpspkchan(1:end-1),tmplfpchan(1:end-1));
            if loc_dist == 0
                chanCnt = chanCnt+1;
                plvall{1,2}(chanCnt,:) = spk2lfpLF.spk2LFPCouplingHM(curpair,:);
                plvH{1,2}(chanCnt,:) = spk2lfpLF.spk2LFPCouplingH(curpair,:);
                plvM{1,2}(chanCnt,:) = spk2lfpLF.spk2LFPCouplingM(curpair,:);
                LFchan{chanCnt,1}=tmplfpchan;
                LFreg{chanCnt,1}=tmplfpchan(1:end-1);
            end
        end
        HFchan=unique(HFchan);
        LFchan=unique(LFchan);
        HFreg=unique(HFreg);
        LFreg=unique(LFreg);
        %% Match Local Gamma Channel with Theta Channel
        lfidx=[];
        for curchan=1:length(HFchan)
            lfidx=[lfidx; strmatch(HFchan{curchan,1}(1:end-1),LFchan)];
        end
        if ~isempty(lfidx)
            lfidx=unique(lfidx);
            LFchan=LFchan(lfidx,1);
            [ p2d ] = [ rdsPath, pId, '/' ,em, '/',sesh, '/' ];
            [ lfpDatfN ] = dir( [ p2d, pId, '_',em, '_',sesh, '_lfpDataStimLockedSegmenteddownsampled.mat' ] );
            [ lfpDat ] = load( [ p2d, lfpDatfN.name ] ); % load the LFP data
            [ trlPool, hitIdx, missIdx, trlENC ] = organizeTrlIdxEM( lfpDat );
            [ lfp, ~, delIx, ~, ~, chanLabLFP ] = preprocLFP( lfpDat, trlPool, trlENC, [] );% preprocess LFP
            HFchidx=find(ismember(chanLabLFP,HFchan));
            LFchidx=find(ismember(chanLabLFP,LFchan));
            
            HFdelIx=delIx(1,HFchidx);
            LFdelIx=delIx(1,LFchidx);
            HFchanlabs=chanLabLFP(1,HFchidx);
            LFchanlabs=chanLabLFP(1,LFchidx);
            %% extract phase from LFP
            [ phi, phiTime, phiFreq ] = computeLFPphase( lfp(1,LFchidx), lfpDat.dsFs, [min(lfpDat.dsTrlTime) max(lfpDat.dsTrlTime)],'LF');
            LFrawlfp=lfp(1,LFchidx);% store raw lfp of phase providing channels to check for asymm waveforms
            [ pow, powTime, powFreq ] = computeLFPpower4CFC	( lfp(1,HFchidx), HFchanlabs,lfpDat.dsFs, [min(lfpDat.dsTrlTime) max(lfpDat.dsTrlTime)],'HF');
            phitoi(1)=nearest(phiTime,powTime(1));
            phitoi(2)=nearest(phiTime,powTime(end));
            %% Pair each HF channel with LF channel on a regional basis to extract local Gamma - Theta interaction
             clear lfp;
            lfwind=[5 11];
            lfwidx=find(spk2lfpLF.spk2LFPfreqAx >= lfwind(1) & spk2lfpLF.spk2LFPfreqAx <= lfwind(2));
            hfwind=[50 80];
            hfwidx=find(spk2lfpHF.spk2LFPfreqAx >= hfwind(1) & spk2lfpHF.spk2LFPfreqAx <= hfwind(2));
            %figure;
            for hfr=1:length(HFreg)
                cnt2=0;
                pacpow=[];
                pacphs=[];
                % Match up LF channels and HF channels within a region and
                % calculate PAC for each pairwise combination
                tmphfidx=strmatch(HFreg{hfr,1},HFchanlabs);
                tmplfidx=strmatch(HFreg{hfr,1},LFchanlabs);
                tmphfchans=HFchanlabs(1,tmphfidx);
                tmplfchans=LFchanlabs(1,tmplfidx);
                for lfi=1:length(tmplfchans)
                    s2lfidx=strmatch(tmplfchans{1,lfi},spk2lfpLF.chan_labels(:,3));
                    phsidx=strmatch(tmplfchans{1,lfi},LFchanlabs);
                    rawlfp=LFrawlfp{1,phsidx}(phitoi(1):phitoi(2),:)';
                    lfppcH=mean(spk2lfpLF.spk2LFPCouplingH(s2lfidx,:),1);
                    lfppcM=mean(spk2lfpLF.spk2LFPCouplingM(s2lfidx,:),1);
                    % Extract peak theta frequency for hits and misses
                    fidxH=extract_peakf(lfppcH,lfwidx,7);
                    fidxM=extract_peakf(lfppcM,lfwidx,7);
                    lffreqH=spk2lfpLF.spk2LFPfreqAx(fidxH);
                    lffreqM=spk2lfpLF.spk2LFPfreqAx(fidxM);
                    tmpdelIx{1,1}=LFdelIx(1,phsidx);
                    % Continue on Monday here ...
                    % Idea: Extract phase of dominant theta frequency for
                    % hits and misses separately as the two have different
                    % frequencies; do it in such a way that we have one
                    % variable LFphase which stores the phase
                    % 
                    phstrls=trlPool;
                    phstrls(tmpdelIx{1,1}{1,1})=[];
                    for nt=1:numel(phstrls)
                        htms=ismember(phstrls(nt),hitIdx);
                        if htms == 1
                            phsfidx=fidxH;
                        else
                            phsfidx=fidxM;
                        end
                        LFphase(nt,:,:)=angle(permute(phi{1,phsidx}(nt,1,phsfidx-3:phsfidx+3,phitoi(1):phitoi(2)),[1 3 4 2]));
                    end
                    
                    for hfi=1:length(tmphfchans)
                        cnt2=cnt2+1;
                        LFphase2=LFphase;% make copy of LF phase because trial selection changes depending on deleted trials of HF power
                        rawlfp2=rawlfp;
                        s2hfidx=strmatch(tmphfchans{1,hfi},spk2lfpHF.chan_labels(:,3));
                        powidx=strmatch(tmphfchans{1,hfi},HFchanlabs);
                        hfppcH=mean(spk2lfpHF.spk2LFPCouplingH(s2hfidx,:),1);
                        hfppcM=mean(spk2lfpHF.spk2LFPCouplingM(s2hfidx,:),1);
                        % Extract peak gamma frequency for hits and misses
                        fidxH=extract_peakf(hfppcH,hfwidx,10);
                        fidxM=extract_peakf(hfppcM,hfwidx,10);
                        hfidxH=nearest(powFreq,spk2lfpHF.spk2LFPfreqAx(fidxH));
                        hfidxM=nearest(powFreq,spk2lfpHF.spk2LFPfreqAx(fidxM));
                        hffreqH=powFreq(hfidxH-16:hfidxH+16);
                        hffreqM=powFreq(hfidxM-16:hfidxM+16);
                        % Idea: Same thing as for theta ... 
                        tmpdelIx{1,2}=HFdelIx(1,powidx);
                        powtrls=trlPool;
                        powtrls(tmpdelIx{1,2}{1,1})=[];
                        for nt=1:numel(powtrls)
                            htms=ismember(powtrls(nt),hitIdx);
                            if htms == 1
                                powfidx=hfidxH;
                            else
                                powfidx=hfidxM;
                            end
                            HFpow(nt,:,:)=permute(pow{1,powidx}(nt,1,powfidx-16:powfidx+16,:),[1 3 4 2]);
                        end
                        %% next we need to match trials because different channels will have different artefacts rejected ...
                        trlPool1=trlPool;
                        trlPool1(tmpdelIx{1,1}{1,1})=[];
                        trlPool2=trlPool;
                        trlPool2(tmpdelIx{1,2}{1,1})=[];
                        trlPool3=intersect(trlPool1,trlPool2);% trial IDs that are shared between LF and HF channels
                        tmphidx=find(ismember(trlPool3,hitIdx));
                        tmpmidx=find(ismember(trlPool3,missIdx));
                        % now find all trials in LF and HF channels that
                        % are not shared between the two and delete these
                        tmpdelIx{1,3}=find(~ismember(trlPool1,trlPool3));% for LF
                        tmpdelIx{1,4}=find(~ismember(trlPool2,trlPool3));% for HF
                        trlPool1(tmpdelIx{1,3})=[];
                        trlPool2(tmpdelIx{1,4})=[];
                        LFphase2(tmpdelIx{1,3},:,:)=[];
                        LFphase3=squeeze(LFphase2(:,4,:));
                        rawlfp2(tmpdelIx{1,3},:)=[];
                        HFpow(tmpdelIx{1,4},:,:)=[];
                        chk1=find(trlPool1-trlPool3);% check that selected trials match
                        chk2=find(trlPool2-trlPool3);
                        if isempty(chk1) && isempty(chk2)
                            [powsrtd,phssrtd,lfpsrtd,pactime]=computePACHM(LFphase3,rawlfp2,HFpow,[lffreqH lffreqM],powTime,[2 3],tmphidx,501);
                            [powsrtdH,phssrtdH,lfpsrtdH,pactime]=computePAC(LFphase3(tmphidx,:),...
                                rawlfp2(tmphidx,:),HFpow(tmphidx,:,:),lffreqH,powTime,[2 3],501);
                            [powsrtdM,phssrtdM,lfpsrtdM,pactime]=computePAC(LFphase3(tmpmidx,:),...
                                rawlfp2(tmpmidx,:),HFpow(tmpmidx,:,:),lffreqM,powTime,[2 3],501);
                            edges=linspace(-pi,pi,18);
                            nshuff=200;
                            [miz(:,:,cnt2),mir,mish]=computePAC_MI(LFphase2,HFpow,edges,nshuff);
                            [mizH(:,:,cnt2),mir,mish]=computePAC_MI(LFphase2(tmphidx,:,:),HFpow(tmphidx,:,:),edges,nshuff);
                            [mizM(:,:,cnt2),mir,mish]=computePAC_MI(LFphase2(tmpmidx,:,:),HFpow(tmpmidx,:,:),edges,nshuff);
                            mizHM(:,:,cnt2)=(mizH(:,:,cnt2)+mizM(:,:,cnt2))./2;
                        else
                            display('Error!!! Trials for LF and HF channels dont match');
                        end
                        % check if gamma peaks or troughs align with theta
                        [cchk,~]=corrcoef(powsrtd(17,:),phssrtd);
                        flpchk(cnt2)=abs(cchk(2,1));
                        if sign(cchk(2,1))<0
                            powsrtd=powsrtd.*-1;
                            powsrtdH=powsrtdH.*-1;
                            powsrtdM=powsrtdM.*-1;
                        end
                        %phssrtdi=interp(phssrtd,srnew/srold)
                        pacpow(:,:,cnt2)=powsrtd;
                        pacpowH(:,:,cnt2)=powsrtdH;
                        pacpowM(:,:,cnt2)=powsrtdM;
                        pacphs(cnt2,:)=phssrtd;
                        paclfp(cnt2,:)=lfpsrtd;
                        hffreq=-16:1:16;
                        figure
                        subplot(10,1,1:7);
                        pcolor(pactime,hffreq,powsrtd);colorbar;
                        shading interp;
                        title([LFchanlabs{1,phsidx} '-to-' HFchanlabs{1,powidx}]);
                        subplot(10,1,8:10);
                        plot(pactime,phssrtd);colorbar;
                        xlim([pactime(1) pactime(end)]);
                        figure;
                        pcolor([-3:3],[-16:16],mir(:,:,1));
                        shading interp;
                        clear tmpphs* tmppow* HFpow
                    end
                    clear LFphase
                end
                if cnt2 > 0
                    cnt=cnt+1;
%                     maxidx=find(mizHM(17,4,:)==max(mizHM(17,4,:)));
%                     PAC.pow(:,:,cnt)=pacpow(:,:,maxidx);
%                     PAC.powH(:,:,cnt)=pacpowH(:,:,maxidx);
%                     PAC.powM(:,:,cnt)=pacpowM(:,:,maxidx);
%                     PAC.phs(cnt,:)=pacphs(maxidx,:);
%                     PAC.lfp(cnt,:)=paclfp(maxidx,:);
%                     PAC.mi(:,:,cnt)=miz(:,:,maxidx);
%                     PAC.miH(:,:,cnt)=mizH(:,:,maxidx);
%                     PAC.miM(:,:,cnt)=mizM(:,:,maxidx);
                    if cnt==1
                        PAC.pow=pacpow;
                        PAC.powH=pacpowH;
                        PAC.powM=pacpowM;
                        PAC.phs=pacphs;
                        PAC.lfp=paclfp;
                        PAC.mi=miz;
                        PAC.miH=mizH;
                        PAC.miM=mizM;
                    else
                        PAC.pow=cat(3,PAC.pow, pacpow);
                        PAC.powH=cat(3,PAC.powH, pacpowH);
                        PAC.powM=cat(3,PAC.powM, pacpowM);
                        PAC.phs=cat(1,PAC.phs, pacphs);
                        PAC.lfp=cat(1,PAC.lfp, paclfp);
                        PAC.mi=cat(3,PAC.mi,miz);
                        PAC.miH=cat(3,PAC.miH,mizH);
                        PAC.miM=cat(3,PAC.miM,mizM);
                    end
                 
                    clear flpchk pac* miz mizH mizM mizHM 
                end
            end
            
    
        end
    end
end
% 
% 
% 
% PAC.time=-250:1:250;
% PAC.freq=-16:16;
% cmap=cbrewer('div','RdBu',20);
% figure;
% subplot(10,1,1:7);
% pcolor(PAC.time,PAC.freq,mean(PAC.pow,3));colorbar;
% colormap(cmap);
% shading interp;
% ax=gca;
% ax.XTickLabel={};
% xlim([-250 250]);
% title('CFC all');
% subplot(10,1,8:10);
% plot(PAC.time,mean(PAC.lfp,1));colorbar;
% xlim([-250 250]);
% ax=gca;
% ax.XTick=[-200 -100 0 100 200];
% ax.XTickLabel={'-2';'-1';'0';'1';'2' };
% 
% figure;
% subplot(10,1,1:7);
% pcolor(PAC.time,PAC.freq,mean(PAC.powH,3));colorbar;
% colormap(cmap);
% shading interp;
% ax=gca;
% ax.XTickLabel={};
% xlim([-250 250]);
% title('CFC Hits');
% subplot(10,1,8:10);
% plot(PAC.time,mean(PAC.lfp,1));colorbar;
% xlim([PAC.time(1) PAC.time(end)]);
% xlim([-250 250]);
% ax=gca;
% ax.XTick=[-200 -100 0 100 200];
% ax.XTickLabel={'-2';'-1';'0';'1';'2' };
% 
% figure;
% subplot(10,1,1:7);
% pcolor(PAC.time,PAC.freq,mean(PAC.powM,3));colorbar;
% colormap(cmap);
% shading interp;
% xlim([-250 250]);
% ax=gca;
% ax.XTickLabel={};
% title('CFC Misses');
% subplot(10,1,8:10);
% plot(PAC.time,mean(PAC.lfp,1));colorbar;
% xlim([-250 250]);
% ax=gca;
% ax.XTick=[-250 -100 0 100 250];
% ax.XTickLabel={'-2';'-1';'0';'1';'2' };
% 
% figure;
% subplot(2,3,1);
% pcolor([-3:3],[-16:16],mean(PAC.mi,3));
% shading interp;
% subplot(2,3,2);
% pcolor([-3:3],[-16:16],mean(PAC.miH,3));
% shading interp;
% subplot(2,3,3);
% pcolor([-3:3],[-16:16],mean(PAC.miM,3));
% shading interp;
% 
% subplot(2,3,4);
% pcolor([-3:3],[-16:16],mean(PAC.miH,3)-mean(PAC.miM,3));
% shading interp;

%% Subfunctions
function [fidx]=extract_peakf(ppc,widx,def)


[pks,loc]=findpeaks(ppc);
tmploc=find(loc>=widx(1) & loc<=widx(end));
if ~isempty(tmploc)
    fidx=loc(tmploc(find(pks(tmploc)==max(pks(tmploc)))));
else
    fidx=def;
end

end

end
