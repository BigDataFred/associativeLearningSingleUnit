%%
addpath('~/prj/Bham/code/mcode/utils/');
addpath('~/tbx/CircStat2012a/');
addpath(genpath('/home/rouxf/tbx/eeglab14_1_1b/functions/sigprocfunc/'));
addpath(genpath('~/tbx/chronux_2_11/spectral_analysis/'));
addpath('~/tbx/fieldtrip-20170618/');
ft_defaults;

%%
pId = {'P07'};%{'P02','P04','P05','P07','P08','P09','P03ERL','P22AMS','P23AMS'};%
expMode = {'fVSpEM','cnEM'};
savePath = '~/resultsSpikeFieldOct18II/';

%% start parallel ppol
if isempty(gcp('nocreate'))
    parpool(36,'SpmdEnabled', false); 
end;

%%
movingwinLF = [.5 .01];

paramsLF                  = [];
paramsLF.pad              = 4;
paramsLF.fpass            = [0 30];
paramsLF.tapers           = [1 1];
paramsLF.trialave         = 0;

movingwinHF = [.25 .01];

paramsHF                  = [];
paramsHF.pad              = 5;
paramsHF.fpass            = [30 170];
paramsHF.tapers           = [3 5];
paramsHF.trialave         = 0;

cfgtf                     = [];
cfgtf.method              = 'wavelet';
cfgtf.output              = 'pow';
cfgtf.pad                 = 'nextpow2';
cfgtf.foi                 = 0:100;
cfgtf.width               = 4;
cfgtf.output              = 'fourier';

%%
for curPat = 1:length(pId) % loop over patients
    for curExp = 1:length(expMode) % loop over exp-modes (either fVSpEM or cnEM)
        
        [ tmp ] = dir(['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pId{curPat},'/',expMode{curExp},'/']);% check for existing data
        
        if ~isempty(tmp) % get the timestamp labels for each session
            
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
            sesh(cnt+1:end) = []
            clear tmp;
            
            %%
            for curSesh = 1:length(sesh) % loop over the recording-sesssions
                
                [ p2d ] = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pId{curPat},'/',expMode{curExp},'/',sesh{curSesh},'/']; %path for preproc data
                [ fN ] = dir([p2d,pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_lfpDataStimLockedSegmenteddownsampled.mat']);% filename of preproc data
                
                [ lfpDat ] = load([p2d,fN.name])% load data
                      
                %%
                [ trlPool ] = sort([lfpDat.hitIdx;lfpDat.missIdx]);
                [ hitIdx ]  = lfpDat.hitIdx;
                [ missIdx ] = lfpDat.missIdx;
                [ trlENC ]  = 1:length(lfpDat.trlENC);
                
                if ~all(trlENC' == trlPool) || (length(trlPool) == length(trlENC))~=1 || ~all(hitIdx == trlPool(ismember(trlPool,hitIdx))) || ~all(sort(missIdx) == trlPool(ismember(trlPool,missIdx)))
                    error('error wrong trial assignment');
                end;
                
                %%
                paramsLF.Fs               = lfpDat.dsFs;
                paramsHF.Fs               = lfpDat.dsFs;
                cfgtf.toi                 = lfpDat.dsTrlTime(1):1/lfpDat.dsFs:lfpDat.dsTrlTime(length(lfpDat.dsTrlTime));
                
                %% extract the labels of Behnke-Fried electrodes
                BFlab = cell(length(lfpDat.chanLab),1);
                for curChan = 1:length( lfpDat.chanLab )% loop over microwires
                    tmp = lfpDat.chanLab{curChan};
                    ix = regexp(tmp,'\d{1}');
                    BFlab(curChan) = {tmp(1:ix-1)};% only keep the BF label, ignore MW index
                end;
                BFid = unique(BFlab);% extract unique labels
                
                BFix = cell(length(BFid),1);
                for curBF = 1:length(BFid)                   
                    BFix{curBF} = find(strcmp(BFid{curBF},BFlab));% map micro-wires to BF-label (ie extract indexes)
                end;                                               
                
                %% calculate root mean square energy for micro-wires
                for curBF = 1:length(BFid)
                    rms = zeros(length(BFix{curBF}),1);
                    for curMW = 1:length(BFix{curBF})
                        rms(curMW) = mean(sqrt(mean(lfpDat.LFPseg{BFix{curBF}(curMW)}(lfpDat.dsTrlTime>=-0.5 & lfpDat.dsTrlTime <=5,:).^2,1)));
                    end;
                    rms = (rms-mean(rms))./std(rms);
                    %delete Micro-wires with RMS above threshold
                    lfpDat.LFPseg(BFix{curBF}(abs(rms)>3)) = []; 
                    lfpDat.chanLab(BFix{curBF}(abs(rms)>3)) = [];
                end;                
                clear BF*;
                
                %%
                [ lfp ]   = cell(1,length( lfpDat.LFPseg));% lfp data      
                [ erp ]   = cell(1,length( lfpDat.LFPseg));% erp data  
                [ delIx ] = cell(1,length( lfpDat.LFPseg));% indexes with delete lfp trials (important for spk data!)
                [ n ]     = zeros(1,length( lfpDat.LFPseg));
                
                %% pre-allocate and initialize
                parfor curChan = 1:length( lfpDat.LFPseg )

                    [ tmpLFP ] = lfpDat.LFPseg{curChan}(:,ismember(lfpDat.trlSel,lfpDat.trlENC)); % extract trials corresponding to encoding                  
                    
                    % compute RMS of LFP over time
                    rms = sqrt(mean(tmpLFP(lfpDat.dsTrlTime>=-0.5 & lfpDat.dsTrlTime <=5,:).^2,1));
                    rms = ( rms-mean(rms) )./std(rms);% standardize RMS
                    
                    % compute z-score of LFP
                    m  = ones(size(tmpLFP,1),1)*mean(tmpLFP,1);
                    sd = ones(size(tmpLFP,1),1)*std(tmpLFP,0,1);
                    z = ( tmpLFP-m )./ sd;
                    z = max(abs(z(lfpDat.dsTrlTime>=-0.5&lfpDat.dsTrlTime <=5,:)),[],1);
                    
                    delIx{curChan} = unique([find(abs(rms) >= 4) find(z >= 4)]);% keep indexes of trials where LFP-activity is out of range (outliers)
                    
                    tmpLFP(:,delIx{curChan}) = []; %delete trials where LFP-activity is out of range (outliers)
                    
                    % LFP-data must have at least 25 trials in total (hits + misses) per session to perform spectral analysis
                    tmp = trlPool;
                    tmp(delIx{curChan}) = [];
                    n(curChan) = length(find(ismember(tmp,hitIdx)));
                    if (~isempty( tmpLFP) ) && (n(curChan)>25) 
                        
                        m  = ones(size(tmpLFP,1),1)*mean(tmpLFP,1); % mean of LFP
                        sd = ones(size(tmpLFP,1),1)*std(tmpLFP,0,1); % SD of LFP
                        
                        lfp{curChan} = (tmpLFP-m);%./sd;% mean-centered single trial lfp data
                        erp{curChan} = lfp{curChan}; % evoked poential
                    end;
                end;
                clear LFPseg;
                
                %%
                selIx = [];
                for curChan = 1:length( lfp )
                    if ~isempty(lfp{curChan})
                        selIx = [selIx;curChan];
                    end;
                end;
                chanLabLFP = lfpDat.chanLab(selIx);
                
                lfp = lfp(selIx);
                erp = erp(selIx);
                delIx = delIx(selIx);
                n     = n(selIx);
                
                %%
                if ( ~isempty(lfp) && ~isempty(erp) && ~isempty(selIx) )
                    
                    %%
                    [ dsTrlTime ] = lfpDat.dsTrlTime;
                    [ erpH ] = cell(1,length(erp));
                    [ erpM ] = cell(1,length(erp));
                    for curChan = 1:length( erp )
                        tmp = trlPool;
                        tmp(delIx{curChan}) = [];    
                        if length(find(ismember(tmp,hitIdx))) == n(curChan)
                            erpH{curChan} = mean(erp{curChan}(dsTrlTime>=-.5 & dsTrlTime<=5,ismember(tmp,hitIdx)),2);
                            erpM{curChan} = mean(erp{curChan}(dsTrlTime>=-.5 & dsTrlTime<=5,ismember(tmp,missIdx)),2);
                        end;
                    end;
                    [ dsTrlTime ] = dsTrlTime(dsTrlTime>=-.5 & dsTrlTime<=5);
                    
                    saveName =[pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_ERPTrialSeg_HitsANDmisses.mat'];
                    save([savePath,saveName],'erpH','erpM','dsTrlTime','chanLabLFP');
                    clear tmp erp* dsTrlTime;
                    
                    %% pre-allocate and initialize
                    [ SxxLF ]    = cell(1,length( selIx )); % low-freq pow of sngle trials
                    [ SxxeLF ]   = cell(1,length( selIx ));% low-freq pow of ERP
                    [ SxxHF ]    = cell(1,length( selIx ));% high-freq pow of sngle trials
                    [ SxxeHF ]   = cell(1,length( selIx ));% high-freq pow of ERP    
                    [ itc ]     = cell(1,length( selIx ));% inter-trial coh
                    [ itcH ]     = cell(1,length( selIx ));% inter-trial coh
                    [ itcM ]     = cell(1,length( selIx ));% inter-trial coh
                    %[ phi ]      = cell(1,length( selIx ));% single trial phase
                    
                    %%
                    parfor curChan = 1:length( lfp )
                        fprintf([num2str(curChan),'/',num2str(length( lfp ))]);
                        if ~isempty( lfp{curChan} )
                            % low-frequency power (multi-taper)
                            [SxxLF{curChan},~,~] = mtspecgramc( gradient(lfp{curChan}')', movingwinLF, paramsLF);
                            [SxxeLF{curChan},~,~] = mtspecgramc( gradient(mean(lfp{curChan},2)')', movingwinLF, paramsLF);
                        end;
                        fprintf('\n');
                    end;
                    [~,tx,fx] = mtspecgramc(ones(size(lfp{1},1),1),movingwinLF, paramsLF);
                    
                    %%
%                     errClf = cell(1,length(SxxLF));
%                     errCrandlf = cell(1,length(SxxLF));
%                     for curChan = 1:length(SxxLF)
%                         fprintf([num2str(curChan),'/',num2str(length(SxxLF))]);
%                         fprintf('\n');
%                         
%                         A = squeeze(mean(SxxLF{curChan}(tx-7>=0.25 & tx-7<=1.8,1:4:end,:),1));
%                         B = squeeze(mean(SxxLF{curChan}(tx-7>=2.25 & tx-7<=3.8,1:4:end,:),1));
%                         
%                         [errClf{curChan},errCrandlf{curChan}] = decodingAvsB([A';B']',[zeros(size(A,2),1);ones(size(B,2),1)],100);                        
%                     end;
                    
                    %% split single trial power data into hits and misses and save data to disc separately
                    % low freq hits
                    SxxLFH = cell(1,length(SxxLF));
                    for curChan = 1:length(SxxLF)
                        trlENC2 = 1:length(trlENC);
                        trlENC2(delIx{curChan}) = [];
                        if (~isempty(SxxLF{curChan})) && (any(ismember(trlENC2,hitIdx)))
                            SxxLFH{curChan} = squeeze(mean(SxxLF{curChan}(:,:,ismember(trlENC2,hitIdx)),3));
                        end;
                    end;
                    
                    saveName =[pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spectrogramTrialSeg_lowFreqHits.mat'];
                    save([savePath,saveName],'SxxLFH','tx','fx','chanLabLFP','-v7.3');
                    clear SxxLFH;
                    
                    % low freq misses
                    SxxLFM = cell(1,length(SxxLF));
                    for curChan = 1:length(SxxLF)
                        trlENC2 = 1:length(trlENC);
                        trlENC2(delIx{curChan}) = [];
                        if (~isempty(SxxLF{curChan})) && (any(ismember(trlENC2,missIdx)))
                            SxxLFM{curChan} = squeeze(mean(SxxLF{curChan}(:,:,ismember(trlENC2,missIdx)),3));
                        end;
                    end;
                    
                    saveName =[pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spectrogramTrialSeg_lowFreqMisses.mat'];
                    save([savePath,saveName],'SxxLFM','tx','fx','chanLabLFP','-v7.3');
                    clear SxxLFM;
                    
                    %% average single trial power over trials
                    for curChan = 1:length( SxxLF )
                        SxxLF{curChan} = squeeze(mean( SxxLF{curChan},3));
                    end;
                    
                    %% save average power
                    saveName =[pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spectrogramTrialSeg_lowFreq.mat'];
                    save([savePath,saveName],'SxxLF','tx','fx','chanLabLFP','-v7.3');
                    clear SxxLF;
                    
                    saveName =[pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spectrogramTrialSegERP_lowFreq.mat'];
                    save([savePath,saveName],'SxxeLF','tx','fx','chanLabLFP');
                    clear SxxeLF;
                    
                    %%
                    parfor curChan = 1:length( lfp )
                        fprintf([num2str(curChan),'/',num2str(length( lfp ))]);
                        if ~isempty(lfp{curChan})
                            % LFP-data must have at least 25 trials in total (hits + misses) per session to perform spectral analysis
                            if (~isempty(lfp{curChan})) && (size(lfp{curChan},2)>25)
                                % high-frequency power (multi-taper)
                                [SxxHF{curChan},~,~] = mtspecgramc( gradient(lfp{curChan}')', movingwinHF, paramsHF);
                                [SxxeHF{curChan},~,~] = mtspecgramc( gradient(mean(lfp{curChan},2)')', movingwinHF, paramsHF);
                            end;
                        end;
                        fprintf('\n');
                    end;
                    [~,tx,fx] = mtspecgramc(ones(size(lfp{1},1),1),movingwinHF, paramsHF);
                    
                    %%
%                     errChf = cell(1,length(SxxHF));
%                     errCrandhf = cell(1,length(SxxHF));
%                     for curChan = 1:length(SxxHF)
%                         fprintf([num2str(curChan),'/',num2str(length(SxxHF))])
%                         A = squeeze(mean(SxxHF{curChan}(tx-7>=0 & tx-7<2,1:18:end,:),1));
%                         B = squeeze(mean(SxxHF{curChan}(tx-7>=2 & tx-7<4,1:18:end,:),1));
%                         
%                         [errChf{curChan},errCrandhf{curChan}] = decodingAvsB([A';B']',[zeros(size(A,2),1);ones(size(B,2),1)],500);
%                         fprintf('\n');
%                     end;
                    
                    %% split single trial power data into hits and misses and save data to disc separately
                    % high freq hits
                    SxxHFH = cell(1,length(SxxHF));
                    for curChan = 1:length(SxxHF)
                        trlENC2 = 1:length(trlENC);
                        trlENC2(delIx{curChan}) = [];
                        if (~isempty(SxxHF{curChan})) && (any(ismember(trlENC2,hitIdx)))
                            SxxHFH{curChan} = squeeze(mean(SxxHF{curChan}(:,:,ismember(trlENC2,hitIdx)),3));
                        end;
                    end;
                    
                    saveName =[pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spectrogramTrialSeg_highFreqHits.mat'];
                    save([savePath,saveName],'SxxHFH','tx','fx','chanLabLFP','-v7.3');
                    clear SxxHFH;
                    
                    % high freq misses
                    SxxHFM = cell(1,length(SxxHF));
                    for curChan = 1:length(SxxHF)
                        trlENC2 = 1:length(trlENC);
                        trlENC2(delIx{curChan}) = [];
                        if (~isempty(SxxHF{curChan})) && (any(ismember(trlENC2,missIdx)))
                            SxxHFM{curChan} = squeeze(mean(SxxHF{curChan}(:,:,ismember(trlENC2,missIdx)),3));
                        end;
                    end;
                    
                    saveName =[pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spectrogramTrialSeg_highFreqMisses.mat'];
                    save([savePath,saveName],'SxxHFM','tx','fx','chanLabLFP','-v7.3');
                    clear SxxHFM;
                    
                    %% average single trial power over trials
                    for curChan = 1:length( SxxHF )
                        SxxHF{curChan} = squeeze(mean( SxxHF{curChan},3));
                    end;
                    
                    %% save average power
                    saveName =[pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spectrogramTrialSeg_highFreq.mat'];
                    save([savePath,saveName],'SxxHF','tx','fx','chanLabLFP','-v7.3');
                    clear SxxHF;
                    
                    saveName =[pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spectrogramTrialSegERP_highFreq.mat'];
                    save([savePath,saveName],'SxxeHF','tx','fx','chanLabLFP');
                    clear SxxeHF;
                    
                    %%
                    [ dsTrlTime ] = lfpDat.dsTrlTime;
                    [ Fs ]        = lfpDat.dsFs;
                    [ pval1 ] = zeros(length( lfp ),1);
                    [ pval2 ] = zeros(length( lfp ),1);                    
                    parfor curChan = 1:length( lfp )
                        if ~isempty(lfp{curChan})
                            % dummy structure for fieldtrip
                            dum                     = [];
                            dum.fsample             = Fs;
                            dum.label               = {'dumChan1'};
                            dum.trial               = cell(1,size(lfp{curChan},2));
                            dum.time                = cell(1,size(lfp{curChan},2));
                            for curTrl = 1:size(lfp{curChan},2)
                                dum.trial{curTrl}            = gradient( lfp{curChan}(:,curTrl)' );
                                dum.time{curTrl}             = dsTrlTime;
                            end;
                            
                            % single trial & time resolved spectral analysis (wavelet)
                            %[pow{curChan}] = ft_freqanalysis( cfgtf, dum );
                            
                            % single trial & time resolved spectral analysis (wavelet)
                            %[phi{curChan}] = ft_freqanalysis( cfgtf, dum );% compute phase
                            [phi] = ft_freqanalysis( cfgtf, dum );% compute phase
                            
                            % compute lfp inter-trial phase coherence (ITC)
                            tmp = squeeze(phi.fourierspctrm./abs(phi.fourierspctrm));% keep complex values
                            
                            trlENC2 = 1:length( trlENC );
                            trlENC2( delIx{curChan} ) = [];
                            
                            hIx = find(ismember(trlENC2,hitIdx));
                            mIx = find(ismember(trlENC2,missIdx));
                            tIx = find(phi.time>=0 & phi.time <=4);
                            fIx = find(phi.freq>=0 & phi.freq <=30);
                            
                            dItcR = zeros(200,1);%size(tmp,2),size(tmp,3)
                            for curPer = 1:200
                                fprintf([num2str(curPer),'/',num2str(200)]);
                                rIx = randperm(size(tmp,1));
                                sel1 = rIx(1:length(hIx));
                                rIx(sel1) = [];
                                sel2 = rIx;
                                if (length(sel1) ~= length(hIx)) || ( length(sel2) ~= length(mIx) )
                                    error('warning index selection erroneous');
                                end;
                                rItcH = squeeze(1/size(tmp(sel1,fIx,tIx),1)*abs(nansum(tmp(sel1,fIx,tIx),1)));
                                rItcM = squeeze(1/size(tmp(sel2,fIx,tIx),1)*abs(nansum(tmp(sel2,fIx,tIx),1)));
                                d = rItcH-rItcM;
                                dItcR(curPer) = max(abs(d(:)));
                                fprintf('\n');
                            end;
                            
                            nT = length(phi.time);
                            shf = [250 750 1500];
                            ItcR = zeros(200,1);%size(tmp,2),size(tmp,3)
                            for curPer = 1:200
                                fprintf([num2str(curPer),'/',num2str(200)]);
                                shDat = zeros(size(tmp));
                                for curTrl = 1:size(tmp,1)
                                    fc = randperm(2);
                                    fc2 = randperm(length(shf));
                                    if fc(1) == 1
                                        shIx = [nT-shf(fc2(1)):nT 1:nT-shf(fc2(1))-1];%shift forward
                                    else
                                        shIx = [shf(fc2(1)):nT 1:shf(fc2(1))-1];%shift backward
                                    end;
                                    
                                    shDat(curTrl,:,:) = tmp(curTrl,:,shIx);
                                end;
                                rItc = squeeze(1/size(shDat(:,fIx,tIx),1)*abs(nansum(shDat(:,fIx,tIx),1)));
                                ItcR(curPer) = max(rItc(:));
                                fprintf('\n');
                            end;
                            
                            itch = squeeze(1/size(tmp(hIx,fIx,tIx),1)*abs(nansum(tmp(hIx,fIx,tIx),1)));
                            itcm = squeeze(1/size(tmp(mIx,fIx,tIx),1)*abs(nansum(tmp(mIx,fIx,tIx),1)));
                            dItcE = itch-itcm;
                            dItcE = max(abs(dItcE(:)));                            
                            pval1(curChan) = length(find(dItcR>=dItcE))/length(dItcR);
                            
                            ItcE = squeeze(1/size(tmp(:,fIx,tIx),1)*abs(nansum(tmp(:,fIx,tIx),1)));
                            ItcE = max(ItcE(:));
                            pval2(curChan) = length(find(ItcR>=ItcE))/length(ItcR);
                            
                            [ itc{curChan} ] = squeeze(1/size(tmp,1)*abs(nansum(tmp,1)));
                            [ itcH{curChan} ] = squeeze(1/size(tmp(hIx,:,:),1)*abs(nansum(tmp(hIx,:,:),1)));
                            [ itcM{curChan} ] = squeeze(1/size(tmp(mIx,:,:),1)*abs(nansum(tmp(mIx,:,:),1)));
                        end;
                    end;
                    
                    %% save phaseDat & ITC
                    % get time and freq information for ITC
                    dum                     = [];
                    dum.fsample             = lfpDat.dsFs;
                    dum.label               = {'dumChan1'};
                    dum.trial               = cell(1,1);
                    dum.time                = cell(1,1);
                    dum.trial{1}            = zeros(1,length(dsTrlTime));
                    dum.time{1}             = dsTrlTime;

                    [phi] = ft_freqanalysis( cfgtf, dum );% compute phase
                    if ~isempty(phi)
                        itcTime = phi.time;
                        itcFreq = phi.freq;
                    end;
                    clear lfp;
                    
                    saveName =[pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spectrogramTrialSegITC_lowFreq.mat'];
                    save([savePath,saveName],'itc','itcH','itcM','itcTime','itcFreq','chanLabLFP','pval1','pval2');
                    clear itc phi;
                    
%                     tmp = cell(1,length(phi));
%                     for curChan = 1:length( phi)
%                         tmp{curChan} = phi{curChan}.fourierspctrm;
%                     end;
%                     phi = tmp;clear tmp;
%                     phiTime = itcTime;
%                     phiFreq = itcFreq;
%                     
%                     saveName =[pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_TrialSegPhase_',num2str(round(min(phiFreq))),'-',num2str(round(max(phiFreq))),'.mat'];
%                     save([savePath,saveName],'phi','phiTime','phiFreq','chanLabLFP','-v7.3');
%                     clear phi*;
                    
                end;
            end;
        end;
    end;
end;

