function timeFreqAnalysisScriptEMtask( pId, expMode, savePath, rdsPath )
%%
addpath(genpath('~/tbx/chronux_2_11/spectral_analysis/'));
addpath('~/tbx/fieldtrip-20170618/');
addpath('~/tbx/eeglab14_1_1b/functions/sigprocfunc/');
ft_defaults;

%% start parallel ppol
if isempty(gcp('nocreate'))
    parpool(36,'SpmdEnabled', false); 
end;

%%
movingwinLF = [1 .01];

paramsLF                  = [];
paramsLF.pad              = 4;
paramsLF.fpass            = [2 30];
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
        
        [ tmp ] = dir([rdsPath,pId{curPat},'/',expMode{curExp},'/']);% check for existing data
        
        if ~isempty(tmp) % get the timestamp labels for each session
            
            [sesh] = extractSeshLabels(tmp);
            clear tmp;
            
            %%
            for curSesh = 1:length(sesh) % loop over the recording-sesssions
                
                [ p2d ] = [rdsPath,pId{curPat},'/',expMode{curExp},'/',sesh{curSesh},'/']; %path for preproc data
                [ fN ] = dir([p2d,pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_lfpDataStimLockedSegmenteddownsampled.mat']);% filename of preproc data
                
                [ lfpDat ] = load([p2d,fN.name])% load data                      
                
                %%
                paramsLF.Fs               = lfpDat.dsFs;
                paramsHF.Fs               = lfpDat.dsFs;
                cfgtf.toi                 = lfpDat.dsTrlTime(1):1/lfpDat.dsFs:lfpDat.dsTrlTime(length(lfpDat.dsTrlTime));
                
                %%
                [trlPool,hitIdx,missIdx,trlENC] = organizeTrlIdxEM(lfpDat);
                
                [lfp,erp,delIx,selIx,n,chanLabLFP] = preprocLFP( lfpDat, trlPool, trlENC );
                
                %%
                if ( ~isempty(lfp) && ~isempty(erp) && ~isempty(selIx) )
                    
                    %%
                    [ dsTrlTime ] = lfpDat.dsTrlTime;
                    [ erpH ] = cell(1,length(erp));
                    [ erpM ] = cell(1,length(erp));
                    for curChan = 1:length( erp )
                        tmp = trlPool;
                        tmp(delIx{curChan}) = [];    
                        if length(tmp) == n(curChan)
                            x = eegfilt( erp{ curChan }', lfpDat.dsFs, 3, 12, 0, [], [], 'fir1', [] )';                            
                            erpH{curChan} = mean(x(dsTrlTime>=-.5 & dsTrlTime<=5,ismember(tmp,hitIdx)),2);
                            erpM{curChan} = mean(x(dsTrlTime>=-.5 & dsTrlTime<=5,ismember(tmp,missIdx)),2);
                            erp{ curChan } = mean(x,2);
                        end;
                    end;
                    [ dsTrlTime ] = dsTrlTime(dsTrlTime>=-.5 & dsTrlTime<=5);
                    
                    saveName =[pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_ERPTrialSeg_HitsANDmisses.mat'];
                    save([savePath,saveName],'erp','erpH','erpM','dsTrlTime','chanLabLFP');
                    clear tmp erp* dsTrlTime;
                    
                    %% pre-allocate and initialize
                    [ SxxLF ]    = cell(1,length( selIx )); % low-freq pow of sngle trials
                    [ SxxeLF ]   = cell(1,length( selIx ));% low-freq pow of ERP
                    [ SxxHF ]    = cell(1,length( selIx ));% high-freq pow of sngle trials
                    [ SxxeHF ]   = cell(1,length( selIx ));% high-freq pow of ERP    
                    [ itc ]     = cell(1,length( selIx ));% inter-trial coh
                    [ itcHvsM ]     = cell(1,length( selIx ));% inter-trial coh
                    
%                     [ itcH ]     = cell(1,length( selIx ));% inter-trial coh
%                     [ itcM ]     = cell(1,length( selIx ));% inter-trial coh
                    
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
                    
                    %% split single trial power data into hits and misses and save data to disc separately                    
                    SxxLFH = cell(1,length(SxxLF));
                    SxxLFM = cell(1,length(SxxLF));
                    for curChan = 1:length(SxxLF)
                        trlENC2 = 1:length(trlENC);
                        trlENC2(delIx{curChan}) = [];
                        if (~isempty(SxxLF{curChan})) && (any(ismember(trlENC2,hitIdx)))
                            [hIx] = find(ismember(trlENC2,hitIdx));
                            [mIx] = find(ismember(trlENC2,missIdx));
                            [hIx,mIx] = stratisfyCondIndexes( hIx, mIx );
                            x = [];
                            for curSamp = 1:size( ix1,1)
                                [ x(curSamp,:) ] = squeeze(mean(SxxLF{curChan}(:,:,hIx(curSamp,:)),3));
                            end;                                                                                       
                            SxxLFH{curChan} = mean(x,1);                            
                             x = [];
                            for curSamp = 1:size( ix1,1)
                                [ x(curSamp,:) ] = squeeze(mean(SxxLF{curChan}(:,:,mIx(curSamp,:)),3));
                            end;                             
                            SxxLFM{curChan} = mean(x,1);
                        end;
                    end;
                    
                    % low freq hits
                    saveName =[pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spectrogramTrialSeg_lowFreqHits.mat'];
                    save([savePath,saveName],'SxxLFH','tx','fx','chanLabLFP','-v7.3');
                    clear SxxLFH;
                    
                    % low freq misses                    
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
                    
                    %% split single trial power data into hits and misses and save data to disc separately
                    SxxHFH = cell(1,length(SxxHF));
                    SxxHFM = cell(1,length(SxxHF));
                    for curChan = 1:length(SxxHF)
                        trlENC2 = 1:length(trlENC);
                        trlENC2(delIx{curChan}) = [];
                        if (~isempty(SxxHF{curChan})) && (any(ismember(trlENC2,hitIdx)))
                            [hIx] = find(ismember(trlENC2,hitIdx));
                            [mIx] = find(ismember(trlENC2,missIdx));
                            [hIx,mIx] = stratisfyCondIndexes( hIx, mIx );
                            x = [];
                            for curSamp = 1:size( ix1,1)
                                [ x(curSamp,:) ] = squeeze(mean(SxxHF{curChan}(:,:,hIx(curSamp,:)),3));
                            end;                                                                                       
                            SxxHFH{curChan} = mean(x,1);                            
                             x = [];
                            for curSamp = 1:size( ix1,1)
                                [ x(curSamp,:) ] = squeeze(mean(SxxHF{curChan}(:,:,mIx(curSamp,:)),3));
                            end;                             
                            SxxHFM{curChan} = mean(x,1);
                        end;
                    end;
                    
                    % high freq hits                    
                    saveName =[pId{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spectrogramTrialSeg_highFreqHits.mat'];
                    save([savePath,saveName],'SxxHFH','tx','fx','chanLabLFP','-v7.3');
                    clear SxxHFH;
                    
                    % high freq misses                    
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
%                     [ pval1 ] = zeros(length( lfp ),1);
%                     [ pval2 ] = zeros(length( lfp ),1);                    
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
                            
%                             dItcR = zeros(200,1);%size(tmp,2),size(tmp,3)
%                             for curPer = 1:200
%                                 fprintf([num2str(curPer),'/',num2str(200)]);
%                                 rIx = randperm(size(tmp,1));
%                                 sel1 = rIx(1:length(hIx));
%                                 rIx(sel1) = [];
%                                 sel2 = rIx;
%                                 if (length(sel1) ~= length(hIx)) || ( length(sel2) ~= length(mIx) )
%                                     error('warning index selection erroneous');
%                                 end;
%                                 rItcH = squeeze(1/size(tmp(sel1,fIx,tIx),1)*abs(nansum(tmp(sel1,fIx,tIx),1)));
%                                 rItcM = squeeze(1/size(tmp(sel2,fIx,tIx),1)*abs(nansum(tmp(sel2,fIx,tIx),1)));
%                                 d = rItcH-rItcM;
%                                 dItcR(curPer) = max(abs(d(:)));
%                                 fprintf('\n');
%                             end;
                            
%                             nT = length(phi.time);
%                             shf = [250 500 750 1000 1250 1500 1750 2500 2750 3000];
%                             ItcR = zeros(200,1);%size(tmp,2),size(tmp,3)
%                             for curPer = 1:200
%                                 fprintf([num2str(curPer),'/',num2str(200)]);
%                                 shDat = zeros(size(tmp));
%                                 for curTrl = 1:size(tmp,1)
%                                     fc = randperm(2);
%                                     fc2 = randperm(length(shf));
%                                     if fc(1) == 1
%                                         shIx = [nT-shf(fc2(1)):nT 1:nT-shf(fc2(1))-1];%shift forward
%                                     else
%                                         shIx = [shf(fc2(1)):nT 1:shf(fc2(1))-1];%shift backward
%                                     end;
%                                     
%                                     shDat(curTrl,:,:) = tmp(curTrl,:,shIx);
%                                 end;
%                                 rItc = squeeze(1/size(shDat(:,fIx,tIx),1)*abs(nansum(shDat(:,fIx,tIx),1)));
%                                 ItcR(curPer) = max(rItc(:));
%                                 fprintf('\n');
%                             end;
                            
                            [hIx,mIx] = stratisfyCondIndexes( hIx, mIx );
                                                                                    
                            x = [];
                            for curSamp = 1:size( hIx,1)
                                [ x(curSamp,:) ] = squeeze(1/size(tmp(hIx(curSamp,:),fIx,tIx),1)*abs(nansum(tmp(hIx(curSamp,:),fIx,tIx),1)));
                            end;
                            [ itch ] = mean(x,1);
                            
                            x = [];
                            for curSamp = 1:size( mIx,1)
                                [ x(curSamp,:) ] = squeeze(1/size(tmp(mIx(curSamp,:),fIx,tIx),1)*abs(nansum(tmp(mIx(curSamp,:),fIx,tIx),1)));
                            end;
                            [ itcm ] = mean(x,1);
                            
%                             dItcE = itch-itcm;
                            
%                             dItcE = max(abs(dItcE(:)));                            
%                             pval1(curChan) = length(find(dItcR>=dItcE))/length(dItcR);
                            
%                             ItcE = squeeze(1/size(tmp(:,fIx,tIx),1)*abs(nansum(tmp(:,fIx,tIx),1)));
%                             ItcE = max(ItcE(:));
%                             pval2(curChan) = length(find(ItcR>=ItcE))/length(ItcR);
                            
                            [ itc{curChan} ] = squeeze(1/size(tmp,1)*abs(nansum(tmp,1)));
                            [ itcHvsM{curChan} ] = itch-itcm;
%                             [ itcH{curChan} ] = squeeze(1/size(tmp(hIx,:,:),1)*abs(nansum(tmp(hIx,:,:),1)));
%                             [ itcM{curChan} ] = squeeze(1/size(tmp(mIx,:,:),1)*abs(nansum(tmp(mIx,:,:),1)));
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
                    save([savePath,saveName],'itc','itcHvsM','itcTime','itcFreq','chanLabLFP');%'itcH','itcM','pval1','pval2'
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

