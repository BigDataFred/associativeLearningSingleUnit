function analyze_EM_data2(pID,expMode)

%% set data labels
if nargin ==0
    pID = 'P04ERL';
    expMode = 'fVSpEM';
end;

%% set Matlab path environment
addpath('/home/rouxf/tbx/releaseDec2015/');
addpath(genpath('/home/rouxf/tbx/osort-v3-rel/'));
addpath(genpath('/home/rouxf/tbx/wave_clus/'));
addpath(genpath('/home/rouxf/tbx/chronux_2_11/'));
addpath('/home/rouxf/tbx/fieldtrip-20170618/');
ft_defaults;

addpath('/home/rouxf/prj/Bham/code/mcode/params/');
addpath('/home/rouxf/prj/Bham/code/mcode/helper/');
addpath('/home/rouxf/prj/Bham/code/mcode/helper/logfile_readers/');
addpath('/home/rouxf/prj/Bham/code/mcode/utils/spikes/');
addpath('/home/rouxf/prj/Bham/code/mcode/utils/LFP/');
addpath('/home/rouxf/prj/Bham/code/mcode/utils/TTL/');
addpath('/home/rouxf/prj/Bham/code/mcode/project_EM/');

%% recruit workers
% Create a parallel pool if none exists
mode = 'debug';
%mode = '';
if strcmp(mode,'')
    if isempty(gcp('nocreate'))
        parpool(36,'SpmdEnabled',false);
    end
else
    if isempty(gcp('nocreate'))
        parpool(1,'SpmdEnabled',false);
    end
end;

%% set data parameters
Fs =32e3;
timeStampsPerSec = (1/Fs)*1e6;

[par] = set_parameters_Bham(Fs);

%% parameters to read the TTL data
FieldSelection1 = [];
FieldSelection1(1) = 1;%timestamps
FieldSelection1(2) = 0;% EventIDs
FieldSelection1(3) = 1;%TTLs
FieldSelection1(4) = 0;% Extras
FieldSelection1(5) = 0;% Event strings

%% search for existing data
[p2LogDat] = ['/media/rouxf/rds-share/Archive/MICRO/',pID,'/log/',expMode,'/'];
%[p2LogDat] = ['/media/rouxf/rds-share/Archive/MICRO/',pID,'/log/',expMode,'/_EMtask_logs/'];

if strcmp(mode,'')
    [rpath]    = ['/home/rouxf/in/',expMode,'/'];
    [savepath] = ['/home/rouxf/out/',expMode,'/'];
else
    [rpath] = ['/media/rouxf/rds-share/Archive/MICRO/',pID,'/',expMode,'/'];
    %[rpath]    = ['/media/rouxf/rds-share/Archive/MICRO/',pID,'/P04ERL_EMtask_MacroMicro/'];
    [savepath] = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID,'/',expMode,'/'];
end;
rpath
savepath

chck = dir(savepath);
if isempty(chck)
    mkdir(savepath);
end;

[sesh] = getSeshLabelsEM(rpath)

[origSesh] = getSeshLabelsEM(['/media/rouxf/rds-share/Archive/MICRO/',pID,'/',expMode,'/']);
%[origSesh] = getSeshLabelsEM(['/media/rouxf/rds-share/Archive/MICRO/',pID,'/',pID,'_EMtask_MacroMicro/']);

%% loop over recording sessions
for seshSel = 1:length(sesh)
    
    [p2d] = [rpath,sesh{seshSel},'/']
    [CSCfiles] = dir([p2d,'Micro_*.ncs']);
    if isempty(CSCfiles)
        [CSCfiles] = dir([p2d,'*.ncs']);
    end;
    
    %% read the Logfile
    [logDat]= getLogFileDataEMtask(pID,p2LogDat,origSesh,find(strcmp(origSesh,sesh{seshSel})));
    
    if ( ~isempty(logDat.LogDat1.log) && ~isempty(logDat.LogDat2.log) )
        
        %% save the Logfile
        logFileName = [ pID ,'_', sesh{seshSel}, '_LogFile_EMtask_LogDat.mat' ];
        chck = dir([ savepath,sesh{seshSel}]);
        if isempty(chck)
            mkdir( [ savepath,sesh{seshSel}] );
        end;
        save( [ savepath,sesh{seshSel},'/',logFileName ] , 'logDat' );
        
        %% load the logfile data
        [blockSz] = diff(logDat.LogDat1.idx,[],2)+1;
        
        ix = 1:blockSz(1);
        dum = {};
        for it = 1:length(blockSz)
            dum{it} = [ix ix+length(ix)];
            if it < length(blockSz)
                ix = [dum{it}(end)+1:dum{it}(end)+blockSz(it+1)];
            end;
        end;
        
        trlENC = {};
        trlRET = {};
        for it = 1:length(dum)
            trlENC{it} = dum{it}(1:blockSz(it));
            trlRET{it} = dum{it}(blockSz(it)+1:2*blockSz(it));
        end;
        [trlENC] = [trlENC{:}];
        [trlRET] = [trlRET{:}];
        
        if ( length(trlENC) ~= length(logDat.RTs) ) || ( length(trlRET) ~= length(logDat.RTs) )
            error('trial assignment must be equal');
        end;
        if ( length(trlENC) ~= length(trlRET) )
            error('trial assignment must be equal');
        end;
        if (length(trlENC) ~= size(logDat.LogDat1.log,1) ) || (length(trlRET) ~= size(logDat.LogDat2.log,1) )
            error('trial assignment must be equal');
        end;
        
        %% read the TTL event data
        [EVfile] = dir([p2d,'*.nev']);
        
        [TimeStampsTTL, ttls,~] = Nlx2MatEV_v3( [p2d,EVfile.name], FieldSelection1, 1, 1, [] );
        
        [events] = zeros(size(ttls,2),2);
        events(:,1) = TimeStampsTTL';
        events(:,2) = ttls';
        firstTstamp = events(1,1);
        
        events(:,1) = (events(:,1)-firstTstamp)./1e6;% convert from micro-s to seconds
        
        TTLix = (TimeStampsTTL-firstTstamp)/timeStampsPerSec+1;
        TTLix = double(uint64(TTLix));
        
        [selIx7] = find(ttls ==7);
        
        if ( length(selIx7) ~= sum([length(trlENC) length(trlRET)]) )
            error('Trigger and trial indexes must have equal length');
        end;
        
        [selIx7E] = selIx7(trlENC);
        [selIx7R] = selIx7(trlRET);
        
        [selIx255] = find(ttls==255);
        
        eventsE = events(selIx7E,:);
        eventsR = events(selIx7R,:);
        
        [selIx0E] = findTTLoffsetIdx(ttls,selIx7E,[1 2]);
        [selIx0R] = findTTLoffsetIdx(ttls,selIx7R,[1 2]);
        [selIx0II] = findTTLoffsetIdx(ttls,selIx255,[3 2]);
        
        %%
        sortedSpikes = cell(1,length(CSCfiles));
        wltCoeffs = cell(1,length(CSCfiles));
        chanLab = cell(1,length(CSCfiles));
        LFPsig = cell(1,length(CSCfiles));
        LFPsig2 = cell(1,length(CSCfiles));
        
        %%
        tic;
        parfor it = 1:length(CSCfiles)
            %for it = 1%17%[6:7]%2%1:length(CSCfiles)
            
            fprintf([num2str(it),'/',num2str(length( CSCfiles) )]);
            
            %%
            dum = CSCfiles(it).name;
            dum(regexp(dum,'_')) = [];
            dum(regexp(dum,'CSC'):regexp(dum,'CSC')+2) = [];
            dum(regexp(dum,'.ncs'):end) = [];
            chanLab{it} = dum;
            
            %% parameters to read the CSC data
            FieldSelection2 = [];
            FieldSelection2(1) = 1;%timestamps
            FieldSelection2(2) = 0;
            FieldSelection2(3) = 0;%sample freq
            FieldSelection2(4) = 0;
            FieldSelection2(5) = 1;%samples
            
            %%
            [~, dataSamples,hdr] = Nlx2MatCSC_v3([p2d,CSCfiles(it).name], FieldSelection2, 1, 1, []);
            
            chck = regexp(hdr,'ADBitVolts');
            selIdx = [];
            for jt = 1:length(chck);selIdx(jt) = ~isempty(chck{jt});end;
            selIdx = find(selIdx~=0);
            scalef = str2double(hdr{selIdx}(min(regexp(hdr{selIdx},'\d')):end));
            
            chck = regexp(hdr,'SamplingFrequency');
            selIdx = [];
            for jt = 1:length(chck);selIdx(jt) = ~isempty(chck{jt});end;
            selIdx = find(selIdx~=0);
            Fs = str2double(hdr{selIdx}(min(regexp(hdr{selIdx},'\d')):end));
            
            %% flatten
            [dataSamples] = double(dataSamples(:))';
            [dataSamples] = dataSamples.*scalef.*1e6;
            
            %% remove TTL artefact
            ntrl = [1 2 4 8];
            tw = [1 2 3 4 5];
            
            if strcmp(mode,'')
                if ~isempty(selIx7E) && ~isempty(selIx0E)
                    [ix] = findParametersForTTLRemoval(ntrl,tw,dataSamples,TTLix,selIx7E,selIx0E,Fs);
                    [dataSamples] = cleanARTIFACTfromLFP(dataSamples,TTLix(selIx7E),[32*0.5 32*tw(ix(2))],ntrl(ix(1)),2,Fs);
                    [dataSamples] = cleanARTIFACTfromLFP(dataSamples,TTLix(selIx0E),[32*0.5 32*2*tw(ix(4))],ntrl(ix(3)),2,Fs);
                end;
                if ~isempty(selIx7R) && ~isempty(selIx0R)
                    [ix] = findParametersForTTLRemoval(ntrl,tw,dataSamples,TTLix,selIx7R,selIx0R,Fs);
                    [dataSamples] = cleanARTIFACTfromLFP(dataSamples,TTLix(selIx7R),[32*0.5 32*tw(ix(2))],ntrl(ix(1)),2,Fs);
                    [dataSamples] = cleanARTIFACTfromLFP(dataSamples,TTLix(selIx0R),[32*0.5 32*2*tw(ix(4))],ntrl(ix(3)),2,Fs);
                end;
                if ~isempty(selIx255) && ~isempty(selIx0II)
                    [ix] = findParametersForTTLRemoval(ntrl,tw,dataSamples,TTLix,selIx255,selIx0II,Fs);
                    [dataSamples] = cleanARTIFACTfromLFP(dataSamples,TTLix(selIx255),[32*0.5 32*tw(ix(2))],ntrl(ix(1)),2,Fs);
                    [dataSamples] = cleanARTIFACTfromLFP(dataSamples,TTLix(selIx0II),[32*0.5 32*tw(ix(4))],ntrl(ix(3)),2,Fs);
                end;
            end;
            %         % use interpolation method in case cleaning unsucessful
            %         %[dataSamples] = interpLFP(dataSamples,TTLix(selIx7),[0.002 0.006],Fs,'linear');
            %         %[dataSamples] = interpLFP(dataSamples,TTLix(selIx0),[0.002 0.006],Fs,'linear');
            
            %% store the LFP data for later
            LFPsig{it} = dataSamples;
            
            %% spike detection (1st pass)
            par2 = par;
            par2.fnamespc = [par2.fnamespc,num2str(it)];
            par2.fname_in = [par2.fname_in,num2str(it)];
            par2.segments = 6;%floor(length(dataSamples)/(5*Fs));
            par2.interpolation = 'y';
            
            nsmp = floor(length(dataSamples)/(par2.segments));%Fs;
            ix = 1:nsmp;
            spikeTimestamps = [];
            for kt = 1:par2.segments
                
                [~,~, ~, spikeTimestampsTmp,~,~,~] = amp_detect(dataSamples(ix),par2);
                
                spikeTimestamps = [spikeTimestamps ix(spikeTimestampsTmp)];
                
                if kt < par2.segments-1
                    ix = ix+nsmp;
                else
                    ix = ix(end)+1:length(dataSamples);
                end;
            end;
            
            %% interpolate LFP around spike times
            if ~isempty( spikeTimestamps )
                %use linear interpolation to remove spike shape from LFP
                [dataSamples] = interpLFP(dataSamples,spikeTimestamps,Fs,'linear',[0.002 0.006]);
                
                % use template subtration to remove spike shape from LFP
                %[dataSamples] = cleanARTIFACTfromLFP(dataSamples,spikeTimestamps,[32*0.5 32*2.5],3,2,Fs);
            end;
            
            %%
            LFPsig2{it} = dataSamples; % store the spike interpolated LFP signal for later use
            fprintf('\n');
            
        end;
        
        %     %% extract indexes of Microwires for each BF shank
        %     chck = regexp(chanLab,'\d{1,8}');
        %     chck  = [chck{:}];
        %
        %     BF = cell(1,length(chanLab));
        %     parfor it = 1:length(chanLab)
        %         BF(it) = {chanLab{it}(1:chck(it)-1)};
        %     end;
        %     BFid = unique(BF);
        %
        %     %% spike detection and sorting
        %     [b_detect,a_detect] = ellip(par.detect_order,0.1,40,[par.detect_fmin par.detect_fmax]*2/par.sr);
        %     logTxt = cell(1,length( CSCfiles));
        %     LFPsig3 = cell(1,length( CSCfiles));
        %     for it = 1:length( BFid )
        %
        %         fprintf([num2str(it),'/',num2str(length( CSCfiles) )]);
        %
        %         %% find the best reference signal to bring down noise in the spike band ie which maximizes kurtosis
        %         sel = find(strcmp(BF,BFid(it)));
        %
        %         refKurt=zeros(length(sel),1);
        %         for kt = 1:length(sel)
        %             refKurt(kt) = kurtosis( filtfilt( b_detect, a_detect, LFPsig{sel(kt)} ) );
        %         end;
        %         trsh = mean( refKurt )+2*std(refKurt);
        %
        %         D=zeros(length(sel),1);
        %         for kt = 1:length(sel)
        %             sel2 = setdiff(sel,kt);
        %             d=zeros(length(sel2),1);
        %             for lt = 1:length(sel2)
        %                 d(lt) = kurtosis( filtfilt( b_detect, a_detect, LFPsig{sel2(lt)}-LFPsig2{sel(kt)} ) );
        %             end;
        %             D(kt) = mean(d);% average kurtosis in the spike band for given MW-reference
        %         end;
        %
        %         if any(D>trsh)
        %             [~,mIx] = max(D);
        %             for jt = 1:length(sel)
        %                 [LFPsig3{sel(jt)}] = LFPsig{sel(jt)} -LFPsig2{sel(mIx)};
        %             end;
        %             logTxt(it) = {['re-referencing channel ',chanLab{it},' to channel ',chanLab{sel(mIx)}]};
        %         else
        %             logTxt(it) = {['']};
        %             for jt = 1:length(sel)
        %                 [LFPsig3{sel(jt)}] = LFPsig{sel(jt)};
        %             end;
        %         end;
        %         fprintf('\n');
        %
        %     end;
        %     LFPsig2 = LFPsig3;
        %     clear LFPsig3 d refKurt sel;
        %
        %     %%
        %     fid = fopen([savepath,sesh{seshSel},'/',pID,'_',expMode,'_',sesh{seshSel},'_SpkrefList.txt'],'w');
        %     for it = 1:length(logTxt)
        %         if ~isempty(logTxt{it})
        %             fprintf(fid,'%s\n',logTxt{it});
        %         else
        %             fprintf(fid,'%s\n','');
        %         end;
        %     end;
        %     clear fid logTxt;
        
        [ LFPsig2 ] = LFPsig;
        
        %%
        parfor it = 1:length( CSCfiles )
            %for it = 1
            fprintf([num2str(it),'/',num2str(length( CSCfiles) )]);
            [dataSamples] = LFPsig2{it};
            
            %% spike detection (2nd pass)
            par2 = par;
            par2.fnamespc = [par2.fnamespc,num2str(it)];
            par2.fname_in = [par2.fname_in,num2str(it)];
            par2.segments = 6;%floor(length(dataSamples)/(5*Fs));
            par2.interpolation = 'n';
            
            nsmp = floor(length(dataSamples)/(par2.segments));%Fs;
            ix = 1:nsmp;
            noiseTraces = [];
            spikeWaveforms = [];
            spikeTimestamps = [];
            noiseSDTmp = [];
            
            for kt = 1:par2.segments
                
                [~,spikeWaveformsTmp, ~, spikeTimestampsTmp,~,~,noiseTracesTmp] = amp_detect(dataSamples(ix),par2);
                
                spikeWaveforms  = [spikeWaveforms;spikeWaveformsTmp];
                spikeTimestamps = [spikeTimestamps ix(spikeTimestampsTmp)];
                noiseSDTmp      = [noiseSDTmp mean(std(noiseTracesTmp,0,2),1)];
                
                if size(noiseTracesTmp,1)>1
                    noiseTracesTmp2=[];
                    noiseTracesTmp2(1:size(noiseTracesTmp,1) ,1) = ones( size(noiseTracesTmp,1), 1 )*1i;
                    noiseTracesTmp2(1:size(noiseTracesTmp,1),2:size(noiseTracesTmp,2)+1)=noiseTracesTmp;
                    noiseTraces = [noiseTraces; noiseTracesTmp2];
                    noiseTracesTmp2 = [];
                    noiseTracesTmp = [];
                end;
                if kt < par2.segments-1
                    ix = ix+nsmp;
                else
                    ix = ix(end)+1:length(dataSamples);
                end;
            end;
            
            %% sanity check
            % there should now be no more TTL pulses in the timestamp data
            if any(ismember(TTLix,spikeTimestamps))
                delIx = find( ismember(spikeTimestamps,TTLix) );
                spikeWaveforms(delIx,:) = [];
                spikeTimestamps(delIx) = [];% if there are any left throw them away
                if any(ismember(TTLix,spikeTimestamps))
                    error('Spiketimes may be contaminated by TTL pulses');
                end;
            end;
            
            %% decorrelate and upsample spike-waveforms
            if ( ~isempty(spikeWaveforms) ) && ( size(spikeWaveforms,1)>1 ) && ( size(spikeWaveforms,2) >1 )
                
                [trans] = posthocWhiten(noiseTraces, spikeWaveforms,[]);
                [trans] = [trans(:,1)*ones(1,2) trans trans(:,end)*ones(1,2)];
                
                par2.interpolation = 'y';
                [trans] = int_spikes(trans,par);
                
                [spikeWaveforms] = [spikeWaveforms(:,1)*ones(1,2) spikeWaveforms spikeWaveforms(:,end)*ones(1,2)];
                [spikeWaveforms] = int_spikes(spikeWaveforms,par);
                
            end;
            
            %%
            waveclus                            = [];
            waveclus.spikes                     = spikeWaveforms;
            waveclus.index                      = spikeTimestamps;
            
            [dim] = size(spikeWaveforms);
            
            %% do spike sorting
            par2.filename = [CSCfiles(it).name];
            
            dum                     = [];
            dum.newSpikeTimes       = [];
            dum.assignedCluster     = [];
            dum.wavf                = [];
            dum.num_clus            = [];
            
            dum2                    = [];
            
            %%
            if dim(1) > dim(2)% the number of spike events must larger than the number of samples in the waveform
                [dum,dum2] = doSpikeSorting_waveclus( waveclus , par2 );
                
                dum.wavf_decorr         = trans;%spikeWaveforms2;%trans;%
                dum.SD                  = mean(noiseSDTmp(noiseSDTmp>0));
                trans = [];
                
                %             %% force cluster-membership
                %             [selIx1] = find(dum.assignedCluster ==0);
                %             [selIx2] = find(dum.assignedCluster ~=0);
                %
                %             [class_out] = force_membership_wc(dum2(selIx2,:), dum.assignedCluster(selIx2), dum2(selIx1,:), par);
                %             dum.assignedCluster(selIx1(class_out~=0)) = class_out(class_out~=0);
                
                %             %% delete noise clusters
                %             [delIx] = find(dum.assignedCluster==0);
                %
                %             dum.assignedCluster( delIx ) = [];
                %             dum.newSpikeTimes( delIx ) = [];
                %             dum.wavf(delIx,:) = [];
                %             dum.wavf_decorr(delIx,:) = [];
                %             dum2(delIx,:) = [];
                
            end;
            
            %% save the clustered data
            sortedSpikes{it} = dum;
            wltCoeffs{it} = dum2;
            fprintf('\n');
            
        end;
        clear LFPsig2;
        
        %%
        parfor it = 1:length( CSCfiles)
            %for it = 1
            
            fprintf([num2str(it),'/',num2str(length( CSCfiles) ),'\n']);
            
            %%
            if strcmp(mode,'')
                
                [ dataSamples ] = LFPsig{it};
                
                %% interpolate LFP around spike times
                if ~isempty( sortedSpikes{it}.assignedCluster )
                    %use linear interpolation to remove spike shape from LFP
                    [dataSamples] = interpLFP(dataSamples,sortedSpikes{it}.newSpikeTimes,Fs,'linear',[0.002 0.006]);
                    
                    % use template subtration to remove spike shape from LFP
                    %             cID = unique(sortedSpikes{it}.assignedCluster);
                    %             for kt = 1:length(cID)
                    %                 [dataSamples] = cleanARTIFACTfromLFP(dataSamples,sortedSpikes{it}.newSpikeTimes(sortedSpikes{it}.assignedCluster==cID(kt)),[32*0.5 32*2.5],3,2,Fs);
                    %             end;
                end;
                
                %% clean LFP of line noise
                s = 15;
                [step] = round(length(dataSamples)/(s*Fs));
                
                [dataSamples] = cleanLFPfromLineNoise(dataSamples,Fs,step,s);
                
                %% notch filter LFP
                %         [b,a] = butter(4,[49 51]./(Fs/2),'stop');% apply band-stop for LFP
                %         [dataSamples] = filtfilt(b,a,dataSamples);
                
                %         %% lowpass filter LFP
                %         bpFreq = [300];
                %         [b] = fir1(3*floor(Fs/bpFreq(1)),bpFreq./(Fs/2),'low');% apply band-pass for LFP
                %         [dataSamples] = filtfilt(b,1,dataSamples);
                
                %% save LFP data
                LFPsig{it} = dataSamples;
                dataSamples = [];
            end;
        end;
        
        %%
        toc;
        fprintf('Finished raw data import & preprocessing, moving on\n');
        
        %% time vector for LFP
        [lfpTime] = [1:size(LFPsig{1},2)]./Fs;
        
        %%    extract reaction times from logfile
        [ encRT ] = str2double(logDat.LogDat1.log(:,end));
        [ retRT ] = str2double(logDat.LogDat2.log(:,7));
        [ RTs ]   = [ encRT;retRT ];
        
        %% indexes of trials that were later remebered and forgotten
        [hitIdx]  = logDat.ix{4};%HITS
        [missIdx] = [logDat.ix{5};logDat.ix{6}];%MISSES
        
        %%
        [pre] = 7;% pres stimulus range in s
        [post] = 7;% post-stimulus range in s
        [trlTime] = -pre:1/Fs:post;
        
        [trlSel] = [trlENC trlRET]';
        [ttlIX] = [TTLix(selIx7E) TTLix(selIx7R)]';
        [events] = [eventsE; eventsR];
        
        %% compute the trl matrix based on event time stamps
        %[trl] = computeTRLmatrix(pre,post,eventsE{jt},lfpTime);% sanity
        %check
        
        %% compute the trl matrix based on sample indexes
        [trl] = zeros(length(ttlIX),2);
        for it = 1:length(ttlIX)
            [trl(it,:)] = [ttlIX(it)-pre*Fs ttlIX(it)+post*Fs];
        end;
        
        if any(diff(trl,[],2)./Fs ~= abs(pre)+post)
            error('sample indices are out of range');
        end;
        
        %% segment the LFP data
        [nsmp]    = length(trl(1,1):trl(1,2));
        if length(trlTime) ~= nsmp
            error('inconsistent number of samples');
        end;
        
        [LFPseg] = segmentLFPdata(LFPsig,trl,nsmp);
        
        %% segment the spike data
        [sortedSpikesSEG] = segmentSpikeData(sortedSpikes,pre,post,events,lfpTime);
        
        %%
        savename = [pID,'_',expMode,'_',sesh{seshSel},'_spkDataStimLockedSegmented.mat'];
        save([savepath,sesh{seshSel},'/',savename],'sortedSpikesSEG','wltCoeffs','Fs','chanLab','RTs','trlSel','trlENC','trlRET','hitIdx','missIdx','-v7.3');
        clear sortedSpikesSEG;
        
        %% use fieldtrip to downsample LFP data
        [LFPseg,dsFs] = downsampleLFP(LFPseg,trlTime,1e3);
        [dsTrlTime] = -pre:1/dsFs:post;
        if (size(LFPseg{1},1) ~= length(dsTrlTime)) || (size(LFPseg{1},2)~=length(trlSel))
            error('trial and sample numbers out of range');
        end;
        
        %%
        savename = [pID,'_',expMode,'_',sesh{seshSel},'_lfpDataStimLockedSegmenteddownsampled.mat'];
        save([savepath,sesh{seshSel},'/',savename],'LFPseg','dsTrlTime','dsFs','chanLab','RTs','trlSel','trlENC','trlRET','hitIdx','missIdx','-v7.3');
        
        %%
        ttlIX(1:length(trlENC)) = round(ttlIX(1:length(trlENC))+3*Fs+encRT*Fs);
        ttlIX(length(trlENC)+1:length(trlENC)+length(trlRET)) = round(ttlIX(length(trlENC)+1:length(trlENC)+length(trlRET))+2*Fs+retRT*Fs);
        [trl] = zeros(length(ttlIX),2);
        for it = 1:length(ttlIX)
            [trl(it,:)] = [ttlIX(it)-pre*Fs ttlIX(it)+post*Fs];
        end;
        
        if any(diff(trl,[],2)./Fs ~= abs(pre)+post)
            error('sample indices are out of range');
        end;
        
        %% segment the LFP data
        [nsmp]    = length(trl(1,1):trl(1,2));
        if length(trlTime) ~= nsmp
            error('inconsistent number of samples');
        end;
        
        [LFPseg] = segmentLFPdata(LFPsig,trl,nsmp);
        
        %% segment the spike data
        fprintf(num2str(size( events(length(trlENC)+1:length(trlENC)+length(trlRET)))));
        fprintf('\n');
        fprintf(num2str(size(retRT)));
        fprintf('\n');
        
        events(1:length(trlENC),1) = events(1:length(trlENC),1)+3+encRT;
        events(length(trlENC)+1:length(trlENC)+length(trlRET),1) = events(length(trlENC)+1:length(trlENC)+length(trlRET),1)+2+retRT;
        
        [sortedSpikesSEG] = segmentSpikeData(sortedSpikes,pre,post,events,lfpTime);
        
        %%
        savename = [pID,'_',expMode,'_',sesh{seshSel},'_spkDataRespLockedSegmented.mat'];
        save([savepath,sesh{seshSel},'/',savename],'sortedSpikesSEG','wltCoeffs','Fs','chanLab','RTs','trlSel','trlENC','trlRET','hitIdx','missIdx','-v7.3');
        clear sortedSpikesSEG;
        
        %% use fieldtrip to downsample LFP data
        [LFPseg,dsFs] = downsampleLFP(LFPseg,trlTime,1e3);
        [dsTrlTime] = -pre:1/dsFs:post;
        if (size(LFPseg{1},1) ~= length(dsTrlTime)) || (size(LFPseg{1},2)~=length(trlSel))
            error('trial and sample numbers out of range');
        end;
        
        %%
        savename = [pID,'_',expMode,'_',sesh{seshSel},'_lfpDataRespLockedSegmenteddownsampled.mat'];
        save([savepath,sesh{seshSel},'/',savename],'LFPseg','dsTrlTime','dsFs','chanLab','RTs','trlSel','trlENC','trlRET','hitIdx','missIdx','-v7.3');
        
        %%
        savename = [pID,'_',expMode,'_',sesh{seshSel},'_spkDataContinuous.mat'];
        save([savepath,sesh{seshSel},'/',savename],'sortedSpikes','wltCoeffs','Fs','chanLab','-v7.3');
        clear sortedSpikes wltCoeffs;
        
        %% use fieldtrip to downsample LFP data
        for it = 1:length(LFPsig)
            if size(LFPsig{it},2)>size(LFPsig{it},1)
                LFPsig{it} = LFPsig{it}';
            end;
        end;
        
        [LFPsig,dsFs] = downsampleLFP(LFPsig,lfpTime,1e3);
        [dsLfpTime] = [0:size(LFPsig{1},1)-1]./dsFs;
        
        %%
        savename = [pID,'_',expMode,'_',sesh{seshSel},'_lfpDataContinuousDownsampled.mat'];
        save([savepath,sesh{seshSel},'/',savename],'LFPsig','dsLfpTime','dsFs','chanLab','-v7.3');
        clear LFPsig dsLfpTime dsFs chanLab;
    end;
end;
%% shut down parpool
delete(gcp);
exit;