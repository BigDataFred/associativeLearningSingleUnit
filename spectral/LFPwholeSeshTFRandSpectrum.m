%%
addpath(genpath('~/tbx/chronux_2_11/spectral_analysis/'));
addpath('~/tbx/fieldtrip-20170618/');
ft_defaults;

%%
pID = {'P02','P04','P05','P07','P08','P09','P03ERL','P22AMS','P23AMS'};%,

%%
expMode = {'fVSpEM','cnEM'};

%%
for curPat = 1:length( pID )
    for curExp = 1:length( expMode )
        
        
        p2d = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID{curPat},'/',expMode{curExp},'/'];
        savePath = ['/home/rouxf/resultsSpikeFieldOct18/'];
        
        chck = dir( savePath );
        if isempty( chck )
            mkdir( savePath );
        end;
        
        sesh = dir(p2d);
        if ~isempty(sesh)
            sesh(1:2) = [];
            sesh
            
            %%
            for curSes = 1:length( sesh )
                
                fN = dir([p2d,sesh(curSes).name,'/',pID{curPat},'_',expMode{curExp},'_',sesh(curSes).name,'_lfpDataContinuousDownsampled.mat']);
                
                if ~isempty( fN )
                    
                    [ lfpDat ] = load([p2d,sesh(curSes).name,'/',fN.name]);
                    
                    %%
                    dum = [];
                    dum.label = lfpDat.chanLab;
                    dum.fsample = lfpDat.dsFs;
                    dum.time = {};
                    dum.trial = {};
                    for curChan = 1:length( lfpDat.chanLab )
                        dum.trial{1}(curChan,:) = lfpDat.LFPsig{curChan};
                    end;
                    dum.time{1} = 0:1/dum.fsample:(length(dum.trial{1})-1)/dum.fsample;
                    
                    cfg                     = [];
                    cfg.resamplefs          = 200;
                    
                    [dum] = ft_resampledata( cfg, dum);
                    
                    %%
                    BFlab = cell(1,length(dum.label));
                    for it = 1:length( dum.label )
                        BFlab(it) = {dum.label{it}(1:regexp(dum.label{it},'\d{1}')-1)};
                    end;
                    BFid = unique(BFlab);
                    
                    selIx = cell(1,length( BFid ));
                    for it = 1:length( BFid )
                        selIx{it} = find(strcmp(BFlab,BFid(it)));
                    end;
                    
%                     %% orthogonalize the single MW-lfp against average of BF-lfp
                    [ ntrl ] = length( dum.trial );
                    [ nsmp ] = size( dum.trial{1},2 );
%                     
%                     for curMW = 1:length( BFid )
%                         
%                         avgLFP = zeros(nsmp,ntrl);
%                         for jt = 1:length( selIx{curMW} )
%                             avgLFP = avgLFP  + dum.trial{1}(selIx{curMW}(jt),:)';
%                         end;
%                         avgLFP = avgLFP./length( selIx{curMW} );
%                         
%                         X = avgLFP';
%                         for kt = 1:length( selIx{curMW} )
%                             Y = dum.trial{1}(selIx{curMW}(kt),:);
%                             [Yp] = real(sum(conj(X.*Y),2)/sum(conj(X.*X),2))*X;
%                             [Yo] = Y - Yp;
%                             %Yo = Y-avgLFP';
%                             dum.trial{1}(selIx{curMW}(kt),:) = Yo;
%                         end;
%                         
%                     end;
                    
                    %%
                    movingwin1 = [4 .25];
                    
                    params1                  = [];
                    params1.Fs               = dum.fsample;
                    params1.pad              = 2;
                    params1.fpass            = [3 30];
                    params1.tapers           = [4 7];
                    params1.trialave         = 0;
                    
                    movingwin2 = [.25 .067];
                    
                    params3                  = [];
                    params3.Fs               = dum.fsample;
                    params3.pad              = 2;
                    params3.fpass            = [30 170];
                    params3.tapers           = [2 3];
                    params3.trialave         = 0;
                    
                    T = (nsmp/dum.fsample);
                    W = 0.1;
                    TW = round(T*W);
                    k = 2*TW-1;
                    
                    params2                  = [];
                    params2.Fs               = dum.fsample;
                    params2.pad              = 0;
                    params2.fpass            = [3 30];
                    params2.tapers           = [TW k];
                    params2.trialave         = 0;
                    
                    T = (nsmp/dum.fsample);
                    W = 4;
                    TW = round(T*W);
                    k = 2*TW-1;
                    
                    params4                  = [];
                    params4.Fs               = dum.fsample;
                    params4.pad              = 0;
                    params4.fpass            = [30 170];
                    params4.tapers           = [TW k];
                    params4.trialave         = 0;
                    
                    %%
                    Slf = cell( 1, length( lfpDat.LFPsig ) );
                    Shf = cell( 1, length( lfpDat.LFPsig ) );
                    %             Sflx = cell( 1, length( lfpDat.LFPsig ) );
                    %             Sfhx = cell( 1, length( lfpDat.LFPsig ) );
                    
                    %%
                    for it = 1:length( dum.label )
                        fprintf([num2str(it),'/',num2str(length(dum.label))]);
                        x = dum.trial{1}(it,:);
                        x = (x-mean(x))./std(x);
                        [Slf{it}] = mtspecgramc( x', movingwin1, params1);
                        [Shf{it}] = mtspecgramc( x', movingwin2, params3);
                        %[Sflx{it}] = mtspectrumc( x', params2);
                        %[Sfhx{it}] = mtspectrumc( x', params4);
                        fprintf('\n');
                    end;
                    clear x;
                    
                    %%
                    [~,tx1,fx1] = mtspecgramc( dum.trial{1}(1,:)' ,movingwin1, params1);
                    %             [~,flxx] = mtspectrumc( dum.trial{1}(1,:)' , params2);
                    [~,tx2,fx2] = mtspecgramc( dum.trial{1}(1,:)' ,movingwin2, params3);
                    %             [~,fhxx] = mtspectrumc( dum.trial{1}(1,:)' , params4);
                    
                    %%
                    chanLab = lfpDat.chanLab;
                    
                    saveName = [pID{curPat},'_',expMode{curExp},'_',sesh(curSes).name,'_tfrWholeRun_lowFreq.mat'];
                    save([savePath,saveName],'Slf','tx1','fx1','chanLab','-v7.3');
                    clear Slf tx1 fx1;
                    
                    saveName = [pID{curPat},'_',expMode{curExp},'_',sesh(curSes).name,'_tfrWholeRun_highFreq.mat'];
                    save([savePath,saveName],'Shf','tx2','fx2','chanLab','-v7.3');
                    clear Shf tx2 fx2;
                    
                    %             saveName = [pID{curPat},'_',expMode{curExp},'_',sesh(curSes).name,'_spectrumWholeRun_lowFreq.mat'];
                    %             save([savePath,saveName],'Sflx','flxx','chanLab','-v7.3');
                    %             clear Sflx flxx;
                    
                    %             saveName = [pID{curPat},'_',expMode{curExp},'_',sesh(curSes).name,'_spectrumWholeRun_highFreq.mat'];
                    %             save([savePath,saveName],'Sfhx','fhxx','chanLab','-v7.3');
                    %             clear Sfhx fhxx;
                    
                    %         %%
                    %         for jt =1:length( Slf )
                    %
                    %             %lab = [sesh(curSes).name,fN(jt).name];
                    %             %lab(regexp(lab,'_')) = [];
                    %
                    %             figure;
                    %             subplot(2,4,1:3);
                    %             imagesc(tx1,fx1,(10*log10(Slf{jt}))');
                    %             axis xy;colormap jet;
                    %             %title(lab);
                    %             subplot(2,4,4);
                    %             hold on;
                    %             Y = 10*log10(Sfx{jt})-10*log10(mean(Sfx{jt}));
                    %             plot(Y,fxx);
                    %             %b = regress(Y,[ones(size(fxx')) 10*log10(fxx)']);
                    %             b = robustfit([10*log10(fxx)'],Y);
                    %
                    %             %yp = [ones(size(fxx')) 10*log10(fxx)']*b;
                    %             yp = [ones(size(fxx')) 10*log10(fxx)']*b;
                    %             plot(yp,fxx,'r');
                    %             axis tight;
                    %             subplot(2,4,5:7);
                    %             X = 10*log10(fx1)';
                    %             yc = zeros(size(Slf{jt}));
                    %             for it = 1:length( tx1 )
                    %                 Y2 = 10*log10(Slf{jt}(it,:))';
                    %                 b = regress(Y2,[ones(size(X)) X]);
                    %                 yp2 = [ones(size(X)) X]*b;
                    %                 yc(it,:) = (Y2-yp2);
                    %             end;
                    %             imagesc(tx1,fx1,yc');
                    %             %caxis([0 max(max(10.^(Y-yp)))]);
                    %             axis xy;
                    %             subplot(2,4,8);
                    %             plot((Y-yp),fxx);
                    %             axis tight;
                    %
                    %         end;
                end;
            end;
        end;
    end;
end;
%%


















