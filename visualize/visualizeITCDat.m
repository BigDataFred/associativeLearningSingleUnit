%%
[p2dat] = '/media/rouxf/rds-share/resultsAUG2019/';
[spk2lfpFiles] = dir([p2dat, '*_spk2LFPCoupling_plv_Sorting_alpha0.05_nRand400_lowFreq.mat' ])

%%
chnCnt = 0;
for curFile = 1:length( spk2lfpFiles )
    
    [spk2lfpDat] = load( [p2dat, spk2lfpFiles(curFile).name] );
    if ~isempty( spk2lfpDat.sigIxLFP )
        chnCnt = chnCnt + length( spk2lfpDat.sigIxLFP );
    end;
end;

%%
fS = 1e3;
[b] = fir1(3*fix(fS/30),30/fS);

%%
[ poolITC ] = zeros(chnCnt,100,5501);
[ poolITCh ] = zeros(chnCnt,100,5501);
[ poolITCm ] = zeros(chnCnt,100,5501);
[ poolERP ] = zeros(chnCnt,5501);

%%
chnCnt = 0;
for curFile = 1:length( spk2lfpFiles )
    
    fprintf( spk2lfpFiles(curFile).name );
    
    [spk2lfpDat] = load( [p2dat, spk2lfpFiles(curFile).name] );
    
    if ~isempty( spk2lfpDat.sigIxLFP )
        ix = regexp(spk2lfpFiles(curFile).name,'_');
        pID = spk2lfpFiles(curFile).name( 1:ix(1)-1 );
        %if ~strcmp( pID,'P02')
            expMode = spk2lfpFiles(curFile).name( ix(1)+1:ix(2)-1 );
            sesh = spk2lfpFiles(curFile).name( ix(2)+1:ix(4)-1 );
            
            [itcFile] = dir([ p2dat,pID,'_',expMode,'_',sesh,'_spectrogramTrialSegITC_lowFreq.mat' ]);
            [itcDat] = load( [ p2dat, itcFile.name ] ,'itc','itcTime','itcFreq','chanLabLFP');%'itcH','itcM',
            
            [erpFile] = dir([p2dat,pID,'_',expMode,'_',sesh,'_ERPTrialSeg_HitsANDmisses.mat']);
            [erpDat]  = load( [ p2dat, erpFile.name ], 'erp', 'dsTrlTime');
            
            [ sigIxLFP ] = zeros(length( spk2lfpDat.chanLabLFPsig ),1);
            for curChan = 1:length( spk2lfpDat.chanLabLFPsig )
                sigIxLFP( curChan ) = find(strcmp( itcDat.chanLabLFP, spk2lfpDat.chanLabLFPsig( curChan )));
            end;
            
            if any( sigIxLFP~=spk2lfpDat.sigIxLFP)
                error('!');
            end;
            
            for curChan = 1:length( spk2lfpDat.chanLabLFPsig )
                chnCnt = chnCnt+1;
                poolITC(chnCnt,:,:) = itcDat.itc{ sigIxLFP(curChan) }(:,itcDat.itcTime>=-0.5 & itcDat.itcTime<=5);
                %poolITCh(chnCnt,:,:) = itcDat.itcH{ sigIxLFP(curChan) }(:,itcDat.itcTime>=-0.5 & itcDat.itcTime<=5);
                %poolITCm(chnCnt,:,:) = itcDat.itcM{ sigIxLFP(curChan) }(:,itcDat.itcTime>=-0.5 & itcDat.itcTime<=5);
                poolERP(chnCnt,:) = abs( erpDat.erp{ sigIxLFP(curChan) }(erpDat.dsTrlTime>=-.5 & erpDat.dsTrlTime<=5) ).^2;
            end;
        %end;
    end;
    fprintf('\n');
    
end;

%%
tAx = itcDat.itcTime(itcDat.itcTime>=-0.5 & itcDat.itcTime<=5);
fAx = itcDat.itcFreq;

figure;
subplot(211);
pcolor(tAx,fAx,squeeze(mean(poolITC,1)));
shading interp;colormap jet;ylim([1 30]);
subplot(212);
M = squeeze(mean(mean(poolITC(:,fAx>=4 & fAx<=12,:),2)))';
SE = squeeze(std(mean(poolITC(:,fAx>=4 & fAx<=12,:),2),0,1))';
jbfill(tAx,M-SE,M+SE,[.8 .8 .8],[.8 .8 .8],1,.45);
hold on;
plot(tAx,squeeze(mean(mean(poolITC(:,fAx>=4 & fAx<=12,:),2),1)),'Color',[1 0 0],'LineWidth',3);
axis tight;

