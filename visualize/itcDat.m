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
[ poolITC ] = zeros(chnCnt,100,5501);
[ poolERP ] = zeros(chnCnt,5501);

%%
chnCnt = 0;
for curFile = 1:length( spk2lfpFiles )
    
    fprintf( spk2lfpFiles(curFile).name );
    
    [spk2lfpDat] = load( [p2dat, spk2lfpFiles(curFile).name] );
    
    if ~isempty( spk2lfpDat.sigIxLFP )
        ix = regexp(spk2lfpFiles(curFile).name,'_');
        pID = spk2lfpFiles(curFile).name( 1:ix(1)-1 );
        expMode = spk2lfpFiles(curFile).name( ix(1)+1:ix(2)-1 );
        sesh = spk2lfpFiles(curFile).name( ix(2)+1:ix(4)-1 );
        
        [itcFile] = dir([ p2dat,pID,'_',expMode,'_',sesh,'_spectrogramTrialSegITC_lowFreq.mat' ]);
        [itcDat] = load( [ p2dat, itcFile.name ], 'itc', 'itcTime', 'itcFreq','chanLabLFP' );
        
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
            poolERP(chnCnt,:) = erpDat.erp{ sigIxLFP(curChan) };
        end;
    end;
    fprintf('\n');
    
end;

%%
tAx = itcDat.itcTime(itcDat.itcTime>=-0.5 & itcDat.itcTime<=5);
fAx = itcDat.itcFreq;

figure;
pcolor(tAx,fAx,squeeze(mean(poolITC,1)));
shading interp;colormap jet;
ylim([1 30]);