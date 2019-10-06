%%
p2dat = '/media/rouxf/rds-share/resultsAUG2019/';

[ fN ] = dir([ p2dat, '*_spk2LFPCouplingHitsANDmisses_ppc_Sorting_alpha0.05_nRand400_stratModeOn_lowFreq.mat' ]);

plvH = [];
plvM = [];
chanCnt = 0;
for curFile = 1:length( fN )
    
    tmpDat = load([ p2dat, fN(curFile).name ]);
    
    for curChan = 1:size(tmpDat.spk2LFPCouplingH,1 )
        chanCnt = chanCnt +1;
        plvH(chanCnt,:) = tmpDat.spk2LFPCouplingH(curChan,:);
        plvM(chanCnt,:) = tmpDat.spk2LFPCouplingM(curChan,:);
    end;
    
end;

%%
ix = find( tmpDat.spk2LFPfreqAx >= 0 & tmpDat.spk2LFPfreqAx<=30);
figure;
subplot(121);
hold on;
plot(tmpDat.spk2LFPfreqAx(ix), nanmean(plvH,1), 'LineWidth', 3);
plot( tmpDat.spk2LFPfreqAx(ix), nanmean(plvM,1),'r', 'LineWidth',3);
subplot(122);
plot(tmpDat.spk2LFPfreqAx(ix), nanmean(plvH-plvM,1), 'LineWidth', 3);
