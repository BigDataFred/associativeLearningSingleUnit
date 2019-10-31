%%
addpath('~/associativeLearningSingleUnit/helper/');
p2dat = '/media/rouxf/rds-share/resultsAUG2019/';

%%
spk2lfpMode = {'plv'};

%%
tw = [-.5 0;...
    0 2;...
    2 3;...
    3 5];

twTitle = {'Baseline','Cue','Association','Response'};

%%
cntPlt=0;
figure;
for  curSpk2LFPmode = 1:length( spk2lfpMode )
    for curTW = 1:size( tw,1 )
        
        [ fN ] = dir([ p2dat, '*_spk2LFPCouplingHitsANDmisses_',spk2lfpMode{curSpk2LFPmode},'_Sorting_alpha0.05_nRand400_stratModeOntimeWindow',num2str(tw(curTW,1)),':',num2str(tw(curTW,2)),'s_lowFreq.mat' ]);
        
        plvH = [];
        plvM = [];
        chanCnt = 0;
        for curFile = 1:length( fN )
            
            tmpDat = load([ p2dat, fN(curFile).name ]);
            
            for curChan = 1:size(tmpDat.spk2LFPCouplingH,1 )
                if all(~isnan(tmpDat.spk2LFPCouplingH(curChan,:)))
                    if ~(all(diff(tmpDat.spk2LFPCouplingH(curChan,:))==0)) && ~(all(diff(tmpDat.spk2LFPCouplingM(curChan,:))==0))
                        chanCnt = chanCnt +1;
                        plvH(chanCnt,:) = tmpDat.spk2LFPCouplingH(curChan,:);
                        plvM(chanCnt,:) = tmpDat.spk2LFPCouplingM(curChan,:);
                    else
                        fprintf('TOTO');
                    end;
                end;
            end;
            
        end;
        
        %%
        ix = find( tmpDat.spk2LFPfreqAx >= 0 & tmpDat.spk2LFPfreqAx<=30);
        
        cntPlt = cntPlt+1;
        subplot(4,2,cntPlt);
        hold on;
        %plot(tmpDat.spk2LFPfreqAx(ix), plvH, 'Color',[0 0 1],'LineWidth', 1);
        %plot(tmpDat.spk2LFPfreqAx(ix), plvM, 'Color',[1 0 0],'LineWidth', 1);
        plot( tmpDat.spk2LFPfreqAx(ix), nanmean(plvH,1), 'Color',[0 0 1],'LineWidth', 3);
        plot( tmpDat.spk2LFPfreqAx(ix), nanmean(plvM,1),'r', 'LineWidth',3);
        ylabel([spk2lfpMode{curSpk2LFPmode},' [a.u.]']);
        xlabel('Frequency [Hz]');
        legend('Hits','Misses');
        title(twTitle(curTW));
        
        cntPlt = cntPlt+1;
        subplot(4,2,cntPlt);
        hold on;
        SE = nanstd( plvH-plvM,0,1)/sqrt(size(plvH,1)-1);
        M = nanmean( plvH-plvM,1);
        jbfill(tmpDat.spk2LFPfreqAx(ix),M-SE,M+SE,[0 0 1],[0 0 1],1,.5);
        hold on;
        plot(tmpDat.spk2LFPfreqAx(ix), nanmean(plvH-plvM,1), 'LineWidth', 3);
        plot([0 30],[0 0],'k--','LineWidth',3);
        ylabel(['Hits-Misses [',spk2lfpMode{curSpk2LFPmode},']']);
        xlabel('Frequency [Hz]');
        
    end;
end;