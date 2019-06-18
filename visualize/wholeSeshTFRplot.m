%%
p2d = '/home/rouxf/resultsSpikeFieldOct18II/';
savePath = '/media/rouxf/rds-share/Fred/figures4Cali2018/';
pId = {'P02','P04','P05','P07','P08','P09','P03ERL','P23AMS','P22AMS'};
expMode = {'fVSpEM','cnEM'};

selIx = {};
selIx{1} = [17 40 47];
selIx{2} = [10 27 28 31];
selIx{3} = [29 30 33 39 41];
selIx{4} = [29 30];
%selIx{4,2} =[43 47 50 51 62 64];
selIx{5} = [18];
selIx{6} = [1 3 9 10 25 32 33 34 35 36 37 38 39 40 41:47];
selIx{7} = [1 3:8];
selIx{8} = [4 19];
selIx{9} = [17];

%%
avgPow = cell(1,length( pId ));
lab = cell(1,length(pId));
for curPat = 6%1:length( pId )
    for curExp = 1%:length( expMode )
        
        fN = dir([p2d,pId{curPat},'_',expMode{curExp},'_*_tfrWholeRun_lowFreq.mat']);
        
        if ~isempty(fN)
            
            ix = regexp(fN(1).name,'\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}');
            sesh = fN(1).name(ix:ix+18)
            load([p2d,fN(1).name]);
            
            p2d2 = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pId{curPat},'/',expMode{curExp},'/',sesh,'/'];
            fN2 = dir([p2d2,pId{curPat},'_',expMode{curExp},'_',sesh,'_lfpDataStimLockedContinuousDownsampled.mat']);
            fN2 = fN2.name;
            load([p2d2,fN2]);
            dsLfpTime = [1:size(LFPsig{1},2)]./dsFs;
            dsLfpTime = dsLfpTime';
            
            %%    
            if (curPat == 6)
                
                figure;
                cnt = 0;
                a = [];
                yL = [];
                for curMW = 1:6:length(selIx{curPat})
                    x = LFPsig{selIx{curPat}(curMW)}-ones(1,length(LFPsig{selIx{curPat}(curMW)}))*mean(LFPsig{selIx{curPat}(curMW)});
                    cnt = cnt+1;
                    subplot(4,1,cnt);   
                    a(cnt) = gca;
                    hold on;
                    plot(dsLfpTime,x,'Color',[.5 .5 .5],'LineWidth',2);
                    axis tight;
                    xlim([2588 2599]);
                    axis off;
                    yL(cnt,:) = get(gca,'YLim');                   
                end;
                plot(a(4),[2599-1 2599],[min(min(yL)) min(min(yL))],'k','LineWidth',3);
                plot(a(4),[2599 2599],[min(min(yL)) min(min(yL))+250],'k','LineWidth',3);
                
                set(a,'YLim',[min(min(yL)) max(max(yL))]);
                set(gcf,'Color','w');
                
                saveas(gcf,[savePath,pId{curPat},'_',expMode{curExp},'_',sesh,'_rawLFPtraces.fig']);
                
            end;
            %%
%             for curChan = 1:length(SxxICsL)
%                 figure;
%                 imagesc(txL,fxL,20*log10(SxxICsL{curChan})');
%                 axis xy;
%             end;
                      
            %%               
            avgPow{curPat} = zeros(length(txL),length(fxL));
            for curChan = 1:length(selIx{curPat})
                avgPow{curPat} = avgPow{curPat} + SxxICsL{selIx{curPat}(curChan)};
            end;
            avgPow{curPat} = avgPow{curPat}./length(selIx{curPat});
            
            figure;
            subplot(1,4,1:3);
            aX = gca;
            imagesc(txL./60,fxL,20*log10(avgPow{curPat})');
            axis xy;
            colormap jet;caxis([-30 0]);
            lab{curPat} = [pId{curPat},'_',expMode{curExp},'_',sesh];
            lab{curPat}(regexp(lab{curPat},'_')) = '-';
            text(txL(1)./60,31,lab{curPat},'FontWeight','bold');
            
            set(aX,'Layer','Top');
            for curAx = 1:length(aX)
                box(aX(curAx),'on');
            end;
            ylabel(aX(1),'Frequency [hz]');
            %xlabel(aX(1),'Time [min]');
            %xlabel(aX(2),'LFP-Power [db]');
            for curAx = 1:length(aX)
                xtl = get(aX(curAx),'XLim');
                set(aX(curAx),'XTick',[min(round(xtl))+1:10:max(floor(xtl))-1]);
                set(aX(curAx),'XTickLabel',repmat('',[1 length(get(gca,'XTick'))]));
            end;
            set(aX,'YTick',[5 15 25]);
            set(aX,'Fontsize',14);
            set(aX,'LineWidth',2);            
            set(aX,'Fontsize',14);
            set(gcf,'Color','w');
                
            saveas(gcf,[savePath,pId{curPat},'_',expMode{curExp},'_',sesh,'_tfrWholeRun_lowFreq.fig']);
            
        end;
    end
end;

%%
X = zeros(1,length(fxL));
cnt = 0;
figure;
hold on;
ax = [];
X = [];
for curPat = 1:length(avgPow)    
    y = mean(avgPow{curPat},1);
    if ~isempty(y)        
        cnt = cnt+1;
        subplot(3,3,cnt);
        ax(cnt) = gca;
        y = (y-min(y))./(max(y)-min(y));
        X(cnt,:) = y;
        plot(fxL,y,'c','LineWidth',2);        
        tmp = [pId{curPat},'_',expMode{curExp}];
        tmp(regexp(tmp,'_')) = [];
        title(tmp);
    end;
end;
subplot(3,3,cnt+1);
ax = [ax gca];
SE = std(X,0,1)./sqrt(size(X,1)-1);
jbfill(fxL,mean(X,1)-SE,mean(X,1)+SE,[.75 .75 .75],[.75 .75 .75],1,1);
hold on;
plot(fxL,mean(X,1),'Color',[0 0 0],'LineWidth',2);

set(ax,'Fontsize',14);
set(ax,'LineWidth',2);
for curAx = 1:length( ax )
    axis(ax(curAx),'tight');   
    set(ax(curAx),'XTick',[5 15 25]);
    set(ax(curAx),'XTickLabel',{'','',''});
    if ismember(curAx,[1 4 7])
        ylabel(ax(curAx),'LFP-Power [a.u.]');    
    end;
    set(ax(curAx),'YTick',[min(round(get(ax(curAx),'YLim')*10)/10)+.1 max(floor(get(ax(end),'YLim')*10)/10)-.1]);
    box(ax(curAx),'on');    
end;
xlabel(ax(end),'Frequency [hz]');
set(ax(end),'XTickLabel',{'5','15','25'});
box(ax(end),'on');
axis(ax(end),'tight');
set(gcf,'Color','w');

saveas(gcf,[savePath,expMode{curExp},'_',sesh,'_tfrWholeRun_lowFreq_multiPanelPlotSpectra.fig']);

%%
figure;
subplot(221);
cb = colorbar;
zlab = get(cb,'YLabel');
set(zlab,'String','LFP-Power [dB]');
colormap jet;
caxis([-30 0]);
axis off;
set(gcf,'Color','w');
set(cb,'Fontsize',14);
set(cb,'YTick',[min(get(cb,'YLim')) max(get(cb,'YLim'))]);

saveas(gcf,[savePath,pId{curPat},'_',expMode{curExp},'_',sesh,'_tfrWholeRun_lowFreqColorbar.fig']);