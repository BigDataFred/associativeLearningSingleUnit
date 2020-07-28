%%
winlin=1;% 1= windows PC; 2 = Linux
if winlin == 1
addpath('V:\fred4simon\associativeLearningSingleUnit\helper');
addpath('V:\fred4simon\associativeLearningSingleUnit\preproc');
addpath('V:\fred4simon\associativeLearningSingleUnit\spectral');
addpath('V:\fred4simon\associativeLearningSingleUnit\visualize');
%addpath(genpath('V:\fred4simon\tbx\chronux_2_11\spectral_analysis\'));
addpath('V:\fred4simon\tbx\fieldtrip-master\');
addpath('V:\fred4simon\tbx\cbrewer\cbrewer\');
%ft_defaults;
p2dat = 'V:\resultsAUG2019\SPKs_filt\';
else
addpath('/media/hanslmas/rds-share/fred4simon/associativeLearningSingleUnit/helper');
addpath('/media/hanslmas/rds-share/fred4simon/associativeLearningSingleUnit/preproc');
addpath('/media/hanslmas/rds-share/fred4simon/associativeLearningSingleUnit/spectral');
addpath('/media/hanslmas/rds-share/fred4simon/associativeLearningSingleUnit/visualize');
%addpath(genpath('/media/hanslmas/rds-share/fred4simon/tbx/chronux_2_11/spectral_analysis/'));
addpath('/media/hanslmas/rds-share/fred4simon/tbx/fieldtrip-master/');
addpath('/media/hanslmas/rds-share/fred4simon/tbx/cbrewer/cbrewer/');
p2dat = '/media/hanslmas/rds-share/resultsAUG2019/SPKs_filt/';
end

[ fN ] = dir([ p2dat, 'CFC_sumDat_allCombs_freqadj4*.mat' ]);
load([p2dat,fN.name]);

% Plot Phase-Locked Gamma Power for ALL, Hits and Misses
cmap=cbrewer('div','RdBu',100);
ndat=size(PACs.pow,3);
HFpow=mean(PACs.pow,3);
HFpowH=mean(PACs.powH,3);
HFpowM=mean(PACs.powM,3);

MI=mean(PACs.mi,3);
MIH=mean(PACs.miH,3);
MIM=mean(PACs.miM,3);

phs=mean(PACs.phs,1);
lfp=mean(PACs.lfp,1);
fAxhf=[-16:1:16];
tAx=[-250:1:250];
fAxlf=[-3:1:3];
clim1=[-0.04 0.04];

figure;
subplot(10,1,1:7);
pcolor(tAx,fAxhf,HFpow);caxis(clim1);
shading interp;
title('CFC All');
colormap(cmap)
ax=gca;
ax.XTickLabel={};
subplot(10,1,8:10);
plot(tAx,phs,'r',tAx,lfp,'b');legend phase lfp;


[tmskd,msk]=ttestmat_fdr(PACs.mi,'FDR');

clim2=[0 1.8];
%cmap2=flipud(cmap(1:50,:));
cmap2=flipud(cmap);
figure;
subplot(2,3,1);
pcolor(fAxlf,fAxhf,MI);caxis(clim2);
shading interp
colormap(cmap2);
colorbar
title('Tort MI All (z-val)');
hold on

windy=[-5:5];
windx=[-1:1];
% Plot box to highlight dataselection
xbox=[windx(1,1) windx(end) windx(end) windx(1,1) windx(1,1)];
ybox=[windy(1,1) windy(1,1) windy(end) windy(end) windy(1,1)];
plot(xbox,ybox,'Color',[1 1 1]);
hold off


clim2=[0 1.5];
subplot(2,3,2);
pcolor(fAxlf,fAxhf,MIH);caxis(clim2);
shading interp
colormap(cmap2);
colorbar
title('Tort MI Hits (z-val)');

subplot(2,3,3);
pcolor(fAxlf,fAxhf,MIM);caxis(clim2);
shading interp
colormap(cmap2);
colorbar
title('Tort MI Miss (z-val)');

clim3=[-1 1];
subplot(2,3,4);
pcolor(fAxlf,fAxhf,MIH-MIM);caxis(clim3);
shading interp
%colormap(cmap2);
colorbar
title('Tort MI Diff');
hold on
windy=[-5:5];
windx=[-1:1];
% Plot box to highlight dataselection
xbox=[windx(1,1) windx(end) windx(end) windx(1,1) windx(1,1)];
ybox=[windy(1,1) windy(1,1) windy(end) windy(end) windy(1,1)];
plot(xbox,ybox,'Color',[1 1 1]);
hold off

df=size(PACs.mi,3)-1;
tcrit=tinv(0.95,df);
clim4=[tcrit 3.5];

[tmskd,msk]=ttestmat_fdr(PACs.miH-PACs.miM,'FDR');

subplot(2,3,5);
pcolor(fAxlf,fAxhf,tmskd);caxis(clim4);
shading interp
%colormap(cmap2);
colorbar
title('Tort MI Stats Tvals');
hold on

datH=PACs.miH(17+windy(1,:),4+windx,:);
datH=squeeze(mean(mean(datH,1),2));
datM=PACs.miM(17+windy(1,:),4+windx,:);
datM=squeeze(mean(mean(datM,1),2));
dat=datH-datM;
[~,p,~,stat]=ttest(dat,0,'tail','right');
pwilcx=signrank(dat,0,'tail','right');

subplot(2,3,6);
figure;
vert_rain_plot([datH datM],[0 0 1;1 0 0],50,150);
title(['Hits vs Misses Tval=', num2str(stat.tstat), '---p=',num2str(p)]); 

%% Plot Individual Datasets

figure;
for n=1:ns
    
    subplot(2,3,1);
    pcolor(tAx,fAxhf,PACs.pow(:,:,n));
    shading interp;
    title('CFC All');
    colormap(cmap)
    
    subplot(2,3,2);
    pcolor(tAx,fAxhf,PACs.powH(:,:,n));
    shading interp;
    title('CFC Hits');
    colormap(cmap)
    
    subplot(2,3,3);
    pcolor(tAx,fAxhf,PACs.powM(:,:,n));
    shading interp;
    title('CFC Misses');
    colormap(cmap)
    
    % 2D peak detection
    tmpmi=PACs.mi(:,:,n);
    [pks,locs] = findpeaks(tmpmi(:), 'MinPeakProminence',2.1);
    maxpk=max(pks);
    pkidx=find(pks==maxpk);
    [r,c] = ind2sub(size(tmpmi), locs);
    cent=[r(pkidx) c(pkidx)];
    
    subplot(2,3,4);
    clim4=[0 max(max(PACs.mi(:,:,n)))];
    pcolor(fAxlf,fAxhf,PACs.mi(:,:,n));caxis(clim4);
    shading interp
    colormap(cmap2);
    colorbar
    title('MI all');
    hold on
    plot(fAxlf(cent(1,2)),fAxhf(cent(1,1)),'o', 'MarkerEdgeColor', [1 1 1]);
    hold off
    
    % 2D peak detection
    tmpmi=PACs.miH(:,:,n);
    [pks,locs] = findpeaks(tmpmi(:), 'MinPeakProminence',2.1);
    maxpk=max(pks);
    pkidx=find(pks==maxpk);
    [r,c] = ind2sub(size(tmpmi), locs);
    cent=[r(pkidx) c(pkidx)];
    
    subplot(2,3,5);
    pcolor(fAxlf,fAxhf,PACs.miH(:,:,n));caxis(clim4);
    shading interp
    colormap(cmap2);
    colorbar
    title('MI Hits');
    hold on
    plot(fAxlf(cent(1,2)),fAxhf(cent(1,1)),'o', 'MarkerEdgeColor', [1 1 1]);
    hold off
    
    % 2D peak detection
    tmpmi=PACs.miM(:,:,n);
    [pks,locs] = findpeaks(tmpmi(:), 'MinPeakProminence',2.1);
    maxpk=max(pks);
    pkidx=find(pks==maxpk);
    [r,c] = ind2sub(size(tmpmi), locs);
    cent=[r(pkidx) c(pkidx)];
    
    subplot(2,3,6);
    pcolor(fAxlf,fAxhf,PACs.miM(:,:,n));caxis(clim4);
    shading interp
    %colormap(cmap2);
    colorbar
    title('MI Miss');
    hold on
    plot(fAxlf(cent(1,2)),fAxhf(cent(1,1)),'o', 'MarkerEdgeColor', [1 1 1]);
    hold off
    
end
%% 2D Peak Detection

ns=size(PACs.mi,3);

zcrit=icdf('normal',0.99,0,1);
brdmsk=ones(size(MI,1),size(MI,2));
%brdmsk(:,1)=0;brdmsk(:,end)=0;
%brdmsk(1,:)=0;brdmsk(end,:)=0;
cnt=0;
figure;
for n=1:ns
    tmpmiH=PACs.miH(:,:,n).*brdmsk;
    tmpmiM=PACs.miM(:,:,n).*brdmsk;
    
    % 2D peak detection
    [pksH,locsH] = findpeaks(tmpmiH(:), 'MinPeakProminence',2.1);
    [pksM,locsM] = findpeaks(tmpmiM(:), 'MinPeakProminence',2.1);
    if ~isempty(pksH) & ~isempty(pksM)
        maxpkH=max(pksH);
        maxpkM=max(pksM);
        if maxpkM>=zcrit && maxpkH>=zcrit
            cnt=cnt+1;
            pkidxH=find(pksH==maxpkH);
            [rH,cH] = ind2sub(size(tmpmiH), locsH);
            centH=[rH(pkidxH) cH(pkidxH)];
            
            plot(fAxlf(centH(1,2)),fAxhf(centH(1,1)),'o', 'MarkerEdgeColor', [0 0 1], 'MarkerSize', maxpkH);
            hold on
            
            pkidxM=find(pksM==maxpkM);
            [rM,cM] = ind2sub(size(tmpmiM), locsM);
            centM=[rM(pkidxM) cM(pkidxM)];
            
            plot(fAxlf(centM(1,2)),fAxhf(centM(1,1)),'o', 'MarkerEdgeColor', [1 0 0], 'MarkerSize', maxpkM);
            hold on
            FcentH(cnt,1)=fAxlf(centH(1,2));
            FcentH(cnt,2)=fAxhf(centH(1,1));
            FcentM(cnt,1)=fAxlf(centM(1,2));
            FcentM(cnt,2)=fAxhf(centM(1,1));
            Fdiffscr(cnt,1)=(FcentH(cnt,1)-FcentM(cnt,1))*(FcentH(cnt,2)-FcentM(cnt,2));
            MIpk(cnt,1)=maxpkH;
            MIpk(cnt,2)=maxpkM;
            
        end
    end
end

xlim([-3 3]);
ylim([-16 16]);

pkfLF=FcentH(:,1)-FcentM(:,1);
pkfHF=FcentH(:,2)-FcentM(:,2);

[~,p1,~,stat1]=ttest(pkfLF,0,'tail','right');
[~,p2,~,stat2]=ttest(pkfHF,0,'tail','right');
[~,p3,~,stat3]=ttest(MIpk(:,1)-MIpk(:,2),0,'tail','right');
[~,p4,~,stat4]=ttest(Fdiffscr,0,'tail','right');

figure;
subplot(2,2,1);
horz_rain_plot([FcentH(:,1) FcentM(:,1)],[0 0 1;1 0 0],50,150);
title(['Peak Theta Freq Hits and Misses']); 
subplot(2,2,3);
horz_rain_plot(FcentH(:,1)-FcentM(:,1),[],50,150);
title(['Peak Theta Freq Diff Theta', num2str(stat1.tstat), '---p=',num2str(p1)]); 

subplot(2,2,2);
horz_rain_plot([FcentH(:,2) FcentM(:,2)],[0 0 1;1 0 0],50,150);
title(['Peak Gamma Freq Hits and Misses']); 
subplot(2,2,4);
horz_rain_plot(FcentH(:,2)-FcentM(:,2),[],50,150);
title(['Peak Gamma Freq Diff Theta', num2str(stat2.tstat), '---p=',num2str(p2)]); 

figure;
vert_rain_plot([MIpk(:,1) MIpk(:,2)],[0 0 1;1 0 0],50,150);
title(['Hits vs Misses Tval=', num2str(stat3.tstat), '---p=',num2str(p3)]); 

figure;
vert_rain_plot(MIpk(:,1)-MIpk(:,2),[],50,350);
title(['Peak Difference Hits vs Miss', num2str(stat4.tstat), '---p=',num2str(p3)]);
ylim([-10 15]);

%% Subfunction for stats

function [tmskd,msk]=ttestmat_fdr(dat,mcctype)

for hf=1:size(dat,1)
    for lf=1:size(dat,2)
        difmat=squeeze(dat(hf,lf,:));
        [~,p(hf,lf),~,stats]=ttest(difmat,0,'tail','right');
        tval(hf,lf)=stats.tstat;
    end
end

msk=zeros(hf,lf);

switch mcctype
    case 'FDR' 
        [~,padj,~]=fdr(p(:));
        msk(find(padj<=0.05))=1;
        tmskd=tval.*msk;
    case 'None'
        msk(find(p<=0.05))=1;
        tmskd=tval.*msk;
    otherwise
        warning('Unexpected MCC type. MCCtype can only be FDR or None')
end

end

