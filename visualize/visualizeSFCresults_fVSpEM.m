%%
addpath('~/prj/Bham/code/mcode/utils/');

%%
p2d = '/home/rouxf/resultsSpikeFieldOct18/';
pId = {'P02','P04','P05','P07','P08','P09','P03ERL'};

cntSPKSrt = 0;
cntSameMW = 0; cntSameBF = 0;cntOther=0;
cntABT = 0; poolSFC1 = [];
for curPat = 1:length( pId )

    fN = dir([p2d,pId{curPat},'_fVSpEM_*_spk2LFPCoupling_plv_lowFreq.mat']);
    [ sesh ] = cell(1,length(fN));
    for curSesh = 1:length(fN)
        ix = regexp(fN(curSesh).name,'\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}');
        sesh(curSesh) = {fN(curSesh).name(ix:ix+18)};
    end;
    
    for curSesh = 1:length(sesh)
        
        fNSesh = dir([p2d,pId{curPat},'_fVSpEM_',sesh{curSesh},'_spk2LFPCoupling_plv_lowFreq.mat'])
        SFCdat = load([p2d,fNSesh.name]);
        x =  SFCdat.spk2LFPCoupling;
        if ~isempty(x)
            lfpBFlab = {};
            for curLFP = 1:size(SFCdat.chanLabLFP)
                ix = regexp(SFCdat.chanLabLFP{curLFP},'\d{1}');
                lfpBFlab(curLFP) = {SFCdat.chanLabLFP{curLFP}(1:ix-1)};
            end;
            lfpBFid = unique(lfpBFlab);
            
            for curMW = 1:length( SFCdat.chanLabSPK )
                ix = regexp(SFCdat.chanLabSPK{curMW},'\d{1}');
                spkBFlab = {SFCdat.chanLabSPK{curMW}(1:ix-1)};
                curMW2 = curMW;
                %for curMW2 = 1:length( SFCdat.chanLabLFP)
                    
                    if strcmp(SFCdat.chanLabSPK{curMW},SFCdat.chanLabLFP{curMW2})
                        cntSameMW = cntSameMW+1;
                    elseif any(strcmp(spkBFlab,lfpBFid))
                        cntSameBF = cntSameBF+1;
                    else
                        cntOther = cntOther+1;
                    end;
                    
                %end;
            end;
            
            for curChan = 1:size( x,1 )
                %for curChan2 = 1:size(x,2)
                    cntABT = cntABT+1;
                    cntSPKSrt = cntSPKSrt+1;
                    %poolSFC1(cntABT,:) =squeeze(x(curChan,curChan2,:));
                    poolSFC1(cntABT,:) =squeeze(x(curChan,:));
                %end;
            end;
        end;
               
    end;
    
end;

%%
cntSPKSrt2 = 0; poolSFC2 = []; spk2LFPCH = []; spk2LFPCM = []; spk2LFPCS = [];
cntSameMW2 = 0; cntSameBF2 = 0;cntOther2=0;
cntSPK2LFPH = 0; cntSPK2LFPM = 0;cntSPK2LFPS = 0;

for curPat = 1:length( pId )
    fN = dir([p2d,pId{curPat},'_fVSpEM_*_spk2LFPCouplingSpikeSorting_plv_lowFreq.mat']);
    [ sesh ] = cell(1,length(fN));
    for curSesh = 1:length(fN)
        ix = regexp(fN(curSesh).name,'\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}');
        sesh(curSesh) = {fN(curSesh).name(ix:ix+18)};
    end;
    
    for curSesh = 1:length(sesh)
        
        fNSesh = dir([p2d,pId{curPat},'_fVSpEM_',sesh{curSesh},'_spk2LFPCouplingSpikeSorting_plv_lowFreq.mat'])     
        SFCdat = load([p2d,fNSesh.name]);
        x =  SFCdat.spk2LFPCoupling;
        if ~isempty(x)
            lfpBFlab = {};
            for curLFP = 1:size(SFCdat.chanLabLFP)
                ix = regexp(SFCdat.chanLabLFP{curLFP},'\d{1}');
                lfpBFlab(curLFP) = {SFCdat.chanLabLFP{curLFP}(1:ix-1)};
            end;
            lfpBFid = unique(lfpBFlab);
            
            for curMW = 1:length( SFCdat.chanLabSPK )
                ix = regexp(SFCdat.chanLabSPK{curMW},'\d{1}');
                spkBFlab = {SFCdat.chanLabSPK{curMW}(1:ix-1)};
                curMW2 = curMW;
                %for curMW2 = 1:length( SFCdat.chanLabLFP)
                    
                    if strcmp(SFCdat.chanLabSPK{curMW},SFCdat.chanLabLFP{curMW2})
                        cntSameMW2 = cntSameMW2+1;
                    elseif any(strcmp(spkBFlab,lfpBFid))
                        cntSameBF2 = cntSameBF2+1;
                    else
                        cntOther2 = cntOther2+1;
                    end;
                    
                %end;
            end;
            
            for curChan = 1:size( x,1 )
                %for curChan2 = 1:size(x,2)
                    cntSPKSrt2 = cntSPKSrt2+1;
                    %poolSFC2(cntSPKSrt2,:) =squeeze(x(curChan,curChan2,:));
                    poolSFC2(cntSPKSrt2,:) =squeeze(x(curChan,:));
                %end;
            end;
        end;
        
        fNHvM = [pId{curPat},'_fVSpEM_',sesh{curSesh},'_spk2LFPCoupling_HitsVsMisses_plv_TW-0.5-5s_lowFreq.mat'];
        fNHvM = dir([p2d,fNHvM]);
        if ~isempty(fNHvM)
            [ hVmDat ] = load([p2d,fNHvM.name]);
            x = hVmDat.spk2LFPCouplingHits;
            for curChan = 1:size( x,1 )
                %for curChan2 = 1:size(x,2)
                    cntSPK2LFPH = cntSPK2LFPH+1;
                    %spk2LFPCH(cntSPK2LFPH,:) =squeeze(x(curChan,curChan2,:));
                    spk2LFPCH(cntSPK2LFPH,:) =squeeze(x(curChan,:));
                %end;
            end;
            
            x = hVmDat.spk2LFPCouplingMiss;
            for curChan = 1:size( x,1 )
                %for curChan2 = 1:size(x,2)
                    cntSPK2LFPM = cntSPK2LFPM+1;
                    %spk2LFPCM(cntSPK2LFPM,:) =squeeze(x(curChan,curChan2,:));
                    spk2LFPCM(cntSPK2LFPM,:) =squeeze(x(curChan,:));
                %end;
            end;
            
            x = hVmDat.spk2LFPCouplingSubs;
            for curChan = 1:size( x,1 )
                %for curChan2 = 1:size(x,2)
                    cntSPK2LFPS = cntSPK2LFPS+1;
                    %spk2LFPCS(cntSPK2LFPS,:) =squeeze(x(curChan,curChan2,:));
                    spk2LFPCS(cntSPK2LFPS,:) =squeeze(x(curChan,:));
                %end;
            end;
            
        end;
    end;
    
end;

%%
figure;

subplot(221);
a = gca;
M1 = mean(poolSFC1,1);
SE1 = std(poolSFC1,0,1)./sqrt(size(poolSFC1,1)-1);
M2 = mean(poolSFC2,1);
SE2 = std(poolSFC2,0,1)./sqrt(size(poolSFC2,1)-1);

jbfill(SFCdat.spk2LFPfreqAx,M2-SE2,M2+SE2,[0 1 0],[0 1 0],1,.5)
hold on;
jbfill(SFCdat.spk2LFPfreqAx,M1-SE1,M1+SE1,[253 106 2]./255,[253 106 2]./255,1,.5)
hold on;
plot(SFCdat.spk2LFPfreqAx,mean(poolSFC2,1),'Color',[.5 .5 .5]./25,'Linewidth',2);
plot(SFCdat.spk2LFPfreqAx,mean(poolSFC1,1),'Color',[0 0 0],'Linewidth',2);
set(gca,'Xtick',[5:10:25]);

subplot(222);
a = [a gca];
hold on;
h = [];
% h(1) = bar(1,(cntABT));
% h(2) = bar(2,(cntSPKSrt));
if (cntSameMW/cntSPKSrt > cntSameMW2/cntSPKSrt2)
    h(1) = bar(1,(cntSameMW)./cntSPKSrt);
    h(2) = bar(1,(cntSameMW2)./cntSPKSrt2);           
else
    h(2) = bar(1,(cntSameMW2)./cntSPKSrt2);
    h(1) = bar(1,(cntSameMW)./cntSPKSrt);        
end;

if (cntSameBF/cntSPKSrt > cntSameBF2/cntSPKSrt2)
    h(3) = bar(2,(cntSameBF)./cntSPKSrt);
    h(4) = bar(2,(cntSameBF2)./cntSPKSrt2);           
else
    h(4) = bar(2,(cntSameBF2)./cntSPKSrt2);
    h(3) = bar(2,(cntSameBF)./cntSPKSrt);        
end;

if (cntOther/cntSPKSrt > cntOther2/cntSPKSrt2)
    h(5) = bar(3,(cntOther)./cntSPKSrt);
    h(6) = bar(3,(cntOther2)./cntSPKSrt2);           
else
    h(6) = bar(3,(cntOther2)./cntSPKSrt2);
    h(5) = bar(3,(cntOther)./cntSPKSrt);        
end;

set(h([1 3 5]),'FaceColor',[253 106 2]./255);
set(h([2 4 6]),'FaceColor',[0 1 0]);

%h(3) = bar(2,(cntSameBF)./cntSPKSrt);
%h(4) = bar(3,(cntOther)./cntSPKSrt);

% set(h(1),'FaceColor',[253 106 2]./255);
% set(h(2),'FaceColor',[0 255 0]./255);
%set(h(1),'FaceColor',[0 0 1]);
%set(h(2),'FaceColor',[0 0 1]);
%set(h(3),'FaceColor',[0 0 1]);

set(gca,'Xtick',[1:3]);
%set(gca,'XTickLabel',{repmat('',[1 length(get(gca,'YTick'))])});
set(gca,'XTickLabel',{'mw','bf','o'});
%set(gca,'XTickLabelRotation',-45);

subplot(223);
M1 = nanmean(spk2LFPCH,1);
SE1 = nanstd(spk2LFPCH,0,1)./sqrt(size(spk2LFPCH,1)-1);
M2 = nanmean(spk2LFPCS,1);
SE2 = nanstd(spk2LFPCS,0,1)./sqrt(size(spk2LFPCS,1)-1);

a = [a gca];
hold on;
jbfill(SFCdat.spk2LFPfreqAx,M2-SE2,M2+SE2,[1 0 0],[1 0 0],1,.5)
hold on;
jbfill(SFCdat.spk2LFPfreqAx,M1-SE1,M1+SE1,[0 0 1],[0 0 1],1,.5)
hold on;
plot(SFCdat.spk2LFPfreqAx,nanmean(spk2LFPCH,1),'Color',[0 0 1],'Linewidth',2);
plot(SFCdat.spk2LFPfreqAx,nanmean(spk2LFPCS,1),'Color',[1 0 0],'Linewidth',2);
set(gca,'Xtick',[5:10:25]);

subplot(224);
a = [a gca];
fx = SFCdat.spk2LFPfreqAx;
pval = zeros(1,length(fx));
h = zeros(1,length(fx));
x = spk2LFPCH-spk2LFPCS;
for curFreq = 1:length( fx )
    
    h(curFreq) = lillietest(squeeze(nanmean(x(:,curFreq),2)));
    %if h(curFreq) ==1
        [pval(curFreq)] = signrank( squeeze(nanmean(x(:,curFreq),2)));
    %else
    %    [~,pval(curFreq)] = ttest( squeeze(nanmean(x(:,curFreq),2)));
    %end;
end;
chck = pval<.05/(length(fx));
fLim= [min(fx(find(chck))) max(fx(find(chck)))];
M = nanmean(x,1);
SE = nanstd(x,0,1)./sqrt(size(x,1)-1);
hold on;
if any(chck)
    h = area([fLim],[max(max([M+SE])) max(max([M+SE]))],min(min([M-SE])),'ShowBaseLine','off');
    set(h,'FaceColor',[.85 .85 .85],'EdgeColor',[.85 .85 .85]);
end;
plot([min(fx) max(fx)],[0 0],'k--','LineWidth',2);
hold on;
jbfill(SFCdat.spk2LFPfreqAx,M-SE,M+SE,[0 0 255]./255,[0 0 255]./255,1,.5)
hold on;
plot(SFCdat.spk2LFPfreqAx,nanmean(x,1),'Color',[0 0 255]./255,'Linewidth',2);
set(gca,'Xtick',[5:10:25]);

xlabel(a(1),'Frequency [hz]');
ylabel(a(1),{'Spike-LFP coupling';'[PLV]'});
xlabel(a(3),'Frequency [hz]');
xlabel(a(4),'Frequency [hz]');
ylabel(a(3),{'Spike-LFP coupling';'[PLV]'});
ylabel(a(2),'% Channels');

for curAx = 1:length(a)
    axis(a(curAx),'tight');
    set(a(curAx),'YTick',[min(get(a(curAx),'YTick')) max(get(a(curAx),'YTick'))]);
    box(a(curAx),'off');
end;
xlim(a(2),[0 4]);

set(a,'Layer','top');
set(a,'Fontsize',14);
set(a,'LineWidth',3);
set(gcf,'Color','w');

%%
spk2LFPCH1 = [];  spk2LFPCH2 = [];  spk2LFPCH3 = [];  spk2LFPCH4 = [];  
spk2LFPCS1 = [];  spk2LFPCS2 = [];  spk2LFPCS3 = [];  spk2LFPCS4 = [];
cntSPK2LFPH1 = 0; cntSPK2LFPH2 = 0; cntSPK2LFPH3 = 0; cntSPK2LFPH4 = 0; 
cntSPK2LFPS1 = 0; cntSPK2LFPS2 = 0; cntSPK2LFPS3 = 0; cntSPK2LFPS4 = 0;

for curPat = 1:length( pId )
    
    fN = dir([p2d,pId{curPat},'_fVSpEM_*_spk2LFPCoupling_plv_lowFreq.mat']);
    [ sesh ] = cell(1,length(fN));
    for curSesh = 1:length(fN)
        ix = regexp(fN(curSesh).name,'\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}');
        sesh(curSesh) = {fN(curSesh).name(ix:ix+18)};
    end;
    
    for curSesh = 1:length(sesh)
        
        fNHvM = [pId{curPat},'_fVSpEM_',sesh{curSesh},'_spk2LFPCoupling_HitsVsMisses_plv_TW-0.5-0s_lowFreq.mat'];
        fNHvM = dir([p2d,fNHvM]);
        if ~isempty(fNHvM)
            [ hVmDat ] = load([p2d,fNHvM.name]);
            x = hVmDat.spk2LFPCouplingHits;
            for curChan = 1:size( x,1 )
                %for curChan2 = 1:size(x,2)
                    cntSPK2LFPH1 = cntSPK2LFPH1+1;
                    %spk2LFPCH1(cntSPK2LFPH1,:) =squeeze(x(curChan,curChan2,:));
                    spk2LFPCH1(cntSPK2LFPH1,:) =squeeze(x(curChan,:));
                %end;
            end;
            
            x = hVmDat.spk2LFPCouplingSubs;
            for curChan = 1:size( x,1 )
                %for curChan2 = 1:size(x,2)
                    cntSPK2LFPS1 = cntSPK2LFPS1+1;
                    %spk2LFPCS1(cntSPK2LFPS1,:) =squeeze(x(curChan,curChan2,:));
                    spk2LFPCS1(cntSPK2LFPS1,:) =squeeze(x(curChan,:));
                %end;
            end;
            
        end;
        
        
        fNHvM = [pId{curPat},'_fVSpEM_',sesh{curSesh},'_spk2LFPCoupling_HitsVsMisses_plv_TW0-2s_lowFreq.mat'];
        fNHvM = dir([p2d,fNHvM]);
        if ~isempty(fNHvM)
            [ hVmDat ] = load([p2d,fNHvM.name]);
            x = hVmDat.spk2LFPCouplingHits;
            for curChan = 1:size( x,1 )
                %for curChan2 = 1:size(x,2)
                    cntSPK2LFPH2 = cntSPK2LFPH2+1;
                    %spk2LFPCH2(cntSPK2LFPH2,:) =squeeze(x(curChan,curChan2,:));
                    spk2LFPCH2(cntSPK2LFPH2,:) =squeeze(x(curChan,:));
                %end;
            end;
            
            x = hVmDat.spk2LFPCouplingSubs;
            for curChan = 1:size( x,1 )
                %for curChan2 = 1:size(x,2)
                    cntSPK2LFPS2 = cntSPK2LFPS2+1;
                    %spk2LFPCS2(cntSPK2LFPS2,:) =squeeze(x(curChan,curChan2,:));
                    spk2LFPCS2(cntSPK2LFPS2,:) =squeeze(x(curChan,:));
                %end;
            end;
            
        end;
        
        fNHvM = [pId{curPat},'_fVSpEM_',sesh{curSesh},'_spk2LFPCoupling_HitsVsMisses_plv_TW2-3s_lowFreq.mat'];
        fNHvM = dir([p2d,fNHvM]);
        if ~isempty(fNHvM)
            [ hVmDat ] = load([p2d,fNHvM.name]);
            x = hVmDat.spk2LFPCouplingHits;
            for curChan = 1:size( x,1 )
                %for curChan2 = 1:size(x,2)
                    cntSPK2LFPH3 = cntSPK2LFPH3+1;
                    %spk2LFPCH3(cntSPK2LFPH3,:) =squeeze(x(curChan,curChan2,:));
                    spk2LFPCH3(cntSPK2LFPH3,:) =squeeze(x(curChan,:));
                %end
            end;
            
            x = hVmDat.spk2LFPCouplingSubs;
            for curChan = 1:size( x,1 )
                %for curChan2 = 1:size(x,2)
                    cntSPK2LFPS3 = cntSPK2LFPS3+1;
                    %spk2LFPCS3(cntSPK2LFPS3,:) =squeeze(x(curChan,curChan2,:));
                    spk2LFPCS3(cntSPK2LFPS3,:) =squeeze(x(curChan,:));
                %end;
            end;
            
        end;
        
        fNHvM = [pId{curPat},'_fVSpEM_',sesh{curSesh},'_spk2LFPCoupling_HitsVsMisses_plv_TW3-5s_lowFreq.mat'];
        fNHvM = dir([p2d,fNHvM]);
        if ~isempty(fNHvM)
            [ hVmDat ] = load([p2d,fNHvM.name]);
            x = hVmDat.spk2LFPCouplingHits;
            for curChan = 1:size( x,1 )
                %for curChan2 = 1:size(x,2)
                    cntSPK2LFPH4 = cntSPK2LFPH4+1;
                    %spk2LFPCH4(cntSPK2LFPH4,:) =squeeze(x(curChan,curChan2,:));
                    spk2LFPCH4(cntSPK2LFPH4,:) =squeeze(x(curChan,:));
                %end;
            end;
            
            x = hVmDat.spk2LFPCouplingSubs;
            for curChan = 1:size( x,1 )
                %for curChan2 = 1:size(x,2)
                    cntSPK2LFPS4 = cntSPK2LFPS4+1;
                    %spk2LFPCS4(cntSPK2LFPS4,:) =squeeze(x(curChan,curChan2,:));
                    spk2LFPCS4(cntSPK2LFPS4,:) =squeeze(x(curChan,:));
                %end;
            end;
            
        end;
    end;
end;

%%
figure;

subplot(4,2,1);
M1 = nanmean(spk2LFPCH1,1);
SE1 = nanstd(spk2LFPCH1,0,1)./sqrt(size(spk2LFPCH1,1)-1);
M2 = nanmean(spk2LFPCS1,1);
SE2 = nanstd(spk2LFPCS1,0,1)./sqrt(size(spk2LFPCS1,1)-1);
a = [gca];
jbfill(SFCdat.spk2LFPfreqAx,M2-SE2,M2+SE2,[1 0 0],[1 0 0],1,.5)
hold on;
jbfill(SFCdat.spk2LFPfreqAx,M1-SE1,M1+SE1,[0 0 1],[0 0 1],1,.5)
hold on;
plot(SFCdat.spk2LFPfreqAx,nanmean(spk2LFPCH1,1),'Color',[0 0 1],'Linewidth',2);
plot(SFCdat.spk2LFPfreqAx,nanmean(spk2LFPCS1,1),'Color',[1 0 0],'Linewidth',2);
set(gca,'Xtick',[5:10:25]);

subplot(4,2,3);
M1 = nanmean(spk2LFPCH2,1);
SE1 = nanstd(spk2LFPCH2,0,1)./sqrt(size(spk2LFPCH2,1)-1);
M2 = nanmean(spk2LFPCS2,1);
SE2 = nanstd(spk2LFPCS2,0,1)./sqrt(size(spk2LFPCS2,1)-1);
a = [a gca];
jbfill(SFCdat.spk2LFPfreqAx,M2-SE2,M2+SE2,[1 0 0],[1 0 0],1,.5)
hold on;
jbfill(SFCdat.spk2LFPfreqAx,M1-SE1,M1+SE1,[0 0 1],[0 0 1],1,.5)
hold on;
plot(SFCdat.spk2LFPfreqAx,nanmean(spk2LFPCH2,1),'Color',[0 0 1],'Linewidth',2);
plot(SFCdat.spk2LFPfreqAx,nanmean(spk2LFPCS2,1),'Color',[1 0 0],'Linewidth',2);
set(gca,'Xtick',[5:10:25]);

subplot(4,2,5);
M1 = nanmean(spk2LFPCH3,1);
SE1 = nanstd(spk2LFPCH3,0,1)./sqrt(size(spk2LFPCH3,1)-1);
M2 = nanmean(spk2LFPCS3,1);
SE2 = nanstd(spk2LFPCS3,0,1)./sqrt(size(spk2LFPCS3,1)-1);
a = [a gca];
jbfill(SFCdat.spk2LFPfreqAx,M2-SE2,M2+SE2,[1 0 0],[1 0 0],1,.5)
hold on;
jbfill(SFCdat.spk2LFPfreqAx,M1-SE1,M1+SE1,[0 0 1],[0 0 1],1,.5)
hold on;
plot(SFCdat.spk2LFPfreqAx,nanmean(spk2LFPCH3,1),'Color',[0 0 1],'Linewidth',2);
plot(SFCdat.spk2LFPfreqAx,nanmean(spk2LFPCS3,1),'Color',[1 0 0],'Linewidth',2);
set(gca,'Xtick',[5:10:25]);

subplot(4,2,7);
M1 = nanmean(spk2LFPCH4,1);
SE1 = nanstd(spk2LFPCH4,0,1)./sqrt(size(spk2LFPCH4,1)-1);
M2 = nanmean(spk2LFPCS4,1);
SE2 = nanstd(spk2LFPCS4,0,1)./sqrt(size(spk2LFPCS4,1)-1);
a = [a gca];
jbfill(SFCdat.spk2LFPfreqAx,M2-SE2,M2+SE2,[1 0 0],[1 0 0],1,.5)
hold on;
jbfill(SFCdat.spk2LFPfreqAx,M1-SE1,M1+SE1,[0 0 1],[0 0 1],1,.5)
hold on;
plot(SFCdat.spk2LFPfreqAx,nanmean(spk2LFPCH4,1),'Color',[0 0 1],'Linewidth',2);
plot(SFCdat.spk2LFPfreqAx,nanmean(spk2LFPCS4,1),'Color',[1 0 0],'Linewidth',2);
set(gca,'Xtick',[5:10:25]);

for curAx = 1:length(a)
    axis(a(curAx),'tight');
    set(a(curAx),'YTick',[min(get(a(curAx),'YTick')) max(get(a(curAx),'YTick'))]);
    box(a(curAx),'off');    
end;
xlabel(a(end),'Frequency [hz]');
ylabel(a(end-2),{'Spike-LFP coupling';'[PLV]'});

set(a,'Layer','top');
set(a,'Fontsize',14);
set(a,'LineWidth',3);
set(gcf,'Color','w');

fx = SFCdat.spk2LFPfreqAx;

dSFC = [];
dSFC(1,:,:) = spk2LFPCH1-spk2LFPCS1;
dSFC(2,:,:) = spk2LFPCH2-spk2LFPCS2;
dSFC(3,:,:) = spk2LFPCH3-spk2LFPCS3;
dSFC(4,:,:) = spk2LFPCH4-spk2LFPCS4;

a = zeros(size(dSFC,1),1);

plotIx = 2;
for curCond = 1:4
    
    x = squeeze(dSFC(curCond,:,:));
    pval = zeros(1,length(fx));
    h = zeros(1,length(fx));
    for curFreq = 1:length( fx )
        
        h(curFreq) = lillietest(squeeze(nanmean(x(:,curFreq),2)));
        %if h(curFreq) ==1
            [pval(curFreq)] = signrank( squeeze(nanmean(x(:,curFreq),2)));
        %else
        %    [~,pval(curFreq)] = ttest( squeeze(nanmean(x(:,curFreq),2)));
        %end;
    end;
    
    chck = pval<.001/(length(fx))*4;
    fLim= [min(fx(find(chck))) max(fx(find(chck)))];
    
    any(chck)
    
    subplot(4,2,plotIx);
    M = nanmean(x,1);
    SE = nanstd(x,0,1)./sqrt(size(x,1)-1);
    a(curCond) = gca;
    hold on;
    if any(chck)
        h = area([fLim],[max(max([M+SE])) max(max([M+SE]))],min(min([M-SE])),'ShowBaseLine','off');
        set(h,'FaceColor',[.85 .85 .85],'EdgeColor',[.85 .85 .85]);
    end;
    hold on;
    plot([min(SFCdat.spk2LFPfreqAx) max(SFCdat.spk2LFPfreqAx)],[0 0],'k--','LineWidth',2);
    jbfill(SFCdat.spk2LFPfreqAx,M-SE,M+SE,[0 0 1],[0 0 1],1,.5)
    hold on;
    plot(SFCdat.spk2LFPfreqAx,M,'Color',[0 0 1],'Linewidth',2);
    set(gca,'Xtick',[5:10:25]);
    plotIx = plotIx+2;
end;

for curAx = 1:length(a)
    axis(a(curAx),'tight');
    set(a(curAx),'YTick',[min(get(a(curAx),'YTick')) max(get(a(curAx),'YTick'))]);
    box(a(curAx),'off');    
end;
xlabel(a(end),'Frequency [hz]');
%ylabel(a(end-2),{'Spike-LFP coupling';'[PLV]'});

%set(a,'YLim',[-0.03 0.07]);
set(a,'Layer','top');
set(a,'Fontsize',14);
set(a,'LineWidth',3);
set(gcf,'Color','w');