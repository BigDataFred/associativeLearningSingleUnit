%%
addpath('~/prj/Bham/code/mcode/utils/');

%% path to files
[ p2d ] = '/home/rouxf/resultsSpikeFieldJun2019/';

%% patient IDs
[ pId ] = {'P02','P04','P05','P07','P08','P09','P03ERL','P22AMS','P23AMS'};%

%% experiment Mode
[ expMode ] = {'fVSpEM','cnEM'};%{'fVSpEM'};%%

%% initialize variables
[ mwLab   ]     = {}; [poolPowLF]     = {}; [poolPowHF]     = {}; [ poolITC ]     = {}; [ poolPowH ]    = {}; [ poolPowM ]    = {};
[ poolFRM ]     = []; [ poolHitFRM ]  = []; [ poolMssiFRM ] = []; [ chanTotpPat ] = zeros(1,length(pId));

%%
cntLFP = 0;
cntSPK = 0;
cntMWtot = 0;
for curPat = 1:length(pId)
    
    pId{curPat}
    
    for curExp = 1:length( expMode )
        
        fN1 = dir([p2d,pId{curPat},'_',expMode{curExp},'_*_spkParams_noSorting.mat'])
        
        for curShes = 1:length( fN1 )
            
            if ~isempty(fN1)
                
                ix = regexp(fN1(curShes).name,'\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}');
                sesh = fN1(curShes).name(ix:ix+18)
                
                fN2 = dir([p2d,pId{curPat},'_',expMode{curExp},'_',sesh,'_spectrogramTrialSeg_lowFreq.mat']);
                fN3 = dir([p2d,pId{curPat},'_',expMode{curExp},'_',sesh,'_spectrogramTrialSegITC_lowFreq.mat']);
                fN4 = dir([p2d,pId{curPat},'_',expMode{curExp},'_',sesh,'_spectrogramTrialSeg_lowFreqHits.mat']);
                fN5 = dir([p2d,pId{curPat},'_',expMode{curExp},'_',sesh,'_spectrogramTrialSeg_lowFreqMisses.mat']);
                fN6 = dir([p2d,pId{curPat},'_',expMode{curExp},'_',sesh,'_spectrogramTrialSeg_highFreq.mat']);
                fN7 = dir([p2d,pId{curPat},'_',expMode{curExp},'_',sesh,'_spectrogramTrialSeg_highFreqHits.mat']);
                fN8 = dir([p2d,pId{curPat},'_',expMode{curExp},'_',sesh,'_spectrogramTrialSeg_highFreqMisses.mat']);
                
                %%                               
                if ~isempty(fN2)
                    
                    [spkDat] = load([p2d,fN1(curShes).name]);
                    
                    [lfPowDat] = load([p2d,fN2.name]);   
                    [hfPowDat] = load([p2d,fN6.name]);
                    
                    [itcDat] = load([p2d,fN3.name]);                    
                    [hitDat] = load([p2d,fN4.name]);
                    [missDat] = load([p2d,fN5.name]);
                    
                    [chanLabLFP] = unique(lfPowDat.chanLabLFP);
                    
                    %%
                    selIx1 = [];
                    cnt1 = 0;
                    for curChan = 1:length( itcDat.itcH )
                        if ~isempty(itcDat.itcH{curChan})
                            cnt1 = cnt1+1;
                            tmp = itcDat.itcH{curChan};
                            tIx = [find(itcDat.itcTime>=0.25 & itcDat.itcTime <=0.75),find(itcDat.itcTime>=2.25 & itcDat.itcTime <=2.75)];
                            tmp = tmp( itcDat.itcFreq>=4 & itcDat.itcFreq<=12, tIx );                            

                            if ( max(max(tmp)) > 0.45  )
                                selIx1 = [selIx1 curChan];
                            end;
                        end;
                    end;
                    selIx1 = 1:length( itcDat.itc );
                    
                    cntMWtot = cntMWtot+cnt1;
                    chanTotpPat(curPat) = chanTotpPat(curPat)+length(itcDat.itcH);
                    
                    %%
                    chanLabSPK = unique( spkDat.chanLabSPK );
                    p2spkDat = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pId{curPat},'/',expMode{curExp},'/',sesh,'/'];
                    spkFN = dir([p2spkDat,[pId{curPat},'_',expMode{curExp},'_',sesh,'_spkDataStimLockedSegmented.mat']]);
                    spkFN = spkFN.name;
                     
                    [ spkDat2 ] = load([p2spkDat,spkFN]);
                    dsSPKtime = -7:1e-3:7;
                    
                    LFPlab2 = cell(1,length(chanLabLFP));
                    for curChan = 1:length(chanLabLFP)
                        LFPlab2(curChan) = {chanLabLFP{curChan}(1:regexp(chanLabLFP{curChan},'\d{1}')-1)};
                    end;
                    
                    SPKlab2 = cell(1,length(chanLabSPK));
                    for curChan = 1:length(chanLabSPK)
                        SPKlab2(curChan) = {chanLabSPK{curChan}(1:regexp(chanLabSPK{curChan},'\d{1}')-1)};
                    end;
                    
                    %%                    
                    lfpIx = [];
                    for curChan = 1:length( SPKlab2 )
                        lfpIx = [lfpIx find(strcmp(LFPlab2,SPKlab2{curChan}))];
                    end;                 
                    lfpIx = unique(lfpIx);
                    
                    %%
                    spkIx = [];
                    for curChan = 1:length(lfpIx)
                        tmp = find(strcmp(SPKlab2,LFPlab2(lfpIx(curChan))));
                        spkIx = [spkIx tmp];
                    end;
                    spkIx = unique(spkIx);
                    
                    %%                                       
                    for curChan = 1:length(spkIx)
                        cntSPK = cntSPK+1;
                        poolFRM(cntSPK,:) = spkDat.frM(spkIx(curChan),:);
                        
                        spkH = spkDat.spkTs{spkIx(curChan)}(:,spkDat2.hitIdx);
                        spkM = spkDat.spkTs{spkIx(curChan)}(:,spkDat2.missIdx);
                        for curTrl = 1:size(spkH,2)
                            spkH(:,curTrl) = conv(spkH(:,curTrl),gausswin(251),'same')./0.251;
                        end;
                        for curTrl = 1:size(spkM,2)
                            spkM(:,curTrl) = conv(spkM(:,curTrl),gausswin(251),'same')./0.251;
                        end;
                        
                        poolHitFRM(cntSPK,:) = mean(spkH,2);
                        poolMissFRM(cntSPK,:) = mean(spkM,2);
                        
                    end;                                     
                    
                    %%
                    lfPowDat.tx = lfPowDat.tx-7;
                    hfPowDat.tx = hfPowDat.tx-7;
                    for curChan = 1:length( lfpIx )                       
                        lfPowDat.SxxLF{lfpIx(curChan)} = lfPowDat.SxxLF{lfpIx(curChan)}(lfPowDat.tx >=-0.5 & lfPowDat.tx <=5 ,:);                                                
                        hfPowDat.SxxHF{lfpIx(curChan)} = hfPowDat.SxxHF{lfpIx(curChan)}(hfPowDat.tx>=-0.5 & hfPowDat.tx <=5,:);
                    end;
                    lfPowDat.tx = hfPowDat.tx(lfPowDat.tx >=-0.5 & lfPowDat.tx <=5);
                    hfPowDat.tx = hfPowDat.tx(hfPowDat.tx>=-0.5 & hfPowDat.tx <=5);
                    
                    %%
                    hitDat.tx = hitDat.tx-7;
                    missDat.tx = missDat.tx-7;
                    for curChan = 1:length( lfpIx )            
                        try
                        hitDat.SxxLFH{lfpIx(curChan)} = hitDat.SxxLFH{lfpIx(curChan)}(hitDat.tx >=-0.5 & hitDat.tx <=5 ,:);                                                
                        missDat.SxxLFM{lfpIx(curChan)} = missDat.SxxLFM{lfpIx(curChan)}(missDat.tx>=-0.5 & missDat.tx <=5,:);
                        catch %%%
                        end;
                    end;
                    hitDat.tx = hitDat.tx(hitDat.tx >=-0.5 & hitDat.tx <=5);
                    missDat.tx = missDat.tx(missDat.tx>=-0.5 & missDat.tx <=5);
                    
                    %%
                    for curChan = 1:length( lfpIx )
                        xLF = lfPowDat.SxxLF{lfpIx(curChan)};
                        xHF = hfPowDat.SxxHF{lfpIx(curChan)}; 
                        if ~isempty(xLF)
                            cntLFP = cntLFP+1;
                            mwLab(cntLFP) = {[pId{curPat},'_',sesh]};
                                                         
                            %m = (ones(size(xLF,1),1)*mean(xLF(lfPowDat.tx>=-0.5 & lfPowDat.tx<=0,:),1));
                            poolPowLF{cntLFP} = xLF;%20*log10(xLF);%20*log10(xlf./m);%20*log10(xlf);%
                            clear xLF;
                            
                            %m = (ones(size(xHF,1),1)*mean(xHF(hfPowDat.tx>=-0.5 & hfPowDat.tx<=0,:),1));
                            poolPowHF{cntLFP} = xHF;%20*log10(xHF);%20*log10(xhf./m);%20*log10(xhf);%
                            clear xHF;
                            
                            poolITC{cntLFP} = squeeze(itcDat.itcH{lfpIx(curChan)});
                            poolPowH{cntLFP} = hitDat.SxxLFH{lfpIx(curChan)};
                            poolPowM{cntLFP} = missDat.SxxLFM{lfpIx(curChan)};
                        end;
                    end;
                    
                end;
                
            end;
        end;
    end;
end;

%%
nCnt = 0;
selIx = cell(1,length(pId));
for curPat = 1:length( pId )
    tmp = zeros(length(mwLab),1);
    for curChan = 1:length(mwLab);
        tmp(curChan) = ~isempty(regexp(mwLab{curChan},pId{curPat}));
    end;
    selIx{curPat} = find(tmp==1);
    if ~isempty(selIx{curPat})
        nCnt = nCnt+1;
    end;
end;

%%
[ avgPowLF ] = zeros(length(poolPowLF),size(poolPowLF{1},1),size(poolPowLF{1},2));
[ avgPowHF ] = zeros(length(poolPowLF),size(poolPowHF{1},1),size(poolPowHF{1},2));

for curChan = 1:length(poolPowLF)
    
    x = poolPowLF{curChan};
    avgPowLF(curChan,:,:) = x;
    
    x = poolPowHF{curChan};
    avgPowHF(curChan,:,:) = x;
end;


%%
eV = squeeze( max(max(avgPowLF,[],3),[],2) );
[sV,~] = sort(eV);
thr = prctile(eV,[25 75]);
ix = [max(find(eV < thr(1))) min(find(eV > thr(2)))]
exIx = find(eV>5*(thr(2) - thr(1)));
delIx = exIx;

eV2 = squeeze( max(max(avgPowHF,[],3),[],2) );
[sV2,~] = sort(eV2);
thr2 = prctile(eV2,[25 75]);
ix2 = [max(find(eV2 < thr2(1))) min(find(eV2 > thr2(2)))]
exIx2 = find(eV2>5*(thr2(2) - thr2(1)));
delIx2 = exIx2;

figure;
subplot(121);
hold on;
plot(sort(eV),'bs-');
plot([ix(1) ix(1)],[min(eV) max(eV)],'r--');
plot([ix(2) ix(2)],[min(eV) max(eV)],'r--');
plot(find(ismember(sV,eV(exIx))),sort(eV(exIx)),'r*');
subplot(122);
hold on;
plot(sort(eV2),'bs-');
plot([ix2(1) ix2(1)],[min(eV2) max(eV2)],'r--');
plot([ix2(2) ix2(2)],[min(eV2) max(eV2)],'r--');
plot(find(ismember(sV2,eV2(exIx2))),sort(eV2(exIx2)),'r*');

%%
avgPowLF(unique([delIx;delIx2]),:,:) = [];
avgPowHF(unique([delIx;delIx2]),:,:) = [];
poolPowH(unique([delIx;delIx2])) = [];
poolPowM(unique([delIx;delIx2])) = [];

%%
for curChan = 1:size(avgPowLF,1)
    
    x = avgPowLF(curChan,:,:);
    x = (x-min(min(x)))./(max(max(x))-min(min(x)));
    avgPowLF(curChan,:,:) = x;
    
    x = avgPowHF(curChan,:,:);
    x = (x-min(min(x)))./(max(max(x))-min(min(x)));
    avgPowHF(curChan,:,:) = x;
end;

%%
avgPowH = zeros(length(poolPowH),size(poolPowH{1},1),size(poolPowH{1},2));
for curChan = 1:length( poolPowH)
    m1 = min(min([poolPowH{curChan}(:) poolPowM{curChan}(:)]));
    m2 = max(max([poolPowH{curChan}(:) poolPowM{curChan}(:)]));
    x = [poolPowH{curChan}];
    x = (x-m1)./(m2-m1);
    avgPowH(curChan,:,:) = x;
end;

avgPowM = NaN(length(poolPowM),size(poolPowH{1},1),size(poolPowH{1},2));
for curChan = 1:length( poolPowM)
    if ~isempty(poolPowM{curChan})
        m1 = min(min([poolPowH{curChan}(:) poolPowM{curChan}(:)]));
        m2 = max(max([poolPowH{curChan}(:) poolPowM{curChan}(:)]));
        x = [poolPowM{curChan}];
        x = (x-m1)./(m2-m1);
        avgPowM(curChan,:,:) = x;
    end;
end;

avgPowD = NaN(length(poolPowH),size(poolPowH{1},1),size(poolPowH{1},2));
for curChan = 1:length( poolPowH)
    if ~isempty(poolPowH{curChan})
        if (~isempty(poolPowH{curChan})) && (~isempty(poolPowM{curChan}))
            m1 = min(min([poolPowH{curChan}(:) poolPowM{curChan}(:)]));
            m2 = max(max([poolPowH{curChan}(:) poolPowM{curChan}(:)]));
            x1 = [poolPowH{curChan}];
            x1 = (x1-m1)./(m2-m1);
            x2 = [poolPowM{curChan}];
            x2 = (x2-m1)./(m2-m1);
            x = x1-x2;
            avgPowD(curChan,:,:) = x;
        end;
    end;
end;

%%
avgITC = zeros(size(poolITC{1}));
for curChan = 1:length( poolITC)
    avgITC = avgITC + poolITC{curChan};
end;
avgITC = avgITC./length(poolITC);

patCnt = {};
for curChan = 1:length( mwLab )
    patCnt{curChan} = mwLab{curChan}(1:regexp(mwLab{curChan},'_')-1);
end;

MWpPatCnt = zeros(1,length( pId ));
for curPat = 1:length( pId )
    MWpPatCnt(curPat) = length(find(strcmp(pId(curPat),patCnt)));
end;

%% figure 3
dsTrlTime = -7.001:1e-3:7.001;

zSC = zeros(size(poolFRM));
uIx1 = [];
uIx2 = [];
tIx = find(dsTrlTime >-1 & dsTrlTime<0);
tIx2 = find(dsTrlTime >=-.5 & dsTrlTime<=5);
for curChan = 1:size(poolFRM,1)
    
    xFR = poolFRM(curChan,:);
    mFR = mean(xFR(tIx));
    sdFR = std(xFR(tIx))+0.1;
    zSC(curChan,:) = (xFR-mFR)./sdFR;
    
    [~,mIx] = max(zSC(curChan,tIx2));
    
    if (dsTrlTime(tIx2(mIx)) > 0) & (dsTrlTime(tIx2(mIx)) <2)
        uIx1 = [uIx1 curChan];
    elseif (dsTrlTime(tIx2(mIx)) > 2)
        uIx2 = [uIx2 curChan];
    end;
    
end;
zSC = zSC(:,tIx2);
dsTrlTime = dsTrlTime(tIx2);

M = mean(zSC,1);
SE = std(zSC,0,1)./sqrt(size(zSC,1)-1);

M1 = mean(zSC(uIx1,:),1);
SE1 = std(zSC(uIx1,:),0,1)./sqrt(length(uIx1)-1);

M2 = mean(zSC(uIx2,:),1);
SE2 = std(zSC(uIx2,:),0,1)./sqrt(length(uIx2)-1);

zSCH = zeros(size(poolHitFRM));
zSCM = zeros(size(poolMissFRM));
for curChan = 1:size(poolFRM,1)
    
    xFRH = poolHitFRM(curChan,:);
    mFRH = mean(xFRH(tIx));
    sdFRH = std(xFRH(tIx))+0.1;
    zSCH(curChan,:) = (xFRH-mFRH)./sdFRH;
    
    xFRM = poolMissFRM(curChan,:);
    mFRM = mean(xFRM(tIx));
    sdFRM = std(xFRM(tIx))+0.1;
    zSCM(curChan,:) = (xFRM-mFRM)./sdFRM;
    
end;
zSCH = zSCH(:,tIx2);
zSCM = zSCM(:,tIx2);

MH = nanmean(zSCH,1);
SEH = nanstd(zSCH,0,1)./sqrt(size(zSCH,1)-1);

MM = nanmean(zSCM,1);
SEM = nanstd(zSCM,0,1)./sqrt(size(zSCM,1)-1);

figure;
hold on;
h = area([0 2],max(M+SE)*ones(1,2),min(M-SE));
set(h,'FaceColor',[.95 .95 .95],'EdgeColor',[.95 .95 .95]);
h = area([2 3],max(M+SE)*ones(1,2),min(M-SE));
set(h,'FaceColor',[.85 .85 .85],'EdgeColor',[.75 .75 .75]);
jbfill(dsTrlTime,M-SE,M+SE,[253 162 2]./255,[253 162 2]./255,1,.1);
hold on;
plot([-.5 5],[0 0],'k--','LineWidth',2);
plot(dsTrlTime,M,'Color',[253 162 2]./255,'LineWidth',2);
axis tight;
set(gca,'Layer','Top');
xlim([-.5 5]);
set(gca,'XTick',[0 2 3]);
set(gca,'XTickLabel',repmat('',[1 length(get(gca,'XTick'))]));
set(gca,'LineWidth',3);
box(gca,'off');
set(gca,'Fontsize',14);
set(gca,'YTick',[min(get(gca,'YTick')) max(get(gca,'YTick'))]);
ylabel('Firing Rate [\sigma]');
xlabel('Time rel. to Cue onset');
set(gcf,'Color','w');


figure;
hold on;
h = area([0 2],max(M1+SE1)*ones(1,2),min(M1-SE1));
set(h,'FaceColor',[.95 .95 .95],'EdgeColor',[.95 .95 .95]);
h = area([2 3],max(M1+SE1)*ones(1,2),min(M1-SE1));
set(h,'FaceColor',[.85 .85 .85],'EdgeColor',[.75 .75 .75]);
jbfill(dsTrlTime,M1-SE1,M1+SE1,[11 102 35]./255,[11 102 35]./255,1,.25);
hold on;
jbfill(dsTrlTime,M2-SE2,M2+SE2,[0 0 0],[0 0 0],1,.25);
hold on;
plot([-.5 5],[0 0],'k--','LineWidth',2);
plot(dsTrlTime,M1,'Color',[11 102 35]./255,'LineWidth',2);
plot(dsTrlTime,M2,'Color',[0 0 0],'LineWidth',2);
axis tight;
set(gca,'Layer','Top');
xlim([-.5 5]);
set(gca,'XTick',[0 2 3]);
set(gca,'XTickLabel',repmat('',[1 length(get(gca,'XTick'))]));
set(gca,'LineWidth',3);
box(gca,'off');
set(gca,'Fontsize',14);
set(gca,'YTick',[min(get(gca,'YTick')) max(get(gca,'YTick'))]);
ylabel('Firing Rate [\sigma]');
xlabel('Time rel. to Cue onset');
set(gcf,'Color','w');

savePath = '/media/rouxf/rds-share/Fred/figures4Cali2018/figure3/';
%saveas(1,[savePath,'fVSpEM_avgFiringRate.fig']);
%saveas(2,[savePath,'fVSpEM_avgFiringRate4UnitsSplitByMax.fig']);

figure;
hold on;
h = area([0 2],max([MH+SEH MM+SEM])*ones(1,2),min([MH-SEH MM-SEM]));
set(h,'FaceColor',[.95 .95 .95],'EdgeColor',[.95 .95 .95]);
h = area([2 3],max([MH+SEH MM+SEM])*ones(1,2),min([MH-SEH MM-SEM]));
set(h,'FaceColor',[.85 .85 .85],'EdgeColor',[.75 .75 .75]);
jbfill(dsTrlTime,MH-SEH,MH+SEH,[0 0 1],[0 0 1],1,.25);
hold on;
jbfill(dsTrlTime,MM-SEM,MM+SEM,[1 0 0],[1 0 0],1,.25);
hold on;
plot([-.5 5],[0 0],'k--','LineWidth',2);
plot(dsTrlTime,MH,'Color',[0 0 1],'LineWidth',2);
plot(dsTrlTime,MM,'Color',[1 0 0],'LineWidth',2);
axis tight;
set(gca,'Layer','Top');
xlim([-.5 5]);
set(gca,'XTick',[0 2 3]);
set(gca,'XTickLabel',repmat('',[1 length(get(gca,'XTick'))]));
set(gca,'LineWidth',3);
box(gca,'off');
set(gca,'Fontsize',14);
set(gca,'YTick',[min(get(gca,'YTick')) max(get(gca,'YTick'))]);
ylabel('Firing Rate [\sigma]');
xlabel('Time rel. to Cue onset');
set(gcf,'Color','w');

figure;
X = [nanmean(zSCH(:,dsTrlTime>0 & dsTrlTime<=2),2) nanmean(zSCH(:,dsTrlTime>2 & dsTrlTime<=3),2)  nanmean(zSCH(:,dsTrlTime>3),2);
    nanmean(zSCM(:,dsTrlTime>0 & dsTrlTime<=2),2) nanmean(zSCM(:,dsTrlTime>2 & dsTrlTime<=3),2)  nanmean(zSCM(:,dsTrlTime>3),2)];

MH = nanmean(X(1:size(zSCH,1),:),1);
SEH = nanstd(X(1:size(zSCH,1),:),0,1)./sqrt(size(X(1:size(zSCH,1),:),1)-1);

MM = nanmean(X(size(zSCH,1)+1:size(zSCH,1)+size(zSCM,1),:),1);
SEM = nanstd(X(size(zSCH,1)+1:size(zSCH,1)+size(zSCM,1),:),0,1)./sqrt(size(X(size(zSCH,1)+1:size(zSCH,1)+size(zSCM,1),:),1)-1);

hold on;
errorbar(1:3,MH,SEH,'s-','Color',[0 0 1],'MarkerFaceColor',[0 0 1],'LineWidth',2);
errorbar(1:3,MM,SEM,'s-','Color',[1 0 0],'MarkerFaceColor',[1 0 0],'LineWidth',2);
set(gca,'LineWidth',3);
set(gca,'Layer','Top');
box(gca,'off');
set(gca,'Fontsize',14);
set(gca,'YTick',[min(get(gca,'YTick')) max(get(gca,'YTick'))]);
set(gca,'XTickLabel',{'Cue','Association','Response'});
set(gca,'XTick',[1 2 3]);
ylabel('Average firing rate [\sigma]');
set(gcf,'Color','w');

%saveas(1,[savePath,'fVSpEM_avgFiringRateHitsVSMisses.fig']);
%saveas(2,[savePath,'fVSpEM_avgFiringRateHitsVSMissesErrorBarPlot.fig']);

%% figure2
tx = lfPowDat.tx;
fx = lfPowDat.fx;
tx2 = hfPowDat.tx;
fx2 = hfPowDat.fx;
itcTime = itcDat.itcTime;
itcFreq = itcDat.itcFreq;


figure;
a = gca;
hold on;
imagesc( tx, fx, squeeze(mean(avgPowLF,1))' );
axis xy;colormap jet;
xlim([-.5 5]);%caxis([.05 .225]);

figure;
b = gca;
hold on;
imagesc( tx2, fx2, squeeze(mean(avgPowHF,1))' );
axis xy;colormap jet;
xlim([-.5 5]);%caxis([.05 .225]);
ylabel(b,'Frequency [hz]');
plot(b,[0 0],[min(fx2) max(fx2)],'Color',[0 0 0],'LineWidth',2);
plot(b,[2 2],[min(fx2) max(fx2)],'Color',[0 0 0],'LineWidth',2);
plot(b,[3 3],[min(fx2) max(fx2)],'Color',[0 0 0],'LineWidth',2);
box(b,'on');


figure;
a = [a gca];
hold on;
imagesc(itcDat.itcTime,itcDat.itcFreq,avgITC);
axis xy;colormap jet;
xlim([-.5 5]);caxis([.1 .2]);
for curAx = 1:length(a)
    ylabel(a(curAx),'Frequency [hz]');
    if (curAx ==1)
        plot(a(curAx),[0 0],[min(fx) max(fx)],'Color',[0 0 0],'LineWidth',2);
        plot(a(curAx),[2 2],[min(fx) max(fx)],'Color',[0 0 0],'LineWidth',2);
        plot(a(curAx),[3 3],[min(fx) max(fx)],'Color',[0 0 0],'LineWidth',2);
    else
        xlabel(a(curAx),'Time re. to Cue onset');
        plot(a(curAx),[0 0],[min(itcFreq) max(itcFreq)],'Color',[0 0 0],'LineWidth',2);
        plot(a(curAx),[2 2],[min(itcFreq) max(itcFreq)],'Color',[0 0 0],'LineWidth',2);
        plot(a(curAx),[3 3],[min(itcFreq) max(itcFreq)],'Color',[0 0 0],'LineWidth',2);
    end;
    box(a(curAx),'on');
end;

for curAx = 1:length(a)
    axis(a(curAx),'tight');
    xlim(a(curAx),[-.5 4.9]);
    ylim(a(curAx),[0.5 30]);
end;
for curAx = 1:length(b)
    axis(b(curAx),'tight');
    xlim(b(curAx),[-.5 4.9]);
    ylim(b(curAx),[30 170]);
end;

set([a b],'Layer','Top');
set(a,'YTick',[5 15 25]);
set([a b],'XTick',[0:5]);
set([a b],'XTickLabel',{repmat('',[1 length(get(gca,'XTick'))])});
set([a b],'Fontsize',14);
set([a b],'LineWidth',3);
set([a b],'Color',[.75 .75 .75]);
set(gcf,'Color','w');

savePath = '/media/rouxf/rds-share/Fred/figures4Cali2018/figure2/';

%saveas(1,[savePath,'fVSpEM_avgSpectrogramANDitcForSpikingChannelsPOW.fig']);
%saveas(2,[savePath,'fVSpEM_avgSpectrogramANDitcForSpikingChannelsITC.fig']);

figure;
subplot(441);
cb = colorbar;
caxis([.05 .25]);
colormap jet;
zlab = get(cb,'YLabel');
set(zlab,'String','LFP-power [a.u]');
set(cb,'Fontsize',14);
set(cb,'YTick',[min(get(cb,'YLim')) max(get(cb,'YLim'))]);
axis off;
set(gcf,'Color','w');
%saveas(gcf,[savePath,'fVSpEM_avgSpectrogramANDitcForSpikingChannelsColorbarPow.fig']);

figure;
subplot(441);
cb = colorbar;
caxis([.1 .2]);
colormap jet;
zlab = get(cb,'YLabel');
set(zlab,'String','ITC [a.u]');
set(cb,'Fontsize',14);
set(cb,'YTick',[min(get(cb,'YLim')) max(get(cb,'YLim'))]);
axis off;
set(gcf,'Color','w');
%saveas(gcf,[savePath,'fVSpEM_avgSpectrogramANDitcForSpikingChannelsColorbarITC.fig']);


MWpPatCnt = zeros(1,length(pId));
for curPat = 1:length(pId)
    chck = regexp(mwLab,pId{curPat});
    cnt = 0;
    for curChan = 1:length( chck )
        if ~isempty( chck{curChan} );
            cnt = cnt+1;
        end;
    end;
    MWpPatCnt(curPat) = cnt;
end;

figure;
subplot(223);
a = gca;
h = bar( 1:length(pId), MWpPatCnt./length(mwLab) );
set(gca,'XTick',1:length(pId));
set(gca,'XTickLabel',pId);
set(gca,'XTickLabelRotation',-45);
axis tight;
box off;
set(gca,'YTick',[0 round(max(get(gca,'YLim')).*10)./10-.05]);
ylabel('% Channels');
set(h,'FaceColor','b');
subplot(221);
a = [a gca];
hold on;
h = [];
h(1) = bar(1,sum(MWpPatCnt),.5);
h(2) = bar(2,sum(chanTotpPat),.5);
axis tight;box off;
set(h(1),'FaceColor',[0 0 0]);
set(h(2),'FaceColor',[.75 .75 .75]);
set(gca,'YTick',[sum(MWpPatCnt) sum(chanTotpPat)]);
xlim([0 3]);
set(gca,'XTick',[]);
ylabel('Count');

set(a,'Layer','Top');
set(a,'Fontsize',14);
set(a,'LineWidth',3);
set(gcf,'Color','w');
%saveas(gcf,[savePath,'microWireCountAndSelectedMWratio.fig']);

%%
tw = [-.5 0;0 2;2 3;3 max(tx)];
pnt= 1;
pnt2=2;
figure;
for curTW = 1:size(tw,1)
    
    pval = zeros(1,length(fx));
    for curFreq = 1:length( fx )
        [pval(curFreq)] = signrank(squeeze(nanmean(avgPowM(:,tx>=tw(curTW,1) & tx<tw(curTW,2),curFreq),2)),squeeze(nanmean(avgPowH(:,tx>=tw(curTW,1) & tx <tw(curTW,2),curFreq),2)));
    end;
    chck = pval<1e-3/length(fx)*size(tw,1);
    
    M1 = squeeze(nanmean(nanmean(avgPowH(:,tx>=tw(curTW,1) & tx <tw(curTW,2),:),2),1));
    SE1 = squeeze(nanstd(nanmean(avgPowH(:,tx>=tw(curTW,1) & tx <tw(curTW,2),:),2),0,1))./sqrt(size(avgPowH,1)-1);
    M2 =  squeeze(nanmean(nanmean(avgPowM(:,tx>=tw(curTW,1) & tx <tw(curTW,2),:),2),1));
    SE2 = squeeze(nanstd(nanmean(avgPowM(:,tx>=tw(curTW,1) & tx <tw(curTW,2),:),2),0,1))./sqrt(size(avgPowM,1)-1);
    
    fLim= [min(fx(find(chck))) max(fx(find(chck)))];
    subplot(4,2,pnt);
    hold on;
    %     if any(chck)
    %         h = area([fLim],[max(max([M1+SE1 M2+SE2])) max(max([M1+SE1 M2+SE2]))],min(min([M1-SE1 M2-SE2])));
    %         set(h,'FaceColor',[.85 .85 .85],'EdgeColor',[.85 .85 .85]);
    %     end;
    jbfill(fx,M1'-SE1',M1'+SE1',[0 0 1],[0 0 1],1,.2);
    hold on;
    jbfill(fx,M2'-SE2',M2'+SE2',[1 0 0],[1 0 0],1,.2);
    hold on;
    plot(fx,M1,'Color',[0 0 1],'LineWidth',2);
    plot(fx,M2,'Color',[1 0 0],'LineWidth',2);
    set(gca,'XTick',[5 15 25]);
    set(gca,'XTickLabel',{'','',''});
    if pnt == 3
        ylabel('LFP-power [a.u.]');
    end;
    if pnt == 7
        xlabel('Frequency [hz]');
        set(gca,'XTickLabel',[5 15 25]);
    end;
    axis(gca,'tight');
    
    set(gca,'YTick',[min(get(gca,'YTick')) max(get(gca,'YTick'))]);
    set(gca,'Layer','Top');
    set(gca,'Fontsize',14);
    set(gca,'LineWidth',3);
    pnt = pnt+2;
    
    pval = zeros(1,length(fx));
    h = zeros(1,length(fx));
    for curFreq = 1:length( fx )
        
        h(curFreq) = lillietest(squeeze(nanmean(avgPowD(:,tx>=tw(curTW,1) & tx<tw(curTW,2),curFreq),2)));
        
        [pval(curFreq)] = signrank( squeeze(nanmean(avgPowD(:,tx>=tw(curTW,1) & tx<tw(curTW,2),curFreq),2)));
    end;
    all(h)
    
    chck = pval<.001/length(fx)*size(tw,1);
    
    M = squeeze(nanmean(nanmean(avgPowD(:,tx>=tw(curTW,1) & tx <tw(curTW,2),:),2),1));
    SE = squeeze(nanstd(nanmean(avgPowD(:,tx>=tw(curTW,1) & tx <tw(curTW,2),:),2),0,1))./sqrt(size(avgPowD,1)-1);
    
    fLim= [min(fx(find(chck))) max(fx(find(chck)))];
    subplot(4,2,pnt2);
    hold on;
    if any(chck)
        h = area([fLim],[max(max([M+SE])) max(max([M+SE]))],min(min([M-SE])));
        set(h,'FaceColor',[.85 .85 .85],'EdgeColor',[.85 .85 .85]);
    end;
    hold on;
    plot([min(fx) max(fx)],[0 0],'k--','LineWidth',2);
    jbfill(fx,M'-SE',M'+SE',[0 0 1],[0 0 1],1,.2);
    hold on;
    plot(fx,M,'Color',[0 0 1],'LineWidth',2);
    set(gca,'XTick',[5 15 25]);
    set(gca,'XTickLabel',{'','',''});
    %     if pnt == 4
    %         ylabel('LFP-power [a.u.]');
    %     end;
    if pnt2 == 8
        xlabel('Frequency [hz]');
        set(gca,'XTickLabel',[5 15 25]);
    end;
    axis(gca,'tight');
    
    set(gca,'YTick',[min(get(gca,'YTick')) max(get(gca,'YTick'))]);
    set(gca,'Layer','Top');
    set(gca,'Fontsize',14);
    set(gca,'LineWidth',3);
    pnt2 = pnt2+2;
end;
set(gcf,'Color','w');

%saveas(gcf,[savePath,'fVSpEM_avgSpectrogramSubsequenMemoryEffectLFP-POW.fig']);

%%
 %%
                    %             frC = {};
                    %             cluLab = {};
                    %             mwLab = {};
                    %             xc = {};
                    %             cluCnt = 0;
                    %             for  curMW = 1:length(sortedSpikesSEG )
                    %                 fprintf([num2str(curMW),'/',num2str(length(sortedSpikesSEG ))]);
                    %                 cID = unique(sortedSpikesSEG{curMW}.assignedClusterSeg);
                    %                 for curClus = 1:length(cID)
                    %                     if ( cID(curClus) ~= 0 )
                    %                         cluCnt = cluCnt+1;
                    %                         spkTs2{cluCnt} = zeros(length(dsSPKtime),length(trlENC));
                    %                         cluIx = sortedSpikesSEG{curMW}.assignedClusterSeg == cID(curClus);
                    %                         frC{cluCnt} = zeros(size(spkTs2{cluCnt}));
                    %                         cluLab{cluCnt} = cID(curClus);
                    %                         xcTrl= zeros(length(trlENC),length(dt2));
                    %                         mwLab{cluCnt}  = curMW;
                    %                         for curTrl = 1:length( trlENC )
                    %                             ts = sortedSpikesSEG{curMW}.SpikeTimesSeg(sortedSpikesSEG{curMW}.trl(cluIx) == trlENC( curTrl )).*1e3;
                    %                             [n,~] = hist(ts,dt);
                    %                             spkTs2{cluCnt}(:,curTrl) = n(2:end-1);
                    %                             x = spkTs2{cluCnt}(:,curTrl);
                    %                             frC{cluCnt}(:,curTrl) = conv(x,gausswin(251),'same')./0.251;
                    %                             xcTrl(curTrl,:) = xcorr(n(2:end-1),300);
                    %                         end;
                    %                         xc{cluCnt} = mean(xcTrl,1);
                    %                     end;
                    %                 end;
                    %                 fprintf('\n');
                    %             end;
                    %
                    %             spkSigIx = [];
                    %             for curChan = 1:length(spkTs2)
                    %                 spkCntPost = sum(spkTs2{curChan}(dsSPKtime>2.2 & dsSPKtime < 2.9,:),1);
                    %                 spkCntBase = sum(spkTs2{curChan}(dsSPKtime>-.9 & dsSPKtime <-0.2,:),1);
                    %
                    %                 [~,pval] = ttest(spkCntPost,mean(spkCntBase),0.025,'right');
                    %                 if pval < 5e-2%( median(spkCntPost) > trsh )
                    %                     spkSigIx = [spkSigIx curChan];
                    %                 end;
                    %
                    %             end;
                    
%%
   %             %%
                    %             savePath = '/media/rouxf/rds-share/Fred/figures4Cali2018/figure4/';
                    %             if ~isempty(spkSigIx)
                    %                 for curChan = 1:length( spkSigIx )
                    %                     if ( mean(mean(frC{spkSigIx(curChan)}(dsSPKtime>-.5 & dsSPKtime<5,:),1)) >0.75 )
                    %                         fLab = [pId{curPat},'-fVSpEM-',sesh];
                    %                         fLab(regexp(fLab,'_')) = '_';
                    %                         figure;
                    %                         subplot(4,5,[1 1.5 6 6.5]);
                    %                         hold on;
                    %                         cluIx = sortedSpikesSEG{mwLab{spkSigIx(curChan)}}.oriIx(sortedSpikesSEG{mwLab{spkSigIx(curChan)}}.assignedClusterSeg == cluLab{spkSigIx(curChan)});
                    %                         plot(linspace(0,2,64),sortedSpikesSEG{mwLab{spkSigIx(curChan)}}.wavf(cluIx,:),'r');
                    %                         plot(linspace(0,2,64),mean(sortedSpikesSEG{mwLab{spkSigIx(curChan)}}.wavf(cluIx,:),1),'k');
                    %                         axis off;
                    %                         subplot(4,5,[3 4 5 8 9 10 13 14 15]);
                    %                         a = gca;
                    %                         hold on;
                    %                         for curTrl = 1:size(spkTs2{spkSigIx(curChan)},2)
                    %                             x = spkTs2{spkSigIx(curChan)}(:,curTrl);
                    %                             x = dsSPKtime(x==1);
                    %                             y = curTrl*ones(1,length(x));
                    %                             x = [x;x];
                    %                             y = [y-.5;y+.5];
                    %                             line(x,y,'Color','k');
                    %                         end;
                    %                         xlim([-.5 5]);
                    %                         ylim([0 curTrl+1]);
                    %                         set(gca,'XTick',[]);
                    %                         ylabel('Trial #');
                    %                         set(gca,'YTick',[min(get(gca,'YTick')) max(get(gca,'YTick'))]);
                    %                         text(-.5,curTrl+2,fLab,'FontWeight','bold');
                    %                         subplot(4,5,[11 11.25 16 16.25]);
                    %                         a  = [a gca];
                    %                         xc{spkSigIx(curChan)}(dt2==0)= NaN;
                    %                         h = bar(dt2(dt2>1),xc{spkSigIx(curChan)}(dt2>1));
                    %                         set(h,'FaceColor',[0 0 0]);
                    %                         axis tight;
                    %                         set(gca,'YTick',[0 round(max(get(gca,'YLim'))*10)/10]);
                    %                         ylabel('Coincidences');
                    %                         xlabel('Lag [ms]');
                    %                         box(gca,'off');
                    %                         set(gca,'XTick',[150 300]);
                    %                         subplot(4,5,[18 19 20]);
                    %                         a = [a gca];
                    %                         hold on;
                    %                         plot(dsSPKtime,mean(frC{spkSigIx(curChan)},2),'k','LineWidth',2);
                    %                         plot([0 0],[min(mean(frC{spkSigIx(curChan)},2)) max(mean(frC{spkSigIx(curChan)},2))],'r','LineWidth',2);
                    %                         plot([2 2],[min(mean(frC{spkSigIx(curChan)},2)) max(mean(frC{spkSigIx(curChan)},2))],'r','LineWidth',2);
                    %                         plot([3 3],[min(mean(frC{spkSigIx(curChan)},2)) max(mean(frC{spkSigIx(curChan)},2))],'r','LineWidth',2);
                    %                         axis tight;box off;
                    %                         xlim([-.5 5]);
                    %                         set(gca,'YTick',[min(get(gca,'YTick')) max(get(gca,'YTick'))]);
                    %                         ylabel('Spikes/s^{-1}');
                    %                         xlabel('Time rel. to Cue onset');
                    %                         set(gca,'XTick',[]);
                    %                         set(a,'Fontsize',14);
                    %                         set(a,'LineWidth',3);
                    %                         set(a,'Layer','Top');
                    %                         set(gcf,'Color','w');
                    %                         %saveas(gcf,[savePath,pId{curPat},'_',sesh,'_cluster',num2str(curChan),'_rastergramANDxcorr.fig']);
                    %                     end;
                    %                 end;
                    %             end;