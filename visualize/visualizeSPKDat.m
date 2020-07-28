%% Visualize Firing Rate Effects and Extract Average

%%
winlin = 1;
exportfigs=1;
if winlin == 2
    [p2dat] = '/media/hanslmas/rds-share/resultsAUG2019/SPKs_filt/';
    [p2rds] = '/media/hanslmas/rds-share/iEEG_DATA/MICRO/';
    slash = '/';
    addpath(genpath('/media/fred4simon/associativeLearningSingleUnit'));
else
    [p2dat] = 'V:\resultsAUG2019\SPKs_filt\';
    [p2rds] = 'V:\iEEG_DATA/MICRO\';
    slash = '\';
    addpath(genpath('V:\fred4simon\associativeLearningSingleUnit'));
    addpath(genpath('V:\fred4simon\tbx\cbrewer\cbrewer'));
end

[spkFiles] = dir([p2dat, '*_spkParams_Sorting_scr.mat' ]);

%%
chnCnt = 0;
for curFile = 1:length( spkFiles )
    load( [p2dat, spkFiles(curFile).name] );
    if ~isempty( spkDat.SPKseg)
        for k=1:length(spkDat.SPKseg);
            chnCnt = chnCnt + 1;
            
        end
    end;
end;



%%
[ pool_All ] = zeros(chnCnt,6001).*nan;
[ pool_H ] = zeros(chnCnt,6001).*nan;
[ pool_M ] = zeros(chnCnt,6001).*nan;
bl = [-1000 -125];
tsel=[-1000 5000];
tix=find(spkDat.dt2>=tsel(1) & spkDat.dt2 <= tsel(2));
tAx=spkDat.dt3;
tbl=find(tAx>=bl(1) & tAx <= bl(2));
%find(spkDat.dt2 >= tsel(1) & spkDat.dt2 <= tsel(2)) 
%% Load individ Units and summarize
Wind=[1 2000;2001 3000;3001 5000];
chnCnt = 0;
for curFile = 1:length( spkFiles )
    fprintf( spkFiles(curFile).name );
    load( [p2dat, spkFiles(curFile).name] );
    ix = regexp(spkFiles(curFile).name,'_');
    pID = spkFiles(curFile).name( 1:ix(1)-1 );
    expMode = spkFiles(curFile).name( ix(1)+1:ix(2)-1 );
    sesh = spkFiles(curFile).name( ix(2)+1:ix(4)-1 );
    for k=1:length(spkDat.chanLabSPK)
        chnCnt = chnCnt + 1;
        frM=spkDat.frM(k,tix);
        frHit=spkDat.frHit(k,:);
        frMiss=spkDat.frMiss(k,:);
        sdall=std(frM(1,tbl))+0.1;
        sdhit=std(frHit(1,tbl))+0.1;
        sdmiss=std(frMiss(1,tbl))+0.1;
        pool_All(chnCnt,:)=spkDat.frM(k,tix)./sdall;
        pool_H(chnCnt,:)=spkDat.frHit(k,:)./sdhit;
        pool_M(chnCnt,:)=spkDat.frMiss(k,:)./sdmiss;
        pool_id{chnCnt,1}=pID;% Patient ID
        pool_id{chnCnt,2}=expMode; % Experiment
        pool_id{chnCnt,3}=sesh; % Session
        pool_id{chnCnt,4}=[spkDat.chanLabSPK{k} 'C' num2str(spkDat.clusID(k))]; % Channel Label
        Fano(chnCnt,1)=spkDat.fano(k,1);
        tmpspksH=find(spkDat.spkTs{1,k}(tix,spkDat.hitIdx));
        tmpspksM=find(spkDat.spkTs{1,k}(tix,spkDat.missIdx));
        % Check for minimum of 50 spks per condition
        minSpks(chnCnt,1)=min([length(tmpspksH) length(tmpspksM)]);
    end
    fprintf('\n');
end
%Cueidx=find(CueAss==1); 
%Assidx=find(CueAss==2);
%pool_Cue  = pool_All(Cueidx,:);
%pool_Ass  = pool_All(Assidx,:);

% kick out Nans
nanidx=find(~isnan(pool_All(:,1)));
pool_All=pool_All(nanidx,:);
pool_H=pool_H(nanidx,:);
pool_M=pool_M(nanidx,:);
pool_id=pool_id(nanidx,:);

% kick out Neurons with less than 50 spikes in any
idx=find(minSpks > 50);
pool_All=pool_All(idx,:);
pool_H=pool_H(idx,:);
pool_M=pool_M(idx,:);
pool_id=pool_id(idx,:);

%%

figure;
x=[0 0;2 2;3 3].*1000;
y=[-0.4 0.8];
for ln=1:3
    plot(x(ln,:),y,'linewidth',2,'color',[0.4 0.4 0.4], 'linestyle', '--');
    hold on
end
M1 = mean(pool_H,1);
SE1 = std(pool_H,0,1)./sqrt(size(pool_H,1)-1);
M2 = mean(pool_M,1);
SE2 = std(pool_M,0,1)./sqrt(size(pool_M,1)-1);
jbfill(tAx,M1-SE1,M1+SE1,[.0 0 1],[0.0 0 1],1,.45);
hold on;
jbfill(tAx,M2-SE2,M2+SE2,[1 0 0],[1 0 0],1,.45);
hold on
plot(tAx,M1,'Color',[0 0.0 1.0],'LineWidth',2);
hold on
plot(tAx,M2,'Color',[1.0 0.0 0],'LineWidth',2);
xlim([-500 5000]);
ax=gca;
ax.YTick = [-0.4 0 0.8];
ax.YLim = [-0.4 0.8];
ax.YGrid = 'on';
ax.YTickLabel = {'-0.4', '0', '0.8'};
ax.XTick = [0 1000 2000 3000 4000 5000];
ax.XTickLabel = {'0', '1', '2', '3', '4', '5'};
clear ax;
ylim=[-0.4 0.8];
% Stats Pat

Wind=[1 2000;2001 3000;3001 5000];
pId=unique(pool_id(:,1));
sId=unique(pool_id(:,3));
for n=1:length(pId)
    idx=find(ismember(pool_id(:,1),pId(n,1)));
    %Pool_Pat_H(n,:)=mean(pool_H(idx,:),1);
    %Pool_Pat_M(n,:)=mean(pool_M(idx,:),1);
    for w=1:length(Wind)
        tmpidx=find(tAx>=Wind(w,1) & tAx<=Wind(w,2));
        pool_Pat_H(n,w)=mean(mean(pool_H(idx,tmpidx),1),2);
        pool_Pat_M(n,w)=mean(mean(pool_M(idx,tmpidx),1),2);
        
    end
end

% Stats Sessions

Wind=[1 2000;2001 3000;3001 5000];
for n=1:length(sId)
    idx=find(ismember(pool_id(:,3),sId(n,1)));
    %Pool_Pat_H(n,:)=mean(pool_H(idx,:),1);
    %Pool_Pat_M(n,:)=mean(pool_M(idx,:),1);
    pool_sesid(n,:)=pool_id(idx(1),:);
    for w=1:length(Wind)
        tmpidx=find(tAx>=Wind(w,1) & tAx<=Wind(w,2));
        pool_Ses_H(n,w)=mean(mean(pool_H(idx,tmpidx),1),2);
        pool_Ses_M(n,w)=mean(mean(pool_M(idx,tmpidx),1),2);
    end
end

for w=1:length(Wind)
    tmpidx=find(tAx>=Wind(w,1) & tAx<=Wind(w,2));
    pool_Clus_H(:,w)=mean(pool_H(:,tmpidx),2);
    pool_Clus_M(:,w)=mean(pool_M(:,tmpidx),2);
end


for k=1:size(pool_Ses_H,2)
    [~,Ppval(k)]=ttest([pool_Pat_H(:,k)-pool_Pat_M(:,k)]);
    [~,Spval(k)]=ttest([pool_Ses_H(:,k)-pool_Ses_M(:,k)]);
    [~,Npval(k)]=ttest([pool_Clus_H(:,k)-pool_Clus_M(:,k)]);
end

figure;
SME=pool_Clus_H-pool_Clus_M;
vert_rain_plot(SME,[0.2 0.1 1],50,150);
ax=gca;
ax.YLim=[-4 6];
box off;
clear ax
figure;
SME=pool_Ses_H-pool_Ses_M;
vert_rain_plot(SME,[0.2 0.1 1],50,150);
figure;
SME=pool_Pat_H-pool_Pat_M;
vert_rain_plot(SME,[0.2 0.1 1],50,150);

%% find neurons with highest SME
SME=pool_Clus_H(:,3)-pool_Clus_M(:,3);
neuronID=1:length(SME);
nID_ranked=sortrows([neuronID' SME],2);
sel2plt=nID_ranked(:,1);
for k=1:length(sel2plt)
    actplt=sel2plt(k)
    figN=[p2dat, pool_id{actplt,1}, '_', pool_id{actplt,2}, '_', pool_id{actplt,3}, '_', ...
        pool_id{actplt,4}(1:end-2),'_', pool_id{actplt,4}(end-1:end), '.fig'];
    open(figN);
end


%% Save data for SPSS
%save('SME_Ses_FR.mat','SME_ses','pool_Ses_H','pool_Ses_M');
% save('SME_Clus_FR.mat','SME','pool_Clus_H','pool_Clus_M');


