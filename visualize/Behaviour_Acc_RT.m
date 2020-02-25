%% Extract Behavioural Data from Task
% Accuracy 
winlin=1;% Windows or Linux Machine
[ pId ] = {'P02','P04','P05','P07','P08','P09','P03ERL','P22AMS','P23AMS'};%
[ expMode ] = {'fVSpEM','cnEM'};
if winlin ==1
    [ rdsPath ] = 'V:\iEEG_DATA\MICRO\';
else
    [ rdsPath ] = '/media/hanslmas/rds-share/iEEG_DATA/MICRO/';
end
nses=0;

varNames = {'PatientID','EMode','Session','NHits', 'Nmiss1', 'Nmiss2', 'RThits', 'RTmiss1', 'RTmiss1'};
%varTypes = {'string','string','double', 'int16', 'int16', 'double', 'double'};
%dat = table('Size',[1 length(varNames)],'VariableTypes',varTypes,'VariableNames',varNames)
dat=cell(1,length(varNames));
for n=1:length(pId)
    for em=1:length(expMode)
        p2d=[rdsPath, pId{1,n}, '/',expMode{1,em}]; 
        sesh=dir(p2d);
        if ~isempty(sesh)
            for se=3:length(sesh)
                fID=dir([p2d,'/', sesh(se).name, '/', pId{n}, '_', sesh(se).name, '_LogFile_EMtask_LogDat.mat']);
                if ~isempty(fID)
                    nses=nses+1;
                    cd([p2d,'/', sesh(se).name, '/']);
                    load(fID.name);
                    dat(nses,1)=pId(n);
                    dat(nses,2)=expMode(em);
                    dat(nses,3)={sesh(se).name};
                    dat{nses,4}=numel(logDat.ix{1,4})/length(logDat.LogDat1.log(:,1));
                    dat{nses,5}=numel(logDat.ix{1,5})/length(logDat.LogDat1.log(:,1));
                    dat{nses,6}=numel(logDat.ix{1,6})/length(logDat.LogDat1.log(:,1));
                    %dat{nses,7}=(dat{nses,4}/length(logDat.LogDat1.log(:,1)))*100;
                    dat{nses,7}=median(logDat.RTs(logDat.ix{1,4}));
                    dat{nses,8}=median(logDat.RTs(logDat.ix{1,5}));
                    dat{nses,9}=median(logDat.RTs(logDat.ix{1,6}));
                    
                    
                end
            end
        end
    end
end

%% Summarize by Patient and Session
cmap=cbrewer('qual','Set1',9);

figure;
x = linspace(0.5,1.5);
y = linspace(16.66,16.66);
plot(x,y,'linewidth',1,'color',[0.5 0.5 0.5],'linestyle', '--');
hold on

for n=1:length(pId)
    pidx=find(ismember(dat(:,1),pId(n)));
    tmpacc=cellfun(@double,dat(pidx,4)).*100;
    tmpm1=cellfun(@double,dat(pidx,5));
    tmpm2=cellfun(@double,dat(pidx,6));
    tmpm=(tmpm1+tmpm2).*100;
    jitt=(rand(1,length(tmpm))-0.5).*0.2;
    scatter(1+jitt,tmpacc,70,cmap(n,:),'filled');
    %scatter(1+jitt,tmpacc,'o-','Color', cmap(n,:),'MarkerSize',10,'MarkerFaceColor',cmap(n,:),'MarkerEdgeColor', [1 1 1]);
    hold on
    scatter(2+jitt,tmpm,70,cmap(n,:),'filled');
    hold on
    Acc_P(n,1)=mean(tmpacc);
    Miss_P(n,1)=mean(tmpm);
    Miss1_P(n,1)=mean(tmpm1)*100;
    Miss2_P(n,1)=mean(tmpm2)*100;
end
xlim([0.5 2.5]);

% Plot Box around mean
Mns(1,1)=mean(Acc_P);
Mns(2,1)=mean(Miss_P);
STE=std(Acc_P)./sqrt(8);% N.B. STE is same for Hits and Misses since Misses is 100-Hits;
xb=[0.8 0.8 1.2 1.2];
yb=[Mns(1,1)-STE Mns(1,1)+STE Mns(1,1)+STE Mns(1,1)-STE];
fill(xb,yb,[0 0 1],'FaceAlpha',.3,'EdgeAlpha',.3);
hold on

xb=[0.8 0.8 1.2 1.2]+1;
yb=[Mns(2,1)-STE Mns(2,1)+STE Mns(2,1)+STE Mns(2,1)-STE];
fill(xb,yb,[1 0 0],'FaceAlpha',.3,'EdgeAlpha',.3);
hold on

% Plot solid Lines for means
x = linspace(0.8,1.2);
y = linspace(mean(Acc_P),mean(Acc_P));

plot(x,y,'linewidth',1,'color',[0 0 1],'linestyle', '-','Color', 'b','LineWidth',3);
hold on

y = linspace(mean(Miss_P),mean(Miss_P));

plot(x+1,y,'linewidth',1,'color',[0 0 1],'linestyle', '-','Color', 'r','LineWidth',3);
hold on

%% Plot RTs

% Accuracy 
winlin=1;% Windows or Linux Machine
[ pId ] = {'P02','P04','P05','P07','P08','P09','P03ERL','P22AMS','P23AMS'};%
[ expMode ] = {'fVSpEM','cnEM'};
if winlin ==1
    [ rdsPath ] = 'V:\resultsAUG2019\RTs\';
else
    [ rdsPath ] = '/media/hanslmas/rds-share/resultsAUG2019/RTs/';
end
nses=0;

fN=dir([rdsPath,'P*_RT.mat']);

varNames = {'PatientID','Session','RTHits', 'RTmiss'};
%varTypes = {'string','string','double', 'int16', 'int16', 'double', 'double'};
%dat = table('Size',[1 length(varNames)],'VariableTypes',varTypes,'VariableNames',varNames)
dat=cell(size(fN,1),length(varNames));
for s=1:length(fN)
    pix = regexp(fN(s).name,'_');
    dat(s,1) = {fN(s).name(1:pix-1)};
    six(1) = regexp(fN(s).name,'\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}');
    six(2) = six(1)+18;
    dat(s,2) = {fN(s).name(six(1):six(2))};
    load([rdsPath, fN(s).name]);
    dat(s,3)={median(RT_hits)};
    dat(s,4)={median(RT_miss)};
end
    
%% Summarize by Patient and Session
cmap=cbrewer('qual','Set1',9);

figure;
x = linspace(0.5,2.5);
y = linspace(5,5);
plot(x,y,'linewidth',1,'color',[0.5 0.5 0.5],'linestyle', '--');
hold on

for n=1:length(pId)
    pidx=find(ismember(dat(:,1),pId(n)));
    tmprth=cellfun(@double,dat(pidx,3))+2;% We add 2 seconds because RTs are locked to Assoc Stimulus
    tmprtm=cellfun(@double,dat(pidx,4))+2;
    jitt=(rand(1,length(tmprth))-0.5).*0.2;
    scatter(1+jitt,tmprth,70,cmap(n,:),'filled');
    %scatter(1+jitt,tmpacc,'o-','Color', cmap(n,:),'MarkerSize',10,'MarkerFaceColor',cmap(n,:),'MarkerEdgeColor', [1 1 1]);
    hold on
    scatter(2+jitt,tmprtm,70,cmap(n,:),'filled');
    hold on
    RTHits_P(n,1)=mean(tmprth);
    RTMiss_P(n,1)=mean(tmprtm);
end

xlim([0.5 2.5]);
ylim([0 45]);
% Plot Box around mean
Mns(1,1)=nanmean(RTHits_P);
Mns(2,1)=nanmean(RTMiss_P);
STE(1,1)=std(RTHits_P)./sqrt(8);
STE(2,1)=nanstd(RTMiss_P)./sqrt(8);

xb=[0.8 0.8 1.2 1.2];
yb=[Mns(1,1)-STE(1,1) Mns(1,1)+STE(1,1) Mns(1,1)+STE(1,1) Mns(1,1)-STE(1,1)];
fill(xb,yb,[0 0 1],'FaceAlpha',.3,'EdgeAlpha',.3);
hold on

xb=[0.8 0.8 1.2 1.2]+1;
yb=[Mns(2,1)-STE(2,1) Mns(2,1)+STE(2,1) Mns(2,1)+STE(2,1) Mns(2,1)-STE(2,1)];
fill(xb,yb,[1 0 0],'FaceAlpha',.3,'EdgeAlpha',.3);
hold on

% Plot solid Lines for means
x = linspace(0.8,1.2);
y = linspace(mean(RTHits_P),mean(RTHits_P));

plot(x,y,'linewidth',1,'color',[0 0 1],'linestyle', '-','Color', 'b','LineWidth',3);
hold on

y = linspace(nanmean(RTMiss_P),nanmean(RTMiss_P));

plot(x+1,y,'linewidth',1,'color',[0 0 1],'linestyle', '-','Color', 'r','LineWidth',3);
hold on
box off
[H,P,ci,stats]=ttest(RTHits_P - RTMiss_P);
