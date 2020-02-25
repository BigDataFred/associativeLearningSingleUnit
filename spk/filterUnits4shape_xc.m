%% Visualize Firing Rate Effects and Extract Average
%
% set paths
addpath(genpath('/media/hanslmas/rds-share/fred4simon/associativeLearningSingleUnit'));
%%
winlin = 2;
exportfigs=1;
if winlin == 2
    [p2datin] = '/media/hanslmas/rds-share/resultsAUG2019/';
    [p2datout] = '/media/hanslmas/rds-share/resultsAUG2019/SPKs_filt/';
    [p2rds] = '/media/hanslmas/rds-share/iEEG_DATA/MICRO/';
    slash = '/';
else
    [p2datin] = 'V:\resultsAUG2019\';
    [p2datout] = 'V:\resultsAUG2019\SPKs_filt\';
    [p2rds] = 'V:\iEEG_DATA/MICRO\';
    slash = '\';
end

[spkFiles] = dir([p2datin, '*_spkParams_Sorting.mat' ]);


%
bl = [-1000 -125];
tsel=[-1000 5000];
%find(spkDat.dt2 >= tsel(1) & spkDat.dt2 <= tsel(2)) 
%% Screen each Unit for Signal or Noise and save figure and data for later processing

for curFile = 20:length( spkFiles )
    fprintf( spkFiles(curFile).name );
    [spkDat] = load( [p2datin, spkFiles(curFile).name] );
    tix=find(spkDat.dt2>=tsel(1) & spkDat.dt2 <= tsel(2));
    tAx=spkDat.dt2(tix);
    tbl=find(tAx>=bl(1) & tAx <= bl(2));
    ix = regexp(spkFiles(curFile).name,'_');
    pID = spkFiles(curFile).name( 1:ix(1)-1 );
    expMode = spkFiles(curFile).name( ix(1)+1:ix(2)-1 );
    sesh = spkFiles(curFile).name( ix(2)+1:ix(4)-1 );
    if ~isempty(spkDat.spkTs)
        [ lfpDatfN ] = dir( [ p2rds, pID, slash,expMode, slash,sesh, slash, '*_lfpDataStimLockedSegmenteddownsampled.mat' ] );
        [ lfpDat ] = load( [ p2rds, pID, slash,expMode, slash,sesh, slash, lfpDatfN.name ] ); % load the LFP data
        [ trlPool, hitIdx, missIdx, trlENC ] = organizeTrlIdxEM( lfpDat );
        % here we assign cluster IDs to distinguish between clusters on the
        % same Microwire
        chanlabs=unique(spkDat.chanLabSPK);
        clusID=[];
        for nlabs=1:length(chanlabs)
            nclus=length(find(ismember(spkDat.chanLabSPK,chanlabs(nlabs,1))));
            for clus=1:nclus
                clusID=[clusID;clus];
            end
        end
        if ~isempty(missIdx) && ~isempty(hitIdx)
            for k=1:length(spkDat.chanLabSPK)
                tmpfr=spkDat.frM(k,tix);
                tmpbl=mean(tmpfr(tbl));
                tmpspks=spkDat.spkTs{1,k}';
                tmpwavf=spkDat.wavf{k,1};
                tmpxc=spkDat.xc(k,:);
                xct=spkDat.dt;
                %tmpspks=zeros(size(tmpspks,1),size(tmpspks,2));% test rasterplot function
                %tmpspks(hitIdx,7000)=1;
                %tmpspks(missIdx,8000)=1;
                for trl=1:length(trlPool)
                    tmpfrs(trl,:) = conv(tmpspks(trl,2:end-1),gausswin(251),'same')./0.251;% compute smoothed firing rate using gaussian kernel
                    tspks = find(tmpspks(trl,:)==1);
                    ISI = tspks(2:end)-tspks(1:end-1);
                    fano(trl,1)=std(ISI)/mean(ISI);% fano factor (1 = Poisson; >1 = bursty);
                end
                
                tmpfrs2=tmpfrs(:,tix);
                frbl=mean(tmpfrs2(:,tbl),2);
                stdall=std(frbl);
                frall=mean(tmpfrs2,1)-mean(frbl,1);
                                
                tmpfrH=mean(tmpfrs2(hitIdx,:),1);
                tmpblH=mean(tmpfrs2(hitIdx,tbl),2);
                
                tmpfrM=mean(tmpfrs2(missIdx,:),1);
                tmpblM=mean(tmpfrs2(missIdx,tbl),2);
                
                frH=tmpfrH-mean(tmpblH,1);
                frM=tmpfrM-mean(tmpblM,1);
                
                tmpfrs=[];
                info=[pID, '-', expMode, '-', sesh, '-', spkDat.chanLabSPK{k}, 'FanoF:', num2str(nanmean(fano))];
                plotSMEspks(tAx,frH,frM,tmpspks(hitIdx,tix),tmpspks(missIdx,tix),info,...
                    tmpwavf,tmpxc,xct);
                UorN(k,1)=input('Unit or noise? Type 1 for unit and 0 for noise');
                if UorN(k) == 1
                    figname=[p2datout,pID, '_',expMode, '_',sesh,'_',spkDat.chanLabSPK{k,1},'_C', num2str(clusID(k,1)),'.fig'];
                    saveas(gcf,figname);
                    close all
                else
                    figname=[p2datout,'Noise',slash, pID, '_',expMode, '_',sesh,'_',spkDat.chanLabSPK{k,1},'_C',num2str(clusID(k,1)),'.fig'];
                    saveas(gcf,figname);
                    close all
                end
                spkDat.frHit(k,:)=frH;
                spkDat.hitIdx=hitIdx;
                spkDat.missIdx=missIdx;
                spkDat.frMiss(k,:)=frM;
                spkDat.fano(k,1)=nanmean(fano,1);
                spkDat.clusID(k,1)=clusID(k,1);
                spkDat.dt3=tAx;
            end
            keepi=find(UorN);
            spkDat.SPKseg=spkDat.SPKseg(1,keepi);
            spkDat.spkTs=spkDat.spkTs(1,keepi);
            spkDat.frM=spkDat.frM(keepi,:);
            spkDat.frSD=spkDat.frSD(keepi,:);
            spkDat.frHit=spkDat.frHit(keepi,:);
            spkDat.frMiss=spkDat.frMiss(keepi,:);
            spkDat.xc=spkDat.xc(keepi,:);
            spkDat.chanLabSPK=spkDat.chanLabSPK(keepi,1);
            spkDat.fano=spkDat.fano(keepi,1);
            spkDat.wavf=spkDat.wavf(keepi,1);
            spkDat.clusID=spkDat.clusID(keepi,1);
            fNout=[p2datout,spkFiles(curFile).name(1:end-4), '_scr.mat'];
            save(fNout,'spkDat');
            clear UorN tmpfrs
        end
    end
    fprintf('\n');
end