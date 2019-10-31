
%%
addpath(genpath('/home/rouxf/associativeLearningSingleUnit'));

%%
rdsPath1 = '/media/rouxf/rds-share/iEEG_DATA/MICRO/';
rdsPath2 = '/media/rouxf/rds-share/resultsAUG2019/';

%%
pID = {'P02','P04','P05','P07','P09','P03ERL','P22AMS','P22AMS','P23AMS'};
expMode = {'fVSpEM','cnEM'};

%%
cntChan = 0;
poolITC   = zeros(500,100,5501);
poolERP   = zeros(500,5501);
poolPowLF = zeros(500,550,246);
poolPowHF = zeros(500,550,1147);

%%
for curPat = 1:length( pID )
    for curExp = 1:length( expMode )
        
        [ tmp ] = dir( [ rdsPath1, pID{curPat},'/', expMode{curExp}, '/' ] );
        
        if ~isempty(tmp)
            
            [ sesh ] = extractSeshLabels( tmp );
            clear tmp;
            
            for curSesh = 1:length( sesh )
                
                [ spk2lfpFn ] = dir( [ rdsPath2, pID{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spk2LFPCoupling_plv_Sorting_alpha0.05_nRand400_lowFreq.mat' ] );
                
                if ~isempty(spk2lfpFn)
                    [ spk2lfpDat ] = load( [ rdsPath2,spk2lfpFn.name ] );                                        
                    
                    if ~isempty(spk2lfpDat.sigIxLFP)                                                
                        
                        [ itcN ] = dir( [ rdsPath2, pID{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spectrogramTrialSegITC_lowFreq.mat' ] );
                        [ itcDat ] = load([ rdsPath2, itcN.name ]);
                        
                        [ erpN] = dir( [ rdsPath2, pID{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_ERPTrialSeg_HitsANDmisses.mat'] );
                        [ erpDat] = load( [ rdsPath2, erpN.name ] );
                        
                        [ selIxITC ] =[];
                        for curChan= 1:length( spk2lfpDat.chanLabLFPsig )
                            [ selIxITC(curChan) ] = find( strcmp( itcDat.chanLabLFP,spk2lfpDat.chanLabLFPsig(curChan) ) );
                        end;
                        
                        [ lfN1 ] = dir([ rdsPath2, pID{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spectrogramTrialSeg_lowFreqHits.mat' ]);
                        [ lfN2 ] = dir([ rdsPath2, pID{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spectrogramTrialSeg_lowFreqMisses.mat' ]);
                         
                        [ hfN1 ] = dir([ rdsPath2, pID{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spectrogramTrialSeg_highFreqHits.mat' ]);
                        [ hfN2 ] = dir([ rdsPath2, pID{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spectrogramTrialSeg_highFreqMisses.mat' ]);
                        
                        [tfrDat1h] = load([rdsPath2,lfN1.name]);
                        [tfrDat1m] = load([rdsPath2,lfN2.name]);
                        
                        [tfrDat2h] = load([rdsPath2,hfN1.name]);
                        [tfrDat2m] = load([rdsPath2,hfN2.name]);
                        
                        [ selIxPOW ] =[];
                        for curChan= 1:length( spk2lfpDat.chanLabLFPsig )
                            [ selIxPOW(curChan) ] = find( strcmp( tfrDat1h.chanLabLFP,spk2lfpDat.chanLabLFPsig(curChan) ) );
                        end;
                        
                        if (selIxITC ~= selIxPOW)
                            error('channel indexes must be equal');
                        end;
                        
                        for curChan = 1:length( selIxITC )
                            
                            if ~isempty( tfrDat1h.SxxLFH{selIxITC(curChan)} ) && ~isempty( tfrDat1m.SxxLFM{selIxITC(curChan)} )
                                cntChan = cntChan +1;
                                
                                poolITC(cntChan,:,:) = itcDat.itc{ selIxITC(curChan) }(:,itcDat.itcTime>=-.5 & itcDat.itcTime<=5);
                                
                                %poolERP(cntCHan,:) = erpDat.erpH{ selIxITC(curChan) }(:,erpDat.dsTrlTime>=-.5 & erpDat.dsTrlTime<=5);
                                
                                x1 = tfrDat1h.SxxLFH{selIxITC(curChan)}(tfrDat1h.tx-7>=-.5 & tfrDat1h.tx-7<=5,:);
                                %x2 = tfrDat1m.SxxLFM{selIxITC(curChan)}(tfrDat1m.tx-7>=-.5 & tfrDat1m.tx-7<=5,:);
                                poolPowLF(cntChan,:,:) = (x1-min(min(x1)))./(max(max(x1))-min(min(x1)));
                                
                                x1 = tfrDat2h.SxxHFH{selIxITC(curChan)}(tfrDat2h.tx-7>=-.5 & tfrDat2h.tx-7<=5,:);
                                %x2 = tfrDat2m.SxxHFM{selIxITC(curChan)}(tfrDat2m.tx-7>=-.5 & tfrDat2m.tx-7<=5,:);
                                poolPowHF(cntChan,:,:) = (x1-min(min(x1)))./(max(max(x1))-min(min(x1)));
                            end;
                        end;
                    end;
                else
                    fprintf([ pID{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spk2LFPCoupling_plv_Sorting_alpha0.05_nRand400_lowFreq.mat' ]);
                    fprintf('\n');
                end;
            end;
        end;
    end;
end;

%%
poolITC(cntChan+1:size(poolITC,1),:,:) = [];
poolPowLF(cntChan+1:size(poolPowLF,1),:,:) = [];
poolPowHF(cntChan+1:size(poolPowHF,1),:,:) = [];

%%
figure;

subplot(221);
hold on;
plot(tfrDat1h.fx,squeeze(mean(poolPowLF,2)),'Color',[.75 .75 .75]);
plot(tfrDat1h.fx,squeeze(mean(mean(poolPowLF,2),1)),'k','LineWidth',3);

subplot(222);
x = squeeze(mean(poolPowLF,2));
mF = zeros(size(x,1),1);
for it = 1:size( x,1 )
    [~,mIx] = max( x(it,:) );
    mF(it) = tfrDat1h.fx( mIx );
end;
[n,x] = hist( mF,0:1:30 );
plot(x,n,'k','LineWidth',3);

subplot(223);
hold on;
plot(tfrDat1h.fx,squeeze(mean(poolPowLF(mF >=2 & mF <=20,:,:),2)),'Color',[.75 .75 .75]);
plot(tfrDat1h.fx,squeeze(mean(mean(poolPowLF(mF >=2 & mF <=20,:,:),2),1)),'r','LineWidth',3);

subplot(224);
hold on;
plot(tfrDat2h.fx,squeeze(mean(poolPowHF(mF >=2 & mF <=20,:,:),2)),'Color',[.75 .75 .75]);
plot(tfrDat2h.fx,squeeze(mean(mean(poolPowHF(mF >=2 & mF <=20,:,:),2),1)),'r','LineWidth',3);

%%
figure;
imagesc( tfrDat1h.tx-7, tfrDat1h.fx, squeeze( mean( poolPowLF,1 ) )' );
axis xy;
xlim([-.5 5]);
colormap jet;

figure;
imagesc( tfrDat2h.tx-7, tfrDat2h.fx, squeeze( mean( poolPowHF,1 ) )' );
axis xy;
xlim([-.5 5]);
colormap jet;

%%
figure;
imagesc(itcDat.itcTime(itcDat.itcTime>=-.5 & itcDat.itcTime<=5),itcDat.itcFreq, squeeze( mean( poolITC,1 )));
axis xy;
colormap jet;