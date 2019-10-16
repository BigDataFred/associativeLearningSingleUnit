
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
poolPow1 = zeros(947,550,246);
poolPow2 = zeros(947,550,1147);

for curPat = 1:length( pID )
    for curExp = 1:length( expMode )
        
        [ tmp ] = dir( [ rdsPath1, pID{curPat},'/', expMode{curExp}, '/' ] );
        
        if ~isempty(tmp)        
            
            [ sesh ] = extractSeshLabels( tmp );
            clear tmp;
            
            for curSesh = 1:length( sesh )             
                
                [ itcN ] = dir( [ rdsPath2, pID{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spectrogramTrialSegITC_lowFreq.mat' ] );
                
                [ itcDat ] = load([rdsPath2,itcN.name]);
                
                %if any( itcDat.pval2 <= 1e-2 ) %0.05/length( itcDat.pval1
                    
                    [ selIx ] = 1:length(itcDat.pval2);%find( itcDat.pval2 <= 1e-2 );
                    
                    [ lfN ] = dir([ rdsPath2, pID{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spectrogramTrialSeg_lowFreq.mat' ]);
                    [ hfN ] = dir([ rdsPath2, pID{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spectrogramTrialSeg_highFreq.mat' ]);
                    
                    [tfrDat1] = load([rdsPath2,lfN.name]);
                    [tfrDat2] = load([rdsPath2,hfN.name]);                                        
                    
                    for curChan = 1:length( selIx )
                        cntChan = cntChan +1;
                        
                        x = tfrDat1.SxxLF{selIx(curChan)};
                        x = x(tfrDat1.tx-7>=-.5 & tfrDat1.tx-7<=5,:);
                        x = (x-min(min(x)))./(max(max(x))-min(min(x)));
                        poolPow1(cntChan,:,:) = x;
                        
                        x = tfrDat2.SxxHF{selIx(curChan)};
                        x = x(tfrDat2.tx-7>=-.5 & tfrDat2.tx-7<=5,:);
                        x = (x-min(min(x)))./(max(max(x))-min(min(x)));
                        poolPow2(cntChan,:,:) = x;
                    end;
                %end;
            end;
        end;
    end;
end;
poolPow1(cntChan+1:size(poolPow1,1),:,:) = [];
poolPow2(cntChan+1:size(poolPow2,1),:,:) = [];

%%
figure;

subplot(221);
hold on;
plot(tfrDat1.fx,squeeze(mean(poolPow1,2)),'Color',[.75 .75 .75]);
plot(tfrDat1.fx,squeeze(mean(mean(poolPow1,2),1)),'k','LineWidth',3);

subplot(222);
x = squeeze(mean(poolPow1,2));
mF = zeros(size(x,1),1);
for it = 1:size( x,1 )
    [~,mIx] = max( x(it,:) );
    mF(it) = tfrDat1.fx( mIx );
end;
[n,x] = hist( mF,0:1:30 );
plot(x,n,'k','LineWidth',3);

subplot(223);
hold on;
plot(tfrDat1.fx,squeeze(mean(poolPow1(mF >=2 & mF <=20,:,:),2)),'Color',[.75 .75 .75]);
plot(tfrDat1.fx,squeeze(mean(mean(poolPow1(mF >=2 & mF <=20,:,:),2),1)),'r','LineWidth',3);

subplot(224);
hold on;
plot(tfrDat2.fx,squeeze(mean(poolPow2(mF >=2 & mF <=20,:,:),2)),'Color',[.75 .75 .75]);
plot(tfrDat2.fx,squeeze(mean(mean(poolPow2(mF >=2 & mF <=20,:,:),2),1)),'r','LineWidth',3);

%%
figure;
imagesc( tfrDat1.tx-7, tfrDat1.fx, squeeze( mean( poolPow1(mF >=2 & mF <=20,:,:),1 ) )' );
axis xy;
xlim([-.5 5]);
colormap jet;

figure;
imagesc( tfrDat2.tx-7, tfrDat2.fx, squeeze( mean( poolPow2(mF >=2 & mF <=20,:,:),1 ) )' );
axis xy;
xlim([-.5 5]);
colormap jet;