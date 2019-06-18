%%
p2d = '/home/rouxf/resultsSpikeFieldOct18II/';

fNs = dir([p2d,'*_fVSpEM_*_spkParams.mat'])

lfpDat.dsTrlTime = -7:1e-3:7;
nP = 0;
pval = 1;
for curFile = 1:length( fNs )
    
    fprintf([num2str(curFile),'/',num2str(length( fNs ))]);
    
    spkParams = load([p2d,fNs(curFile).name]);
    
    ix = min(regexp(fNs(curFile).name,'_'))-1;
    pId = fNs(curFile).name(1:ix);    
    
    ix = regexp(fNs(curFile).name,'\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}');
    sesh = fNs(curFile).name(ix:ix+18);
    
    lfpDat = load([p2d,pId,'_fVSpEM_',sesh,'_spectrogramTrialSegERP_lowFreq.mat']);
    
%     [ spkSelIx2 ] = [];
%     for curMW = 1:length( spkParams.spkTs )
%         x = spkParams.spkTs{curMW}(lfpDat.dsTrlTime >=-0.5 & lfpDat.dsTrlTime<=5,:);
%         x(:,delIx{curMW}) = [];
%         if ( sum(x(:)) > 50 )
%             spkSelIx2 = [spkSelIx2 curMW];
%         end;
%     end;  clear x;
    nSPK = length(spkParams.chanLabSPK);
    nLFP = 0;
    for curLFP = 1:length( lfpDat.chanLabLFP )
        if ~isempty(lfpDat.Ser{curLFP})
            nLFP = nLFP+1;
        end;
    end;
    tmp = nLFP*nSPK;
    p   = 0.05/tmp;
    
    if p<pval
        pval = p;
    end;
    
    nP = nP + tmp;
    fprintf('\n');
    
end;