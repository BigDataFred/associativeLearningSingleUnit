function [sigIxLFP,sigIxSPK] = computeSpk2LFPCouplingSig( phi, spkTs, nRand, fIx, tIx, delIx,  spk2LFPmode, alpha )

%%
[ spk2lfp ] = zeros(length(phi),length(spkTs));
[ spk2lfpRnd ] = zeros(length(phi),length(spkTs),nRand);

%%
for curMW = 1:length(phi)
    fprintf([num2str(curMW),'/',num2str(length(phi))]);
    
    if ~isempty(phi{curMW})
        [ tmp ] = squeeze(phi{curMW}(:,:,fIx,tIx));
        [ phiTrl ] = size(tmp,1);
        
        [ mxRnd ] = zeros(length(spkTs),nRand);
        [ mX ] = zeros(1,length(spkTs));
        
        for curMW2 = 1:length( spkTs )
            ts = spkTs{curMW2}(tIx,:);
            ts(:,delIx{curMW}) = [];
            ix = find(ts==1);
            ts = [];
            
            [tmp2] = computeSPK2LFPcoupling( tmp, ix, spk2LFPmode );
            mX(curMW2) = max( tmp2 );
            
            [tmp2] = shuffleSpk2LFPCoupling( tmp, ix,  phiTrl, nRand, spk2LFPmode );
            [mxRnd(curMW2,:)] = max( tmp2,[],1 );
        end;
        [spk2lfp(curMW,:)] = mX;            mX      = [];
        [spk2lfpRnd(curMW,:,:)] = mxRnd;    mxRnd   = [];
        
    end;
    fprintf('\n');
end;

%%
[ Pval ] = zeros(size(spk2lfp));
for curMW = 1:size(spk2lfp,1)
    for curMW2 = 1:size(spk2lfp,2)
        Pval(curMW,curMW2) = length(find(spk2lfpRnd(curMW,curMW2,:) >= spk2lfp(curMW,curMW2)))/nRand;
    end;
end;
[sigIxLFP,sigIxSPK] = find( Pval < (alpha/length(Pval(:))) );