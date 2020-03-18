%% This function computes the modulation index as described in Tort et al. (2010) J Neurophysiol
% input: 
% lfphase = 3 dimensional Matrix containing phase angles organized as [trls, freq, time];
% hfpow   = 3 dimensional Matrix containing power organized as [trls, freq, time];
% edges   = vector contain left limits for phase bins; i.e. linspace(-pi,pi, 18)
% nshuff  = number of shuffles (i.e. 200)

function [miz,mir,mish]=computePAC_MI(lfphase,hfpow,edges,nshuff);

nphs=size(lfphase,2);
npow=size(hfpow,2);
ntrls=size(lfphase,1);

% Calculate MI for real data
for psi=1:nphs
    for poi=1:npow
        tmpphs=squeeze(lfphase(:,psi,:));
        tmppow=squeeze(hfpow(:,poi,:));
        [mir(poi,psi),~]=ModIndex_v2(tmpphs, tmppow, edges(1:end-1));
    end
end

% Calculate MI for shuffled data; NB only trials for phase will be shuffled
mish=zeros(npow,nphs,nshuff);
parfor ns=1:nshuff
    shidx=randperm(ntrls);
    for psi=1:nphs
        for poi=1:npow
            tmpphs=squeeze(lfphase(shidx,psi,:));% shuffle phase trials
            tmppow=squeeze(hfpow(:,poi,:));
            [mish(poi,psi,ns),~]=ModIndex_v2(tmpphs, tmppow, edges(1:end-1));
        end
    end
end

% Calculate MI as z-value compared to shuffled data
miz=(mir-mean(mish,3))./std(mish,0,3);

