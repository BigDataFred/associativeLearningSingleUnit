function [spk2LFPCoupling] = computeSPK2LFPcoupling( phi, spkIx, spk2LFPmode)

[ spk2LFPCoupling ] = NaN(size( phi,2 ),1);

switch spk2LFPmode
    case 'plv'
        parfor curFreq = 1:size( phi,2 )
            ph = squeeze(phi(:,curFreq,:))';
            ph = ph(:);
            ph = ph(spkIx);
            ph(isnan(ph)) = [];
            ph = ph./abs(ph);
            spk2LFPCoupling(curFreq) = abs(nansum(ph))/length(ph);
        end;
    case 'ppc'
        N = length(spkIx);
        p = zeros(N*(N-1)/2,2);
        cnt = 0;
        for curTrl = 1:N-1
            for curTrl2 = curTrl+1:N
                cnt= cnt+1;
                p(cnt,:) = [curTrl curTrl2];
            end;
        end;
        
        parfor curFreq =  1:size( phi,2 )
            ph = squeeze(phi(:,curFreq,:))';
            ph = ph(:);
            ph = ph(spkIx);
            d = diff(angle(ph(p)),[],2);
            spk2LFPCoupling(curFreq) = nansum(cos(d))/size(p,1);
        end;
end;
