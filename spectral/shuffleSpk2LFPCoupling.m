function [spk2LFPRand] = permutationStats4Spk2LFPcoupling( phi, spkIx,  phiTrl, fIx, nRand, spk2LFPmode )
 
[ spk2LFPRand ] = NaN( length( fIx), nRand );
parfor randIter = 1:nRand
    
    dum = NaN(1,length(fIx));
    rIx = randperm(phiTrl);
    
    switch spk2LFPmode
        
        case 'plv'
            for curFreq = 1:length( fIx)
                ph = squeeze(phi(:,curFreq,:))';
                ph = ph(:,rIx);
                ph = ph(:);
                ph = ph(spkIx);
                ph(isnan(ph)) = [];
                ph = ph./abs(ph);
                dum(curFreq) = abs(nansum(ph))/length(ph);
            end;
            spk2LFPRand(:,randIter) = dum;
            
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
                        
            for curFreq = 1:length( fIx)
                ph = squeeze(phi(:,curFreq,:))';
                ph = ph(:,rIx);
                ph = ph(:);
                ph = ph(spkIx);
                d = diff(angle(ph(p)),[],2);
                dum(curFreq) = nansum(cos(d))/size(p,1);
            end;
            spk2LFPRand(:,randIter) = dum;
            
    end;
    
end;
