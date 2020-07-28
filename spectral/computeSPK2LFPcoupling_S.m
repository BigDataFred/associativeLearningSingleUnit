%% This function computes spike-field coupling; 
% Input: phi = NxM matrix with rows as spike times and columns as
% frequencies; Note: it is essential that the phase - spike matching is
% done outside of this script
% spk2LFPmode = 'ppc' for pairwise phase consistency; 'plv' for phase
% locking value

function [spk2LFPCoupling] = computeSPK2LFPcoupling_S( phi, spk2LFPmode)

[ spk2LFPCoupling ] = NaN(size( phi,2 ),1);
% phi = phi./abs(phi);

switch spk2LFPmode

    case 'plv'
        parfor curFreq = 1:size( phi,2 )
            ph = squeeze(phi(:,curFreq));
            ph(isnan(ph)) = [];
            spk2LFPCoupling(curFreq) = abs(nansum(ph))/length(ph);
        end

        
    case 'ppc'
        
        N = size(phi,1);
        %h1 = waitbar(0, 'COMPUTING SPIKE TO LFP COUPLING');
        
        parfor curFreq =  1:size( phi,2 )
            %waitbar(curFreq/size( phi,2 ), h1)
            ph = phi(:,curFreq);
            ph = angle(ph); % take the angle here, we don't need the vectors anymore
            dum = zeros(N-1, 1);
            for curspk = 1:N-1 
                % hybrid of loop & indexing here to save on RAM & CPU
                dum(curspk) = nansum(cos(ph(curspk) - ph(curspk+1:N)));
            end
            spk2LFPCoupling(curFreq) = nansum(dum) / (N*(N-1)/2);
        end
        %close(h1)
end
%toc
