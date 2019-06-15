function [phi,phiTime,phiFreq] = computeLFPphase( lfp, Fs, toi )

addpath('~/tbx/fieldtrip-20170618/');
ft_defaults;

%%
cfgtf                     = [];
cfgtf.method              = 'wavelet';
cfgtf.output              = 'pow';
cfgtf.pad                 = 'nextpow2';
cfgtf.foi                 = 0:30;
cfgtf.width               = 4;
cfgtf.toi                 = toi(1):1/Fs:toi(2);

%%
cfgtf2                     = cfgtf;
cfgtf2.output              = 'fourier';

%%
[ phi ] = cell(1,length( lfp) );
for curChan = 1:length( lfp )
    dum                     = [];
    dum.fsample             = Fs;
    dum.label               = {'dumChan1'};
    dum.trial               = cell(1,size(lfp,2));
    dum.time                = cell(1,size(lfp,2));
    for curTrl = 1:size(lfp,2)
        dum.trial{curTrl}            = gradient( lfp(:,curTrl)' );
        dum.time{curTrl}             = toi(1):1/Fs:toi(2);
    end;
    [phi{curChan}] = ft_freqanalysis( cfgtf2, dum );
end;

%%
for curChan = 1:length(phi)
    if ~isempty(phi{curChan})
        phiTime = phi{curChan}.time;
        phiFreq = phi{curChan}.freq;
        break;
    end;
end;