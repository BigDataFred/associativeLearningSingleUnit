function [spkSelIx] = filterUnits4SPKrate(spkTs,frM, tw,trsh,trlTime)

[ spkSelIx ] = [];
parfor curMW = 1:length( spkTs )
    fprintf([num2str(curMW),'/',num2str(length( spkTs ))]);
    frM2 = frM;
    x = spkTs{curMW}(trlTime >=tw(1) & trlTime<=tw(2),:);
    if ( sum(x(:)) > trsh ) && ( mean(frM2(curMW,trlTime >=tw(1) & trlTime<=tw(2))) > 1 )
        spkSelIx = [spkSelIx curMW];
    end;
    fprintf('\n');
end;