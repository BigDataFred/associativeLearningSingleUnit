function [spkSelIx] = filterUnits4SPKrate(spkTs,frM, tw,trsh,trlTime)

[ spkSelIx ] = [];
parfor curMW = 1:length( spkTs )
    fprintf([num2str(curMW),'/',num2str(length( spkTs ))]);
    frM2 = frM;
    tw2 = tw;
    x = spkTs{curMW}(trlTime >= tw2(1) & trlTime<=tw2(2),:);
    if ( sum(x(:)) > trsh ) && ( mean(frM2(curMW,trlTime >=tw2(1) & trlTime<=tw2(2))) > 1 )
        spkSelIx = [spkSelIx curMW];
    end;
    fprintf('\n');
end;