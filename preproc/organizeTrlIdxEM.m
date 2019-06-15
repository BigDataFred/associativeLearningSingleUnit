function [trlPool,hitIdx,missIdx,trlENC] = organizeTrlIdxEM(lfpDat)

[ trlPool ] = sort([lfpDat.hitIdx;lfpDat.missIdx]);
[ hitIdx ]  = lfpDat.hitIdx;
[ missIdx ] = lfpDat.missIdx;
[ trlENC ]  = 1:length(lfpDat.trlENC);

if ~all(trlENC' == trlPool) || (length(trlPool) == length(trlENC))~=1 || ~all(hitIdx == trlPool(ismember(trlPool,hitIdx))) || ~all(sort(missIdx) == trlPool(ismember(trlPool,missIdx)))
    error('error wrong trial assignment');
end;