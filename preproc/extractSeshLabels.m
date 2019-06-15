function [sesh] = extractSeshLabels(tmp)

cnt = 0;
[ sesh ] = cell(1,length(tmp));
for curSesh = 1:length(sesh)
    if ~isempty(regexp(tmp(curSesh).name,'\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}'))
        cnt = cnt+1;
        ix = regexp(tmp(curSesh).name,'\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}');
        ix2 = ix+18;
        sesh(cnt) = {tmp(curSesh).name(ix:ix2)};
    end;
end;
sesh(cnt+1:end) = [];
sesh