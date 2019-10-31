function [ix1,ix2] = stratisfyCondIndexes( ix1, ix2 )


if ~isempty(ix1) && ~isempty(ix2)
    
    if ( length( ix1 ) > length( ix2 ) )
        rIx = randperm( length( ix1 ) );% random permute trial indexes
        cnt = 0;
        nsamp = fix( length( ix1 )/ length( ix2 ) );
        tmpIx = [];
        for curSamp = 1:nsamp
            if ~isempty( rIx ) && ( length( rIx) >= length( ix2 ) )
                cnt = cnt+1;
                tmpIx(cnt,:) = ix1( rIx( 1:length( ix2 ) ) );% equate number of spks 4 hits and misses
                rIx(1:length( ix2 )) = [];
            end;
        end;
        ix1 = tmpIx; clear tmpIx;
    elseif ( length( ix1 ) < length( ix2 ) )
        rIx = randperm( length( ix2 ) );% random permute trial indexes
        cnt = 0;
        nsamp = fix( length( ix2 )/ length( ix1 ) );
        tmpIx = [];
        for curSamp = 1:nsamp
            if ~isempty( rIx ) && ( length( rIx) >= length( ix1 ) )
                cnt = cnt+1;
                tmpIx(cnt,:) = ix2( rIx( 1:length( ix1 ) ) );% equate number of spks 4 hits and misses
                rIx(1:length( ix1 )) = [];
            end;
        end;
        ix2 = tmpIx;clear tmpIx;
    end;    
    
end;