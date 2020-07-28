function [pow,phs,rlfp,timeint]=computePAC(lfphs,rawlfp,hfpow,lffreq,time,toi,nbins);

ntrls1=size(lfphs,1);
ntrls2=size(hfpow,1);
tidx=find(time>=toi(1) & time<=toi(2));
cyclen=1000/lffreq;
ncyc=5;
winbound=round(cyclen*ncyc/2);
wind=-winbound:winbound;
pow=[];
phs=[];
rlfp=[];
cnt=0;
timeint=-((nbins-1)/2):1:(nbins-1)/2;
if ntrls1 ~= ntrls2
    display('Error!!! Trial numbers must match');
else
    
    for nt=1:ntrls1
        [tmppks,tmplocs]=findpeaks(lfphs(nt,:),'MinPeakHeight',2.9);
        tmpidx=find(ismember(tmplocs,tidx));
        selpks=tmplocs(tmpidx);
        for np=1:length(tmpidx)
            cnt=cnt+1;
            selwind=selpks(np)+wind;
            tmppow=squeeze(hfpow(nt,:,selwind));
            tmpowmn=repmat(mean(tmppow,2),[1 length(wind)]);
            tmpowsd=repmat(std(tmppow,[],2), [1 length(wind)]);
            tmpowz=(tmppow-tmpowmn)./tmpowsd;
            tmpphs=cos(lfphs(nt,selwind));
            tmplfp=rawlfp(nt,selwind);
            tmppowzi=imresize(tmpowz,[size(tmpowz,1) nbins]);
            tmpphsi=imresize(tmpphs,[1 nbins]);
            tmplfpi=imresize(tmpphs,[1 nbins]);
            pow=sum(cat(3,pow,tmppowzi),3);
            phs=sum([phs;tmpphsi],1);
            rlfp=sum([rlfp;tmplfpi],1);
        end
    end
end

pow=pow./cnt;
phs=phs./cnt;
rlfp=rlfp-mean(rlfp);
rlfp=rlfp./abs(min(rlfp));

                        

