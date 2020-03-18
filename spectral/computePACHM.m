% Compute Phase sorted power spectra for hit and miss trials; because hit
% and miss trials have different dominant frequencies we need to adjust the
% phase and amplitude providing frequency for each trial;

function [pow,phs,rlfp,timeint]=computePACHM(lfphs,rawlfp,hfpow,lffreq,time,toi,hitIdx,nbins);

ntrls1=size(lfphs,1);
ntrls2=size(hfpow,1);
if ntrls1 ~= ntrls2
    display('Error!!! Trial numbers must match');
else
    timeint=-((nbins-1)/2):1:(nbins-1)/2;
    htms=zeros(ntrls1,1);
    htms(hitIdx,1)=1;
    tidx=find(time>=toi(1) & time<=toi(2));
    cyclen(1,1)=1000/lffreq(1);
    cyclen(1,2)=1000/lffreq(2);
    ncyc=5;
    pow=[];
    phs=[];
    rlfp=[];
    cnt=0;
    for nt=1:ntrls1
        if htms(nt,1) == 1
            winbound=round(cyclen(1,1)*ncyc/2);
        else
            winbound=round(cyclen(1,2)*ncyc/2);
        end
        wind=-winbound:winbound;
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

                        

