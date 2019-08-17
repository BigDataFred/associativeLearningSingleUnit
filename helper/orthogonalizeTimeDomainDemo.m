

%%
Fs = 1e3;
t = 0:1/Fs:1;

x1 = zeros(100,length(t));
x2 = zeros(100,length(t));
for it = 1:100
    x1(it,:) = sin(2*pi*5.*t) + 1.*randn(1,length(t));% pseudo common
    x2(it,:) = cos(2*pi*5.*t) + 1.*randn(1,length(t));% pseudo signal of interest
end;

x1(55,250) = x1(55,250)+1e4;
x1(65,500) = x1(65,500)+1e4;
x1(75,750) = x1(75,750)+1e4;

mxS = x2+x1;

Yo2 = zeros(100,length(t));
for it = 1:size(x1,1)
    [Yo2(it,:)] = orthogonalizeTimeDomain(mean(x1,1) , x2(it,:) );
end

%%
figure;
subplot(231);
imagesc(t,1:size(x2,1),mxS);
colorbar;
subplot(232);
imagesc(t,1:size(x2,1),mxS-ones(size(x1,1),1)*mean(x1,1));
colorbar;
subplot(233);
imagesc(t,1:size(x2,1),Yo2);
colorbar;

subplot(234);
plot(t,mean(mxS,1));
subplot(235);
plotyy(t,mean(mxS-ones(size(x1,1),1)*mean(x1,1),1),t,mean(x2,1));
subplot(236);
hold on;
%plot(t,mean(x2,1),'r');
plotyy(t,mean(Yo2,1),t,mean(x2,1));
