

%%
Fs = 1e3;
t = 0:1/Fs:1;

x1 = zeros(100,length(t));
x2 = zeros(100,length(t));
for it = 1:100
    x1(it,:) = sin(2*pi*5.*t);% + 1.*randn(1,length(t));% pseudo common
    x2(it,:) = cos(2*pi*5.*t) + 1.*randn(1,length(t));% pseudo signal of interest
end;

%x2 = x2+x1;

%Yo = zeros(100,length(t));
%for it = 1:size(x1,1)
%    [Yo(it,:)] = orthogonalizeTimeDomain(mean(x1,1) , x2(it,:) );
%end

[Yo2] = orthogonalizeTimeDomain( ones(size(x1,1),1)*mean(x1,1) , x2 );

%%
figure;
subplot(231);
imagesc(t,1:size(x2,1),x2);
colorbar;
subplot(232);
imagesc(t,1:size(x2,1),x2-ones(size(x1,1),1)*mean(x1,1));
colorbar;
subplot(233);
imagesc(t,1:size(x2,1),Yo2);
colorbar;

subplot(234);
plot(t,mean(x2,1));
subplot(235);
plot(t,mean(x2-ones(size(x1,1),1)*mean(x1,1),1));
subplot(236);
plot(t,mean(Yo2,1));
