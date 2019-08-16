function [Yo] = orthogonalizeTimeDomain( X , Y )

%% proof of principle
% Fs = 1000;
% t = 0:1/Fs:1;
% 
% X = sin(2*pi*10.*t);
% Y1 = real(hilbert(X).*exp(1i*(-pi*1.5)));
% Y2 = Y1+X;
%[Yp] = real(sum(conj(X.*Y2),2)/sum(conj(X.*X),2))*X;
% figure;
% subplot(121);
% plot(t,X);
% subplot(122);
% hold on;
% plot(t,Y1);
% plot(t,Y2-Yp,'rs');
% axis tight;

%%
[Yp] = real(sum(conj(X.*Y),2)/sum(conj(X.*X),2))*X;
[Yo] = Y - Yp;

