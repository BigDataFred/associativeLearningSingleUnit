function [Yo,chk] = orthogonalizeTimeDomain( X , Y )

%% This function orthogonalizes the two vectors x and y based on the Gram-Schmidt equation 
% as specified in https://en.wikipedia.org/wiki/Gram–Schmidt_process

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
[Yp] = (sum(X.*Y )./sum(X.*X ))*X; % project signal Y orthogonally onto the line spanned by signal X.
[Yo] = Y - Yp; % subtract Yp to recover signal Yo which now is orthogonal to signal x (i.e. the dot product of Yo and X vectors now is 0)
chk=dot(Yo,X); % This value should be zero (or very close to zero)

