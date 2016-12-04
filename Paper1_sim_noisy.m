%Field detection from location unaware sensors on a discrete grid
%Clustering using EM algorithm in the presence of Gaussian noise
%Code optimised for large bandwidths (upto b=10) and sample sizes
tic;
clc
close all
clear all
b=10; %Signal Bandwidth
kb=2*b+1;


%% Initialising Distribution of Sample locations
v=1:1:2*b+1;
p=(v.^2)*3/((b+1)*(2*b+1)*(4*b+3));

%% Field Estimation from noisy samples
sig=0.05; %Standard deviation of noise
a=zeros(1,2*b+1);
a(b+1)=1;
err=zeros(3,10000);
dist=zeros(3,10000);
M=zeros(3,10000);
iter=zeros(3,10000);
for k=3:5
    n=10^k;
    S=0;
    %10000 monte-carlo simulations
    for i=1:10000
        a(1:b)=rand(1,b)+(rand(1,b))*1i;
        %a(b+2:2*b+1)=flipud(conj(a(1:b))')';
        a(b+2:end)=conj(a(b:-1:1));
        xs=real(fft(ifftshift(a')));
        samples=datasample(xs,n,'Replace',true,'Weights',p);
        samples=samples+sig*randn(size(samples));
        [yb, I]=segment_em(samples,kb,sig);
        %disp(I) %Number of iterations for convergence
        a1=(fftshift(ifft(yb)))';
        S=S+((norm(a-a1).^2)/(norm(a).^2));
        iter(k-2,i)=I;
        dist(k-2,i)=min(pdist(xs).^2); %Min. pairwise squared distance b/w cluster means
        err(k-2,i)=S/i; %Average error
        M(k-2,i)=(norm(a-a1).^2)/(norm(a).^2);
    end
    figure
    hist(M(k-2,:),100);
    title(strcat('n=',num2str(n)));
    figure
    hist(dist(k-2,:),100);
    title(strcat('n=',num2str(n)));
    % figure
    % plot(err);
    % title(strcat('n=',num2str(n)));
    disp(n); %Sample size
    disp(S/10000); %Average distortion
    disp(min(M(k-1,:))); %Minimum distortion
end
toc;