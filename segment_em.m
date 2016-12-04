%Author: Ankur Mallick
function [yb, I]=segment_em(samples,kb,sig)
%Clustering noisy samples via the EM algorithm
n=length(samples);
u=zeros(1,kb);
u(1)=samples(1);
%K-Means++ algorithm for initial guess of cluster means
for i=2:kb
    v=u(1:i-1);
    d1=(samples*ones(1,length(v))-ones(n,1)*v).^2;
    M=min(d1,[],2);
    u(i)=datasample(samples,1,'Weights',M.^2);
%     [M1,I]=max(M);
%     u(i)=samples(I);
end
w=ones(1,kb)/kb;
C=sig*sig*ones(1,kb);
th=10^-5;
s1=samples*ones(1,kb);
u1=ones(n,1)*u;
w1=ones(n,1)*w;
C1=ones(n,1)*(C);
d=(s1-u1).^2;
for i=1:100
%      disp(i);
%      disp(u);
%      disp(C);
%      disp(w);
    g=exp(d./(-2*C1));
    g=g./sqrt(2*pi*C1);
    g=g.*w1;
    %disp(g);
    S=sum(g,2)*ones(1,kb);
    %disp(S);
    mems=g./S;
%     disp(sum(mems));
    u_old=u;
    u=(sum(s1.*mems))./(sum(mems));
    u1=ones(n,1)*u;
    d=(s1-u1).^2;
    %C=(sum(d.*mems))./(sum(mems));
    w=sum(mems)/n;
    w1=ones(n,1)*w; 
    %C1=ones(n,1)*(C);
    if((norm(u-u_old)^2)<th)
        break;
    end
end
% disp('Means');
% disp(u);
% disp('Weights')
% disp(w);
% disp('Variances')
% disp(C);
I=i;
[s1,I1]=sort(w);
yb=u(I1)';
%u is a column vector
%w is a row vector
end
