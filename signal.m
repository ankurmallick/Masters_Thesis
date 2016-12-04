function [g,a]=signal(b)
%Function to generate a real valued signal 'g' of bandwidth 'b'
%'a' is the vector of Fourier Series Coefficients
a=zeros(1,2*b+1);
a(1:b)=rand(1,b)+(rand(1,b))*1i;
a(b+1)=1;
a(b+2:2*b+1)=flipud(conj(a(1:b))')';
t=(0:1/1000:0.999)';
th=0+((-b:1:b)*pi)*2i;
fs=exp(t*th);
rep=ones(1000,1);
g=sum((rep*a).*fs,2);
% g=zeros(1,1000);
% for j=1:1000
%     t=(j-1)/1000;
%     th=0+((-b:1:b)*pi*t)*2i;
%     g(j)=sum(a.*exp(th));
% end
g=real(g);
end