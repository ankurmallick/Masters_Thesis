%Author: Ankur Mallick
%To estimate the number of samples required for 1% detection error probability
close all
clear all
clc
bw=[3, 5, 10, 20];
N1=zeros(1,4);
Eplot=zeros(1,4);
for j=1:4
    b=bw(j);
    v=0:1:2*b;
    v=v+1;
    p=(v.^2)*3/((b+1)*(2*b+1)*(4*b+3));
    str=strcat('BW',num2str(b),'.mat');
    load(str);
    xs=real(fft(ifftshift(a')));
    E=err(1,:);
    I=find(E==0.01,1);
    if(isempty(I))
        d=E-0.01;
        [M1, i1]=min(d(d>0));
        [M2, i2]=max(d(d<0));
        n1=n(d>0);
        m1=n1(i1);
        n2=n(d<0);
        m2=n2(i2);
        disp([m1, m2])
        L=m1;
        U=m2;
        while (L<U)
            m=round((L+U)/2);
            pe=0;
            for i=1:10000
                [samples,idx]=datasample(xs,m,'Replace',true,'Weights',p);
                [h,x]=hist(samples,unique(samples));
                [h1,I1]=sort(h);
                gb=x(I1);
                if(length(gb)<2*b+1)
                    %The case where there are no samples at 1 location
                    pe=pe+1;
                else
                    %a1=(phi1*gb);
                    if(norm((gb-xs))>0)
                        %The case where the number of samples do not follow the
                        %expected ordering thus leading to errors
                        pe=pe+1;
                    end
                end
            end
            pe=pe/10000;
            if(pe>0.015)
                L=m+1;
            elseif(pe<0.005)
                U=m-1;
            else
                N1(j)=m;
                Eplot(j)=pe;
                disp(pe);
                disp(m);
                break;
            end
        end
    end
end
plot(bw,N1);
xlabel('Bandwidth');
ylabel('Sample Size');
save('Plot1e.mat','bw','N1','Eplot');