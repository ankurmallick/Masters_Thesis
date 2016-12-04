%Field detection from location unaware sensors on a discrete grid
%Comparing performance for different distributions on the sensor locations
%No measurement noise
%Code optimised for large bandwidths (upto b=20) and sample size(upto 10^6)
close all
clear all
b=20;
%b=3; %Signal bandwidth
sb=1/(2*b+1);

%% Generating the Signal
a=zeros(1,2*b+1);
a(b+1)=1;
a(1:b)=rand(1,b)+(rand(1,b))*1i;
a(b+2:end)=conj(a(b:-1:1));
xs=real(fft(ifftshift(a')));
% [g,a]=signal(b);
% a=conj(a');
% figure
% plot(g);
% title('Original Signal');
% K=(1+round((0:1:2*b)*sb*1000))';
% xs=g(K); %Actual values at sampling locations

%% Monte Carlo Simulations on Sample Locations
n=[100:100:1000, 2000:1000:10000, 20000:10000:100000, 2*10^5:1*10^5:10^6]; %Change this to change the number of samples at each iteration
err=zeros(4,length(n)); %err(k,j)=Probability of error for the kth distrubution with number of samples=n(j) 
P=zeros(4,2*b+1); %P(k,:)=kth distribution vector
%probs=zeros(1,length(g));
v=0:1:2*b;
%pos=1+round(v*sb*1000);
v=v+1;
for k=1:4
    if k==1
        %Optimal Distribution
        p=(v.^2)*3/((b+1)*(2*b+1)*(4*b+3));
        cdf=cumsum(p);
    elseif k==2
        %Linear Distribution
        p=v./sum(v);
    elseif k==3
        %Cubic Distribution
        v1=v.^3;
        p=v1./sum(v1);
    else
        %Random Distribution
        p=rand(1,2*b+1);
        p=p./sum(p);
        p=sort(p);
        cdf=cumsum(p);
    end
    P(k,:)=p;
    %probs(pos)=p;
    for j=1:length(n)
        disp(j);
        for i=1:10000
            %gb=sample(g,b,cdf,n(j)); %Estimated signal
            [samples,idx]=datasample(xs,n(j),'Replace',true,'Weights',p);
            [h,x]=hist(samples,unique(samples));
            [h1,I1]=sort(h);
            gb=x(I1);
            if(length(gb)<2*b+1)
                %The case where there are no samples at 1 location
                err(k,j)=err(k,j)+1;
            else
                %a1=(phi1*gb);
                if(norm((gb-xs))>0)
                    %The case where the number of samples do not follow the
                    %expected ordering thus leading to errors
                    err(k,j)=err(k,j)+1;
                end
            end
        end
    end
    disp(k);
end
err=err/10000;
figure
loglog(n,err(1,:),'-k',n,err(2,:),'--k',n,err(3,:),':k',n,err(4,:),'-.k');
xlabel('Number of Samples');
ylabel('Empirical detection error probability');
legend('optimal','linear','cubic','random');
str=strcat('BW',num2str(b),'.jpg');
print(str,'-djpeg');
str=strcat('BW',num2str(b),'.mat');
save(str,'a','n','err');
%toc;
