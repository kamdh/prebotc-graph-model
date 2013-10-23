clear VV spikeV lenABovN dfAbovN SpikesNum 
clear AbovN spNNdfAbovN spikeNN spikeNaux spikeTime dfAbovN dfAbovN
clear NN ti tinrv
thresh=0.0;
minInterSpikeDurN=1;
r1=10; %number of active neurons 

%two neutons have 10 variables 
for n=1:1:2
    s=1+10*(n-1);
    VV(:,n)=yy(:,s);
end
%eight neurons have 8 variables

for n1=3:1:10
    s1=-3+8*n1;
    VV(:,n1)=yy(:,s1);
end


%______________________________________________________

for n=1:1:r1
    clear V;
    clear spikeV lenABovN dfAbovN SpikesNum 
    clear AbovN spNNdfAbovN spikeNN spikeNaux spikeTime dfAbovN dfAbovN

    V=VV(:,n);
    AbovN=find(V>thresh);
    lenAbovN=length(AbovN);
    dfAbovN=AbovN(2:lenAbovN)-AbovN(1:lenAbovN-1);
    spNNdfAbovN=find( dfAbovN>minInterSpikeDurN);


    SpikesNum=length(spNNdfAbovN);



    for i=1:SpikesNum-1
        [spikeV(i) spikeNaux(i)]=max(V(AbovN(spNNdfAbovN(i)+1):AbovN(spNNdfAbovN(i+1))));
    end
    spikeNN=spikeNaux(:)+AbovN(spNNdfAbovN(1:SpikesNum-1)+1)-1;

    spikeTime=t(spikeNN);

    % subplot(2,1,1);
    % plot(t,V,'r');
    % hold on;
    % plot(spikeTime,spikeV,'g+');

    tinrv=40*(t(3)-t(1));
    ti=0:tinrv:t(length(t));
    % subplot(2,1,2)
    for k=2:1:length(ti)
        SpikePerInt=find (spikeTime > ti(k-1) & spikeTime < ti(k));
        NSpikePerInt(n,k)=length(SpikePerInt);
        %     hold on;
    
    %Total number of spikes in populations
    %Ntotal(1)=0;
    %Ntotal=Ntotal(n)+NSpikePerInt
    
    %     plot(ti,NSpikePerInt,'b*');
    %     hold on;
    %    % clear NSpikePerInt;
    
end

for j=1:1:length(ti)
    NN(j)=sum(NSpikePerInt(:,j));
end
%plot(ti,NN,'.');
%subplot(r,1,5)
figure(5);
bar(ti,NN,'g');
axis([ 0. 20. 2.5 70])
hold on;
%Filter
order=2;
dt=ti(1,2)-ti(1,1);
Ny=1/dt*0.5;
coF=4;
coF1=coF/Ny;
[b c]=butter(order, coF1,'low');
filtNs=filtfilt(b, c, NN(1,:));

