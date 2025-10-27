clear; clc; close all;

Properties1=[0.008 0.99 420*10^6 620*10^6 0.1 200*10^9 0 ];
Properties2=[0 0.99 1400*10^6 1850*10^6 0.06 165*10^9 0];

Strain1=0:0.001:Properties1(5);
for c=1:1:length(Strain1)
Stress1(c)=STEEL(Strain1(c),Properties1);
end

Strain2=0:0.001:Properties2(5);
for c=1:1:length(Strain2)
Stress2(c)=STEEL(Strain2(c),Properties2);
end

plot(Strain1,Stress1,'b',Strain2,Stress2,'r'),xlabel('Strain'),ylabel('Stress (Pa)');
legend('Steel Bar','Steel Strand');