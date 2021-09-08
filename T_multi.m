clear
clc

load('data_80.mat')

X2=real(ZjC);
Y2=imag(ZjC);

X3=real(ZjW);
Y3=imag(ZjW);

Int2=inpolygon(xx,yy,X2,Y2);
Int3=inpolygon(xx,yy,X3,Y3);

for i=10:108
    
    load(['data_',num2str(i),'.mat'])
    Pd=abs(psiA).^2+abs(psiB).^2;
    
    A1=mean(mean(Pd(101:129,1:125)));
    A2=mean(mean(Pd(101:129,1:500)))+sum(sum(Pd.*Int2))./(sum(sum(Int2)));
    
    B1=mean(mean(Pd0(1:200,1:50)));
    B2=mean(mean(Pd0(1:200,1:200)));
    
    T(i)=(A1/A2)/(B1/B2);
end

semilogy(T)