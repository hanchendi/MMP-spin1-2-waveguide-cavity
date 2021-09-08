clear
clc

load('data_80.mat')

X2=real(ZjC);
Y2=imag(ZjC);

X3=real(ZjW);
Y3=imag(ZjW);

Int2=inpolygon(xx,yy,X2,Y2);
Int3=inpolygon(xx,yy,X3,Y3);

for i=80:97
    
    load(['data_',num2str(i),'.mat'])
    Pd=abs(psiA).^2+abs(psiB).^2;
    
    DOS(i)=sum(sum(Pd.*Int2));
    
end

semilogy(DOS)