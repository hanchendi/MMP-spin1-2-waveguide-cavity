clear
clc

load('data_96.mat')

X2=real(ZjC);
Y2=imag(ZjC);

X3=real(ZjW);
Y3=imag(ZjW);

Int2=inpolygon(xx,yy,X2,Y2);
Int3=inpolygon(xx,yy,X3,Y3);
OutC=zeros(500,500);

for i=1:500
    for j=1:500
        if Int2(i,j)==0
            OutC(i,j)=1;
        end
    end
end

Pd=abs(psiA).^2+abs(psiB).^2;

figure();
pcolor(xx, yy,  (Pd).^(1/8)); hold on;
plot(X2, Y2, 'k--');hold on;plot(X3, Y3, 'k--');
shading flat; axis([-2 2 -2 1.5]); 
colorbar

figure();
pcolor(xx, yy,  imag(psiA)); hold on;
plot(X2, Y2, 'k--');hold on;plot(X3, Y3, 'k--');
shading flat; axis([-2 2 -2 1.5]); 
colorbar

figure();
pcolor(xx, yy,  imag(psiB)); hold on;
plot(X2, Y2, 'k--');hold on;plot(X3, Y3, 'k--');
shading flat; axis([-2 2 -2 1.5]); 
colorbar

jx=psiA.*conj(psiB)+psiB.*conj(psiA);
jy=sqrt(-1)*psiA.*conj(psiB)-sqrt(-1)*psiB.*conj(psiA);

jx1=(psiA.*conj(psiB)+psiB.*conj(psiA)).*Int3;
jy1=(sqrt(-1)*psiA.*conj(psiB)-sqrt(-1)*psiB.*conj(psiA)).*Int3;

jz=psiA.*conj(psiA)-psiB.*conj(psiB);

m=10;
figure();
quiver(xx(1:m:500,1:m:500), yy(1:m:500,1:m:500),  jx(1:m:500,1:m:500), jy(1:m:500,1:m:500)); hold on;
plot(X2, Y2, 'k--');hold on;plot(X3, Y3, 'k--');
shading flat; axis([-2 2 -2 1.5]); 
colorbar

miuB=abs(sum(sum((jx.*yy-jy.*xx).*Int2))/sum(sum(abs(Pd).*Int2)));

figure();
pcolor(xx, yy,  jz.*Int3); hold on;
plot(X2, Y2, 'k--');hold on;plot(X3, Y3, 'k--');
shading flat; axis([-2 2 -2 1.5]); 
colorbar

jz_sum=zeros(1,500);

for i=1:500
    jz_sum(i)=sum(jz(101:129,i));
end

figure()
plot(jz_sum)