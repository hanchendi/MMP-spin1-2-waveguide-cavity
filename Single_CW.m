function [xx yy ZjC ZjW ee psiA psiB] = Single_CW(kx,E)

V2=10;
V3=8;

a=0.1;
b=5;
d=0.1;
deviate_W=1+a+d;

s1=sign(E);
s2=sign(E-V2);
s3=sign(E-V3);

k1=abs(E);
k2=abs(E-V2);
k3=abs(E-V3);

%%%%%%%%%%%%%%%%%%%%%%%
% Circle
%%%%%%%%%%%%%%%%%%%%%%%

Rm=0.95;
NmC=4*60;
theta_mC=linspace(0,2*pi,NmC);
ZmC=Rm*cos(theta_mC)+sqrt(-1)*Rm*sin(theta_mC);

Rl=1.05;
NlC=4*64;
theta_lC=linspace(0,2*pi,NlC);
ZlC=Rl*cos(theta_lC)+sqrt(-1)*Rl*sin(theta_lC);

NjC=(NlC+NmC)*5; %%%%%%%%%%%%% boundary
theta_jC=linspace(0,2*pi,NjC);
ZjC=cos(theta_jC)+sqrt(-1)*sin(theta_jC);



%%%%%%%%%%%%%%%%%%%%%%%
% Waveguide
%%%%%%%%%%%%%%%%%%%%%%%

NlW=4*204;
delta_l=0.02;
t1=floor(NlW/2*a/b);

xl(1:NlW/2)=linspace(-b-delta_l,b+delta_l,NlW/2);
yl(1:NlW/2)=(-a-delta_l)*ones(1,NlW/2);

xl=[xl (b+delta_l)*ones(1,t1)];
yl=[yl linspace(-a-delta_l,a+delta_l,t1)];

xl=[xl linspace(b+delta_l,-b-delta_l,NlW/2)];
yl=[yl (a+delta_l)*ones(1,NlW/2)];

xl=[xl (-b-delta_l)*ones(1,t1)];
yl=[yl linspace(a+delta_l,-a-delta_l,t1)];

yl=yl-deviate_W;
ZlW=xl+yl*sqrt(-1);
NlW=length(xl);

NmW=4*200;

delta_m=0.02;
t2=floor(NmW/2*a/b);
theta_l=linspace(0,2*pi,NlW);
xm(1:NmW/2)=linspace(-b+delta_m,b-delta_m,NmW/2);
ym(1:NmW/2)=(-a+delta_m)*ones(1,NmW/2);

xm=[xm (b-delta_m)*ones(1,t2)];
ym=[ym linspace(-a+delta_m,a-delta_m,t2)];

xm=[xm linspace(b-delta_m,-b+delta_m,NmW/2)];
ym=[ym (a-delta_m)*ones(1,NmW/2)];

xm=[xm (-b+delta_m)*ones(1,t2)];
ym=[ym linspace(a-delta_m,-a+delta_m,t2)];

ym=ym-deviate_W;
ZmW=xm+ym*sqrt(-1);
NmW=length(xm);

NjW=(NlW+NmW)*5; %%%%%%%%%%%%% boundary

t3=floor(NjW/2*a/b);
XW(1:NjW/2)=linspace(-b,b,NjW/2);
YW(1:NjW/2)=(-a)*ones(1,NjW/2);

XW=[ XW (b)*ones(1,t3)];
YW=[ YW linspace(-a,a,t3)];

XW=[ XW linspace(b,-b,NjW/2)];
YW=[YW (a)*ones(1,NjW/2)];

XW=[ XW (-b)*ones(1,t3)];
YW=[ YW linspace(a,-a,t3)];

YW=YW-deviate_W;
ZjW=XW+YW*sqrt(-1);
NjW=length(XW);

% Construct matrix

zA=[ZmC ZmW];
zB=ZlC;
zC=ZlW;

[zl, zjl] = meshgrid(zA, ZjC);
Djl = zjl - zl; 
Phi_AC = angle(Djl);
R_AC = abs(Djl);

[zl, zjl] = meshgrid(zB, ZjC);
Djl = zjl - zl; 
Phi_BC = angle(Djl);
R_BC = abs(Djl);

[zl, zjl] = meshgrid(zA, ZjW);
Djl = zjl - zl; 
Phi_AW = angle(Djl);
R_AW = abs(Djl);

[zl, zjl] = meshgrid(zC, ZjW);
Djl = zjl - zl; 
Phi_CW = angle(Djl);
R_CW = abs(Djl);

%% Incidence

L=1;
nL=-L:L;

alpha1=sqrt(kx^2-(E-V3)^2);
alpha2=sqrt(kx^2-E^2);

A_norm=cos(imag(alpha1)*a)/exp(-alpha2*a);

plnAC=A_norm*exp(sqrt(-1)*kx*real(ZjC)).*exp(-alpha2*abs(imag(ZjC)+deviate_W));
plnBC=A_norm*s3*exp(sqrt(-1)*kx*real(ZjC)).*exp(-alpha2*abs(imag(ZjC)+deviate_W));

plnAW = exp(sqrt(-1)*kx*real(ZjW)).*cos(imag(alpha1)*(imag(ZjW)+deviate_W));
plnBW = s3*exp(sqrt(-1)*kx*real(ZjW)).*cos(imag(alpha1)*(imag(ZjW)+deviate_W));

pln=conj([plnAC plnBC plnAW plnBW]');

%% MMP

AL1=[];AL2=[];
AR1=[];AR2=[];

BL1=[];BL2=[];
CR1=[];CR2=[];

for l=1:length(nL)    
    
        HjmL1 = besselh(nL(l), 1, k1*R_AC).*exp(sqrt(-1)*nL(l)*Phi_AC);
        HjmL2 = sqrt(-1)*s1*besselh(nL(l)+1, 1, k1*R_AC).*exp(sqrt(-1)*(nL(l)+1)*Phi_AC);
        
        AL1 = [AL1 HjmL1];
        AL2 = [AL2 HjmL2];
        
        HjlL1 = besselh(nL(l), 1, k2*R_BC).*exp(sqrt(-1)*nL(l)*Phi_BC);
        HjlL2 = sqrt(-1)*s2*besselh(nL(l)+1, 1, k2*R_BC).*exp(sqrt(-1)*(nL(l)+1)*Phi_BC);
        
        BL1 = [BL1  HjlL1];
        BL2 = [BL2  HjlL2];
        
        HjmR1 = besselh(nL(l), 1, k1*R_AW).*exp(sqrt(-1)*nL(l)*Phi_AW);
        HjmR2 = sqrt(-1)*s1*besselh(nL(l)+1, 1, k1*R_AW).*exp(sqrt(-1)*(nL(l)+1)*Phi_AW);
        
        AR1 = [AR1 HjmR1];
        AR2 = [AR2 HjmR2];
        
        HjlR1 = besselh(nL(l), 1, k3*R_CW).*exp(sqrt(-1)*nL(l)*Phi_CW);
        HjlR2 = sqrt(-1)*s3*besselh(nL(l)+1, 1, k3*R_CW).*exp(sqrt(-1)*(nL(l)+1)*Phi_CW);
        
        CR1 = [CR1  HjlR1];
        CR2 = [CR2  HjlR2];
        
end

J1=NjC;L1=length(nL)*(NlC);
J2=NjW;L2=length(nL)*(NlW);
MA=[AL1;AL2;AR1;AR2];
MB=[BL1;BL2;zeros(J2,L1);zeros(J2,L1)];
MC=[zeros(J1,L2);zeros(J1,L2);CR1;CR2];

M=[MA -MB -MC];
C = pinv(M)*pln;
    
ee = norm(M*C - pln)/norm(pln);

N=500;
x_choose=linspace(-2,2,N);
y_choose=linspace(-2,1.5,N);

[xx,yy]=meshgrid(x_choose,y_choose);
zz=xx+sqrt(-1)*yy;

psiA=zeros(N,N);
psiB=zeros(N,N);

X2=real(ZjC);
Y2=imag(ZjC);

X3=real(ZjW);
Y3=imag(ZjW);

Int2=inpolygon(xx,yy,X2,Y2);
Int3=inpolygon(xx,yy,X3,Y3);
Out1=zeros(N,N);

for i=1:N
    for j=1:N
        if Int2(i,j)==0 && Int3(i,j)==0
            Out1(i,j)=1;
        end
    end
end

t=1;
for i=1:length(nL)
    
    for j=1:NmC
        
        Z_p=(zz-ZmC(j));
        r_p=abs(Z_p);
        theta_p=angle(Z_p);
        
        psiA=psiA+C(t)*besselh(nL(i), 1, k1*r_p).*exp(sqrt(-1)*nL(i)*theta_p).*Out1;
        psiB=psiB+s1*sqrt(-1)*C(t)*besselh(nL(i)+1, 1, k1*r_p).*exp(sqrt(-1)*(nL(i)+1)*theta_p).*Out1;
        
        t=t+1;
    end
    
	for j=1:NmW
        
        Z_p=(zz-ZmW(j));
        r_p=abs(Z_p);
        theta_p=angle(Z_p);
        
        psiA=psiA+C(t)*besselh(nL(i), 1, k1*r_p).*exp(sqrt(-1)*nL(i)*theta_p).*Out1;
        psiB=psiB+s1*sqrt(-1)*C(t)*besselh(nL(i)+1, 1, k1*r_p).*exp(sqrt(-1)*(nL(i)+1)*theta_p).*Out1;
        
        t=t+1;
    end
    
end

for i=1:length(nL)
    
    for j=1:NlC
        
        Z_p=(zz-ZlC(j));
        r_p=abs(Z_p);
        theta_p=angle(Z_p);
        
        psiA=psiA+C(t)*besselh(nL(i), 1, k2*r_p).*exp(sqrt(-1)*nL(i)*theta_p).*Int2;
        psiB=psiB+s2*sqrt(-1)*C(t)*besselh(nL(i)+1, 1, k2*r_p).*exp(sqrt(-1)*(nL(i)+1)*theta_p).*Int2;
        
        t=t+1;
    end
    
end

for i=1:length(nL)
    
    for j=1:NlW
        
        Z_p=(zz-ZlW(j));
        r_p=abs(Z_p);
        theta_p=angle(Z_p);
        
        psiA=psiA+C(t)*besselh(nL(i), 1, k3*r_p).*exp(sqrt(-1)*nL(i)*theta_p).*Int3;
        psiB=psiB+s3*sqrt(-1)*C(t)*besselh(nL(i)+1, 1, k3*r_p).*exp(sqrt(-1)*(nL(i)+1)*theta_p).*Int3;
        
        t=t+1;
    end
    
end

psiA=psiA+exp(sqrt(-1)*kx*xx).*cos(imag(alpha1)*yy).*Int3;
psiB=psiB+s3*exp(sqrt(-1)*kx*xx).*cos(imag(alpha1)*yy).*Int3;

psiA=psiA+A_norm*exp(sqrt(-1)*kx*xx).*exp(-alpha2*abs(yy+deviate_W)).*Int2;
psiB=psiB+A_norm*s3*exp(sqrt(-1)*kx*xx).*exp(-alpha2*abs(yy+deviate_W)).*Int2;

end

