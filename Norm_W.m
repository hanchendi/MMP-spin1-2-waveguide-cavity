function [xx,yy,Pd,ee] = Norm_W(kx,E)


%% Eigenmode for waveguide
epsilon=E;
a=0.1;
V=8;

alpha1=sqrt(kx^2-(epsilon-V)^2);
alpha2=sqrt(kx^2-epsilon^2);
A11=exp(alpha1*a);
A12=exp(-alpha1*a);
A21=1./(V-epsilon+kx).*alpha1.*exp(alpha1*a);
A22=1./(V-epsilon+kx).*(-alpha1).*exp(-alpha1*a);
B1=exp(-alpha2*a);
B2=1./(kx-epsilon).*(-alpha2).*exp(-alpha2*a);
T=B1./(A11+A12)-B2./(A21+A22);
A=B1./(A11+A12);


%% MMP

k=epsilon;
q=epsilon-V;

b=5;

Nl=4*204;
delta_l=0.02;
t1=floor(Nl/2*a/b);

xl(1:Nl/2)=linspace(-b-delta_l,b+delta_l,Nl/2);
yl(1:Nl/2)=(-a-delta_l)*ones(1,Nl/2);

xl=[xl (b+delta_l)*ones(1,t1)];
yl=[yl linspace(-a-delta_l,a+delta_l,t1)];

xl=[xl linspace(b+delta_l,-b-delta_l,Nl/2)];
yl=[yl (a+delta_l)*ones(1,Nl/2)];

xl=[xl (-b-delta_l)*ones(1,t1)];
yl=[yl linspace(a+delta_l,-a-delta_l,t1)];
xZl=xl+yl*sqrt(-1);
Nl=length(xl);

Nm=4*200;

delta_m=0.02;
t2=floor(Nm/2*a/b);
theta_l=linspace(0,2*pi,Nl);
xm(1:Nm/2)=linspace(-b+delta_m,b-delta_m,Nm/2);
ym(1:Nm/2)=(-a+delta_m)*ones(1,Nm/2);

xm=[xm (b-delta_m)*ones(1,t2)];
ym=[ym linspace(-a+delta_m,a-delta_m,t2)];

xm=[xm linspace(b-delta_m,-b+delta_m,Nm/2)];
ym=[ym (a-delta_m)*ones(1,Nm/2)];

xm=[xm (-b+delta_m)*ones(1,t2)];
ym=[ym linspace(a-delta_m,-a+delta_m,t2)];
xZm=xm+ym*sqrt(-1);
Nm=length(xm);

Nj=(Nl+Nm)*5;
t3=floor(Nj/2*a/b);
X(1:Nj/2)=linspace(-b,b,Nj/2);
Y(1:Nj/2)=(-a)*ones(1,Nj/2);

X=[ X (b)*ones(1,t3)];
Y=[ Y linspace(-a,a,t3)];

X=[ X linspace(b,-b,Nj/2)];
Y=[Y (a)*ones(1,Nj/2)];

X=[ X (-b)*ones(1,t3)];
Y=[ Y linspace(a,-a,t3)];

Z=X+Y*sqrt(-1);
Nj=length(X);

[zl, zjl] = meshgrid(xZl, Z);
Djl = zjl - zl; 
Phjl = angle(Djl);
Rjl = abs(Djl);

[zm, zjm] = meshgrid(xZm, Z);
Djm = zjm - zm;
Phjm = angle(Djm);
Rjm = abs(Djm);

%% MMP

L = 1; % the maximum angular momentum
nL = -L:L;

%plnA = exp(sqrt(-1)*kx*X).*exp(-alpha2*abs(Y));
%plnB = sign(q)*exp(sqrt(-1)*kx*X).*exp(-alpha2*abs(Y));

plnA = exp(sqrt(-1)*kx*X).*cos(imag(alpha1)*Y);
plnB = sign(q)*exp(sqrt(-1)*kx*X).*cos(imag(alpha1)*Y);

pln = conj([plnA plnB]');

Ai = []; Bi = [];
Ao = []; Bo = [];



for l = 1:length(nL)

        HjmA = besselh(nL(l), 1, abs(k)*Rjm).*exp(sqrt(-1)*nL(l)*Phjm);
        HjmB = sqrt(-1)*sign(k)*besselh(nL(l)+1, 1, abs(k)*Rjm).*exp(sqrt(-1)*(nL(l)+1)*Phjm);
        
        Ai = [Ai HjmA];
        Bi = [Bi HjmB];
        
        HjlA =besselh(nL(l), 1, abs(q)*Rjl).*exp(sqrt(-1)*nL(l)*Phjl);
        HjlB = sqrt(-1)*sign(q)*besselh(nL(l)+1, 1, abs(q)*Rjl).*exp(sqrt(-1)*(nL(l)+1)*Phjl);
        
        Ao = [Ao HjlA];
        Bo = [Bo HjlB];
        
end

ABi = [Ai; Bi];  
ABo = [Ao; Bo];
M = [ABi, -ABo];

C = pinv(M)*pln;
    
ee = norm(M*C - pln)/norm(pln);

N=200;
x_choose=linspace(-2,2,N);
y_choose=linspace(-0.1,0.1,N);

[xx,yy]=meshgrid(x_choose,y_choose);
zz=xx+sqrt(-1)*yy;

psiA=zeros(N,N);
psiB=zeros(N,N);

psi1=zeros(N,N);
psi2=zeros(N,N);
Ind=inpolygon(xx,yy,X,Y);
Oud=zeros(N,N);

for i=1:N
    for j=1:N
        if Ind(i,j)==0
            Oud(i,j)=1;
        end
    end
end

t=length(nL)*Nm+1;

for i=1:length(nL)

    
    for j=1:Nl
        
        Z_p=(zz-xZl(j));
        r_p=abs(Z_p);
        theta_p=angle(Z_p);
        
        psiA=psiA+C(t)*besselh(nL(i), 1, abs(q)*r_p).*exp(sqrt(-1)*nL(i)*theta_p).*Ind;
        psiB=psiB+sign(q)*sqrt(-1)*C(t)*besselh(nL(i)+1, 1, abs(q)*r_p).*exp(sqrt(-1)*(nL(i)+1)*theta_p).*Ind;
        
        t=t+1;
    end
end


psiA=psiA+exp(sqrt(-1)*kx*xx).*cos(imag(alpha1)*yy).*Ind;
psiB=psiB+sign(q)*exp(sqrt(-1)*kx*xx).*cos(imag(alpha1)*yy).*Ind;

Pd=abs(psiA).^2+abs(psiB).^2;


end

