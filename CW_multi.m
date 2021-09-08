clear
clc

load('data_spectral.mat')

for i=1:79
    kx=kx_choose(i);
    E=E_choose(i);
    
    [xx yy ZjC ZjW ee psiA psiB] = Single_CW(kx,E);
    
    [xx0,yy0,Pd0,ee0] = Norm_W(kx,E);
    
    save(['data_',num2str(i),'.mat'])
    disp(i)
end
