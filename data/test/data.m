clear all
close all
clc

% L = 2*pi/3.0619; N = 10000;
% eps0=1; wp = 1;
% qe = -wp^2/(N/L); me = -qe;

spec = importdata('record_readme.out');
Nt = spec(1); Nmod = spec(2); dt = spec(3); Nx = spec(4); Nv = spec(5); Lx = spec(6); Lv = spec(7);
Nr = floor(Nt/Nmod);

fileID = fopen('xg.bin');
xg = fread(fileID,Nx,'double');
fileID = fopen('vg.bin');
vg = fread(fileID,2*Nv+1,'double');
[V,X] = meshgrid(vg,xg);
dx = Lx/Nx; dv = Lv/Nv;

fileID = fopen('E.bin');
E = fread(fileID,Nx*Nr,'double');
E = reshape(E,[Nx,Nr]);
fileID = fopen('rho.bin');
rho = fread(fileID,Nx*Nr,'double');
rho = reshape(rho,[Nx,Nr]);
fileID = fopen('phi.bin');
phi = fread(fileID,Nx*Nr,'double');
phi = reshape(phi,[Nx,Nr]);

fileID = fopen('KE.bin');
KE = fread(fileID,Nr,'double');
fileID = fopen('PE.bin');
PE = fread(fileID,Nr,'double');
fileID = fopen('TE.bin');
TE = fread(fileID,Nr,'double');

time = dt*Nmod*(1:Nr);
figure(1)
plot(time,KE,'-r',time,PE,'-b',time,TE,'-k');

fileID = fopen(strcat('src.bin'));
src = fread(fileID,Nx*(2*Nv+1),'double');
src = reshape(src,[Nx,(2*Nv+1)]);
figure(2)
surf(X,V,src,'LineStyle','none');
colorbar;
axis([0 Lx -Lv Lv]);

fileID = fopen(strcat('f0.bin'));
f0 = fread(fileID,Nx*(2*Nv+1),'double');
f0 = reshape(f0,[Nx,(2*Nv+1)]);
figure(3)
surf(X,V,f0,'LineStyle','none');
colorbar;
axis([0 Lx -Lv Lv]);

%%
close all

for i=1:Nr
    fileID = fopen(strcat('f/',num2str(i),'.bin'));
    f = fread(fileID,Nx*(2*Nv+1),'double');
    f = reshape(f,[Nx,(2*Nv+1)]);
    
    figure(1)
    surf(X,V,f,'LineStyle','none');
    colorbar;
    axis([0 Lx -Lv Lv]);
    xlabel('$X$','Interpreter','Latex');
    ylabel('$V$','Interpreter','Latex');
    title('Plasma distribution','Interpreter','Latex');
    set(gca,'fontsize',25,'ticklabelinterpreter','latex');
    view([0 0 1]);
    
    figure(2)
    surf(X,V,f0,'LineStyle','none');
    colorbar;
    axis([0 Lx -Lv Lv]);
    xlabel('$X$','Interpreter','Latex');
    ylabel('$V$','Interpreter','Latex');
    title('Initial distribution','Interpreter','Latex');
    set(gca,'fontsize',25,'ticklabelinterpreter','latex');
    view([0 0 1]);

%     %videoclip
%     frame = getframe(gcf);
%     writeVideo(writerObj,frame);
    
    fclose('all');
    
    pause();

end

% % videoclip close
% close(writerObj);
%

%%

t = linspace(0,Nt*dt,Nr+1); t = t(2:Nr+1);
figure(2)
semilogy(k,histSpectrum(:,1));
title(strcat('FFT for E in space at T=',num2str(dt)),'fontsize',20)
xlabel('k','fontsize',20);
ylabel('$\hat{E}$','fontsize',20,'Interpreter','Latex');
% 
figure(3)
plot(t,PE,'-b',t,KE,'-r',t,TE,'-k');
% axis([0 400 0 0.045]);
title('Energy history','fontsize',20);
xlabel('Time','fontsize',20);
ylabel('Energy','fontsize',20);
legend('PE','KE','TE');
% 
% figure(4)
% plot(xg,E(:,1));
% title(strcat('E field at T=',num2str(dt)),'fontsize',20);
% xlabel('x','fontsize',20);
% ylabel('y','fontsize',20);
% 
figure(5)
semilogy(t,histSpectrum(2,:),linspace(dt,40,1000),0.0001*exp((linspace(dt,40,1000))/2/sqrt(2)),'-r');
% axis([0 80 1e-6 1e8]);
title('Growth of most unstable mode $k$','fontsize',20,'Interpreter','Latex');
xlabel('Time','fontsize',20);
ylabel('$\hat{E}$','fontsize',20,'Interpreter','Latex');
h=legend('$\hat{E}(k)$','$e^{\omega t}$');
set(h,'Interpreter','Latex');