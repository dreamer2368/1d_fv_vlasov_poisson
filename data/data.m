clear all
close all
clc

L = 2*pi/3.0619; N = 10000;
eps0=1; wp = 1;
qe = -wp^2/(N/L); me = -qe;

spec = importdata('record_readme.out');
Nt = spec(1); Nr = spec(2); Nx = spec(3); Nv = spec(4); dx = spec(5); dv = spec(6); dt = spec(7);

fileID = fopen('E.bin');
E = fread(fileID,Nx*Nr,'double');
E = reshape(E,[Nx,Nr]);
fileID = fopen('rho.bin');
rho = fread(fileID,Nx*Nr,'double');
rho = reshape(rho,[Nx,Nr]);
fileID = fopen('phi.bin');
phi = fread(fileID,Nx*Nr,'double');
phi = reshape(phi,[Nx,Nr]);
fileID = fopen('xg.bin');
xg = fread(fileID,Nx,'double');
fileID = fopen('vg.bin');
vg = fread(fileID,2*Nv+1,'double');
fileID = fopen('T.bin');
T = fread(fileID,Nr,'double');
fileID = fopen('f.bin');
f = fread(fileID,Nx*(2*Nv+1)*Nr,'double');
f = reshape(f,[Nx*(2*Nv+1), Nr]);

fileID = fopen('KE.bin');
KE = fread(fileID,Nr,'double');
fileID = fopen('PE.bin');
PE = fread(fileID,Nr,'double');
fileID = fopen('TE.bin');
TE = fread(fileID,Nr,'double');

%%
close all

[V,X] = meshgrid(vg,xg);

U = reshape(f(:,1),[Nx,2*Nv+1]);
U = dx*sum( U, 1 ); maxU = max(U);

% %video clip
% writerObj = VideoWriter('limiter.avi');
% writerObj.FrameRate = 40;
% open(writerObj);

histSpectrum = [];
NFFT = 2^nextpow2(Nx);
k = 2*pi*(1/(2*dx))*linspace(0,1,NFFT/2+1);

figure(1)
% mesh(X,V,reshape(f(:,Nr),[Nx,2*Nv+1]) );
% view([0 0 -1]);
mesh(X,V,reshape(f(:,1230),[Nx,2*Nv+1]) );
axis([0 L -1 1 0 1.5*max(max(f))]);
xlabel('$X$','Interpreter','Latex');
ylabel('$V$','Interpreter','Latex');
title('Two-stream Instability (Finite-Volume)','Interpreter','Latex', ...
        'fontsize',25);
set(gca,'fontsize',30,'TickLabelInterpreter', 'latex');
view([0 0 -1]);


%%

for i=1:Nr
%     mesh(X,V,reshape(f(:,i),[Nx,2*Nv+1]) );
%     axis([0 L -1 1 0 1.5*max(max(f))]);
%     xlabel('$X$','Interpreter','Latex');
%     ylabel('$V$','Interpreter','Latex');
%     title('Plasma distribution','Interpreter','Latex');
%     set(gca,'fontsize',25);
%     view([0 0 -1]);
%     U = reshape(f(:,i),[Nx,2*Nv+1]);
%     U = dx*sum( U, 1 );
%     plot( vg, U );
%     axis([min(vg) max(vg) 0 maxU]);
%     plot(xg,E(:,i));
%     axis([0 L -0.2 0.2]);
%     pause();

%     %videoclip
%     frame = getframe(gcf);
%     writeVideo(writerObj,frame);

    Y = fft(E(:,i),NFFT)/length(E(:,i));
    histSpectrum = [histSpectrum 2*abs(Y(1:NFFT/2+1))];

end

% % videoclip close
% close(writerObj);

t = linspace(0,Nt*dt,Nr+1); t = t(2:Nr+1);
figure(2)
semilogy(k,histSpectrum(:,1));
title(strcat('FFT for E in space at T=',num2str(dt)),'fontsize',20)
xlabel('k','fontsize',20);
ylabel('$\hat{E}$','fontsize',20,'Interpreter','Latex');

figure(3)
plot(t,PE,'-b',t,KE,'-r',t,TE,'-k');
% axis([0 400 0 0.045]);
title('Energy history','fontsize',20);
xlabel('Time','fontsize',20);
ylabel('Energy','fontsize',20);
legend('PE','KE','TE');

figure(4)
plot(xg,E(:,1));
title(strcat('E field at T=',num2str(dt)),'fontsize',20);
xlabel('x','fontsize',20);
ylabel('y','fontsize',20);

figure(5)
semilogy(t,histSpectrum(2,:),linspace(dt,40,1000),0.0001*exp((linspace(dt,40,1000))/2/sqrt(2)),'-r');
% axis([0 80 1e-6 1e8]);
title('Growth of most unstable mode $k$','fontsize',20,'Interpreter','Latex');
xlabel('Time','fontsize',20);
ylabel('$\hat{E}$','fontsize',20,'Interpreter','Latex');
h=legend('$\hat{E}(k)$','$e^{\omega t}$');
set(h,'Interpreter','Latex');