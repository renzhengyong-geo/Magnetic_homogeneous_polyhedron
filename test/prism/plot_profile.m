clear all;
close all;
%% load closed-form from prism
xyz    = load('xyz.dat');
x      = xyz(:,1);
y      = xyz(:,2);
%% solution
T      = load('T.dat');                  %% new solutions
T2     = load('old_m3d/T1.dat');         %% previous solutions
T3     = load('health2005/result.dat');  %% health's solutions

lsize=2;
fsize=12;
msize= 4;

%% (xx,xy,xz,yx,yy,yz,zx,zy,zz) in T.dat
%% plot tensor (xx,yy,zz)

figure('Position',[0 0 850 850]); 
subplot(3,2,1);
plot(y, T(:,1),'k-', y, T2(:,1), 'bo', y, T3(:,1), 'rh', 'MarkerSize',  msize, 'LineWidth', lsize);
xlim([-25,25]);
ylim([-5e4,1e4]);
ylabel('T_{xx}(nT/m)','FontSize', fsize);
set(gca,'fontsize',fsize);

subplot(3,2,2);
plot(y, (T(:,1)-T2(:,1))./T2(:,1)*100, 'bo', y, (T(:,1)-T3(:,1))./T3(:,1)*100, 'rh', 'MarkerSize',  msize, 'LineWidth', lsize);
xlim([-25,25]);
ylabel('Relative error of T_{xx} (%)','FontSize', fsize);
set(gca,'fontsize',fsize);

subplot(3,2,3);
plot(y, T(:,5), 'k-', y, T2(:,5), 'bo',  y, T3(:,5), 'rh', 'MarkerSize',  msize, 'LineWidth', lsize);
xlim([-25,25]);
ylim([-5000,7000]);
ylabel({'T_{yy}(nT/m)'},'FontSize', fsize);
set(gca,'fontsize',fsize);

subplot(3,2,4);
plot(y, (T(:,5)-T2(:,5))./T2(:,5)*100, 'bo', y, (T(:,5)-T3(:,5))./T3(:,5)*100, 'rh', 'MarkerSize',  msize, 'LineWidth', lsize);
xlim([-25,25]);
ylabel('Relative error of T_{yy} (%)','FontSize', fsize);
set(gca,'fontsize',fsize);

subplot(3,2,5);
plot(y, T(:,9), 'k-',  y, T2(:,9), 'bo', y, T3(:,9), 'rh', 'MarkerSize',  msize, 'LineWidth', lsize);
xlim([-25,25]);
xlabel('y(m)','FontSize', fsize);
ylabel({'T_{zz}(nT/m)'},'FontSize', fsize);
ylim([0,5e4]);
hl=legend('Our new solution', 'Ren et al.(2017)''s result', 'Heath et al.(2005)''s result' );
set(hl, 'Box', 'off', 'location', 'Best','FontSize', fsize) ;
set(gca,'fontsize',fsize);

subplot(3,2,6);
plot(y, (T(:,9)-T2(:,9))./T2(:,9)*100, 'bo', y, (T(:,9)-T3(:,9))./T3(:,9)*100, 'rh','MarkerSize',  msize, 'LineWidth', lsize);
xlim([-25,25]);
xlabel('y(m)','FontSize', fsize);
ylabel('Relative error of T_{zz} (%)','FontSize', fsize);
set(gca,'fontsize',fsize);



%% plot tensor (xy,xz,yz)
%% plot tensor (xy)
figure('Position',[0 0 850 850]); 
subplot(3,2,1);
plot(y, T(:,2),'k-', y, T2(:,2), 'bo', y, T3(:,2), 'rh', 'MarkerSize',   msize, 'LineWidth', lsize);
xlim([-25,25]);
ylabel({'T_{xy}(nT/m)'},'FontSize', fsize);
%ylim([-1600,1200]);
set(hl, 'Box', 'off', 'location', 'Best','FontSize', fsize) ;
xlabel('y(m)','FontSize', fsize);
set(gca,'fontsize',fsize);

subplot(3,2,2);
T2(11,2) = 0; 
plot(y, (T(:,2)-T2(:,2))./T2(:,2)*100, 'bo', y, (T(:,2)-T3(:,2))./T2(:,3)*100, 'rh', 'MarkerSize',  msize, 'LineWidth', lsize);
xlim([-25,25]);
%ylim([-8e-13,4e-13]);
ylabel('Relative error of T_{xy} (%)','FontSize', fsize);
xlabel('y(m)','FontSize', fsize);
set(gca,'fontsize', fsize);

%% plot tensor (xz)
subplot(3,2,3);
plot(y, T(:,3),'k-', y, T2(:,3), 'bo', y, T3(:,3), 'rh', 'MarkerSize',  msize, 'LineWidth', lsize);
xlim([-25,25]);
ylabel({'T_{xz}(nT/m)'},'FontSize', fsize);
ylim([0,1400]);
set(hl, 'Box', 'off', 'location', 'Best','FontSize', fsize) ;
xlabel('y(m)','FontSize', fsize);
set(gca,'fontsize',fsize);

subplot(3,2,4);
%T2(11,8) = 0; 
plot(y, (T(:,3)-T2(:,3))./T2(:,3)*100, 'bo', y, (T(:,3)-T3(:,3))./T3(:,3)*100, 'rh','MarkerSize',  msize, 'LineWidth', lsize);
xlim([-25,25]);
%ylim([-8e-13,4e-13]);
ylabel('Relative error of T_{xz} (%)','FontSize', fsize);
xlabel('y(m)','FontSize', fsize);
set(gca,'fontsize', fsize);

%% plot tensor (yz)
subplot(3,2,5);
plot(y, T(:,8),'k-', y, T2(:,8), 'bo', y, T3(:,8), 'rh','MarkerSize',  msize, 'LineWidth', lsize);
xlim([-25,25]);
ylabel({'T_{yz}(nT/m)'},'FontSize', fsize);
ylim([-1800,1200]);
hl=legend('Our new solution', 'Ren et al.(2017)''s result', 'Heath et al.(2005)''s result' );
set(hl, 'Box', 'off', 'location', 'Best','FontSize', fsize) ;
xlabel('y(m)','FontSize', fsize);
set(gca,'fontsize',fsize);

subplot(3,2,6);
%T2(11,8) = 0; 
plot(y, (T(:,8)-T2(:,8))./T2(:,8)*100, 'bo', y, (T(:,8)-T3(:,8))./T3(:,8)*100, 'rh', 'MarkerSize',  msize, 'LineWidth', lsize);
xlim([-25,25]);
%ylim([-8e-13,4e-13]);
ylabel('Relative error of T_{yz} (%)','FontSize', fsize);
xlabel('y(m)','FontSize', fsize);
set(gca,'fontsize', fsize);



