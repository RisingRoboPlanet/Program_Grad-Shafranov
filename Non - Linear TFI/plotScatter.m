%% Program grid
clear all;
clc

%% Import boundary file
name = 'inputNONLINEAR.xlsx';
rb = xlsread(name,'A:A');               % R boundary
zb = xlsread(name,'B:B');               % Z boundary
bv = xlsread(name,'C:C');               % Psi in boundary

%% Generate grid in domain using TFI
% number of discretizations along xi and eta axis
mTFI2 = 20;
nTFI2 = 20;
% discretize along xi and eta axis
xiTFI = linspace(0,1,mTFI2) ;
etaTFI = linspace(0,1,nTFI2) ;
wTFI=0;
for i = 1:mTFI2
    XiTFI = xiTFI(i) ;
    for j = 1:nTFI2
        EtaTFI = etaTFI(j) ;
        wTFI=wTFI+1;
        % Transfinite Interpolation 
        XY = (1-EtaTFI)*Rb(XiTFI)+EtaTFI*Rt(XiTFI)+(1-XiTFI)*Rl(EtaTFI)+XiTFI*Rr(EtaTFI)......
            -(XiTFI*EtaTFI*Rt(1)+XiTFI*(1-EtaTFI)*Rb(1)+EtaTFI*(1-XiTFI)*Rt(0)+(1-XiTFI)*(1-EtaTFI)*Rb(0)) ;
        x(wTFI,1) = XY(1) ;
        y(wTFI,1) = XY(2) ;
    end
end

%% Plotting
%psi
 figure(1)
 hold on; box on
 scatter(x,y,10)
 axis image
 plot(rb,zb,'LineStyle','-','Color','r','LineWidth',2)
 colormap 'jet'
 title('400 Nodes in domain')
 xlabel('R/Ro') % x-axis label
 ylabel('Z/Zo') % y-axis label
 
 
%% Generate grid in domain using TFI
% number of discretizations along xi and eta axis
mTFI1 = 10;
nTFI1 = 10;
% discretize along xi and eta axis
xiTFI = linspace(0,1,mTFI1) ;
etaTFI = linspace(0,1,nTFI1) ;
wTFI=0;
for i = 1:mTFI1
    XiTFI = xiTFI(i) ;
    for j = 1:nTFI1
        EtaTFI = etaTFI(j) ;
        wTFI=wTFI+1;
        % Transfinite Interpolation 
        XY = (1-EtaTFI)*Rb(XiTFI)+EtaTFI*Rt(XiTFI)+(1-XiTFI)*Rl(EtaTFI)+XiTFI*Rr(EtaTFI)......
            -(XiTFI*EtaTFI*Rt(1)+XiTFI*(1-EtaTFI)*Rb(1)+EtaTFI*(1-XiTFI)*Rt(0)+(1-XiTFI)*(1-EtaTFI)*Rb(0)) ;
        x(wTFI,1) = XY(1) ;
        y(wTFI,1) = XY(2) ;
    end
end
 figure(2)
 hold on; box on
 scatter(x,y,10)
 axis image
 colormap 'jet'
  plot(rb,zb,'LineStyle','-','Color','r','LineWidth',2)
 title('100 Nodes in domain')
 xlabel('R/Ro') % x-axis label
 ylabel('Z/Zo') % y-axis label
 
 %% Generate grid in domain using TFI
% number of discretizations along xi and eta axis
mTFI3 = 30;
nTFI3 = 30;
% discretize along xi and eta axis
xiTFI = linspace(0,1,mTFI3) ;
etaTFI = linspace(0,1,nTFI3) ;
wTFI=0;
for i = 1:mTFI3
    XiTFI = xiTFI(i) ;
    for j = 1:nTFI3
        EtaTFI = etaTFI(j) ;
        wTFI=wTFI+1;
        % Transfinite Interpolation 
        XY = (1-EtaTFI)*Rb(XiTFI)+EtaTFI*Rt(XiTFI)+(1-XiTFI)*Rl(EtaTFI)+XiTFI*Rr(EtaTFI)......
            -(XiTFI*EtaTFI*Rt(1)+XiTFI*(1-EtaTFI)*Rb(1)+EtaTFI*(1-XiTFI)*Rt(0)+(1-XiTFI)*(1-EtaTFI)*Rb(0)) ;
        x(wTFI,1) = XY(1) ;
        y(wTFI,1) = XY(2) ;
    end
end
 figure(3)
 hold on; box on
 scatter(x,y,10)
 axis image
 colormap 'jet'
  plot(rb,zb,'LineStyle','-','Color','r','LineWidth',2)
 title('900 Nodes in domain')
 xlabel('R/Ro') % x-axis label
 ylabel('Z/Zo') % y-axis label