clear all;
clc

%% Import boundary file
name = 'inputNONLINEAR.xlsx';
rb = xlsread(name,'A:A');               % R boundary
zb = xlsread(name,'B:B');               % Z boundary
bv = xlsread(name,'C:C');               % Psi in boundary

%% Reading matrix file
% adjusted by the number of tab from data on excel
G      = xlsread('G.xlsx','A1:CB80');			
Qtotb  = xlsread('Qtot.xlsx','A1:O80'); % Q in boundary
Qtot   = xlsread('Qtot.xlsx','A1:O980'); % Q in domain and boundary
GG     = xlsread('GG.xlsx','A1:CB900');
fprintf('Complete import boundary file\n');

%% initiation
errorpsi = 1;  %initiation error psi
iteration = 0; %initiation iteration
Level = 4;

%% determine the midpoint, element length, and normal direction for each element
n = length(rb)-1; %n    = number of element
for i=1:n
    rm(i)=0.5*(rb(i)+rb(i+1));
    zm(i)=0.5*(zb(i)+zb(i+1));
    lm(i)=sqrt((rb(i+1)-rb(i))^2+(zb(i+1)-zb(i))^2);
    nr(i)=(zb(i+1)-zb(i))/lm(i);
    nz(i)=(rb(i)-rb(i+1))/lm(i);
end

%% Generate grid in domain using TFI
% number of discretizations along xi and eta axis
mTFI = 30;
nTFI = 30;
% discretize along xi and eta axis
xiTFI = linspace(0.01,0.99,mTFI) ;
etaTFI = linspace(0.01,0.99,nTFI) ;
wTFI=0;
for i = 1:mTFI
    XiTFI = xiTFI(i) ;
    for j = 1:nTFI
        EtaTFI = etaTFI(j) ;
        wTFI=wTFI+1;
        % Transfinite Interpolation 
        XY = (1-EtaTFI)*Rb(XiTFI)+EtaTFI*Rt(XiTFI)+(1-XiTFI)*Rl(EtaTFI)+XiTFI*Rr(EtaTFI)......
            -(XiTFI*EtaTFI*Rt(1)+XiTFI*(1-EtaTFI)*Rb(1)+EtaTFI*(1-XiTFI)*Rt(0)+(1-XiTFI)*(1-EtaTFI)*Rb(0)) ;
        x(wTFI,1) = XY(1) ;
        y(wTFI,1) = XY(2) ;
    end
end

while (errorpsi>(1e-5))
    iteration = iteration + 1;
    fprintf('iteration -:%d \n',iteration);
    if (iteration == 1)
        % Review the boundary area
        for m=1:n
             O=0; %initiation for element matrix columns of the Q matrix by level
             for L=0:Level              %L
                 for M=0:Level          %M
                     if((L+M)<=Level)
                         O=O+1;
                         rlzm(m,O)=((rm(m)/1.32)^L)*(((zm(m))/0.544)^M);
                     end
                 end
             end
             psi1(m,1)=0; % psi initial
             psi2(m,1)=0; % psi after 
             ji(m,1)=f1(rm(m),psi1(m,1),1,0.7,4.78,-0.48);
        end
        fprintf('Finish to make element of matrix Q and G for boundary area\n');
        
        %% Review the domain area
        % form matriks G and Q for domain area,
        w=0;
        E = wTFI; %number of discretization in domain, becomes loop too
        fprintf('%d loop\n',E);
        for i=1:E		
            psi1(i+n,1)=0.6;
            ji(i+n,1)=f1(x(i),psi1(i+n,1),1,0.7,4.78,-0.48);
             O=0;
             for L=0:Level
                 for M=0:Level
                     if((L+M)<=Level)
                         O=O+1;
                         rlzm(i+n,O)=((x(i)/1.32)^L)*(((y(i))/0.544)^M);
                     end
                 end
             end
             fprintf('loop - %d\n',i);
        end
        fprintf('Finish to make matrix element G and Q part for domain area\n');
        [TU,TS,TV]=svd(G,0);
        alpha = rlzm\ji;
    end
    b=Qtotb*alpha; 
	%determine q using svd
    q=TV*((TU'*b)./diag(TS)); 
    %% menghitung matrix bagian dalam
    w=0;
    for i=1:wTFI
            u(i,1)=GG(i,:)*q+Qtot(i+n,:)*alpha;
            psi2(i+n,1)=u(i,1);
            ji(i+n,1)=f1(x(i),psi2(i+n,1),1,0.7,4.78,-0.48);
            jd(i,1)=ji(i+n,1)/(x(i)*1.2566);
    end
	%determine new alpha
    alpha2=rlzm\(ji);
    %%check convergence
    errorpsi=abs((psi2'-psi1')/psi1');
    errorspsi(iteration,1)=errorpsi;
    fprintf('loop ke %d, error psi: %8.7f\n', iteration,errorpsi);
        if(errorpsi>1e-5)
            psi1=psi2;
            alpha=alpha2;
        end
    if (iteration>100)
        break;
    end
end
xlin=linspace(min(x),max(x),300); 
ylin=linspace(min(y),max(y),300); 
[XX,YY]=meshgrid(xlin,ylin); 
 %analitics
 for i=1:wTFI
         u2(i,1)=(x(i)^4/8)+0.08285315656 + -0.1926498500*(x(i)^2) + -0.04905557541*(y(i)^2-x(i)^2*log(x(i)))+-0.04705357733*(x(i)^4-4*x(i)^2*y(i)^2)+...
              0.004892596401*(2*y(i)^4-9*x(i)^2*y(i)^2+3*x(i)^4*log(x(i))-12*x(i)^2*y(i)^2*log(x(i))) + -0.004219788580*(x(i)^6-12*x(i)^4*y(i)^2+8*x(i)^2*y(i)^4) + ...
             -0.0001088519924*(8*y(i)^6-140*x(i)^2*y(i)^4+75*x(i)^4*y(i)^2-15*x(i)^6*log(x(i))+180*x(i)^4*y(i)^2*log(x(i))-120*x(i)^2*y(i)^4*log(x(i)));
         jia= f1(x(i),u2(i,1),1,0.7,4.78,-0.48);
         jda(i,1)=jia/(x(i)*1.2566);
 end
fit_SAR1=griddata(x,y,u,XX,YY);  
fit_SAR2=griddata(x,y,-1.*u2,XX,YY);
%% Plotting
%psi
 figure(1)
 hold on; box on
 contour(XX,YY,fit_SAR1,20);
 contour(XX,YY,fit_SAR2,20,'LineStyle','--');
 plot(rb,zb,'LineStyle','-','Color','r','LineWidth',2)
 axis image
 colormap 'jet'
 title('Polodail Flux Analytics Vs BEM')
 xlabel('R/Ro') % x-axis label
 ylabel('Z/Zo') % y-axis label
 legend('BEM','Analytics','Location','southwest')
 colormap 'jet'
 
fit_SAR1j=griddata(x,y,jd,XX,YY);  
fit_SAR2j=griddata(x,y,jda,XX,YY);
 %current 
 figure(2)
 hold on; box on
 contour(XX,YY,fit_SAR1j,20);
 contour(XX,YY,fit_SAR2j,20,'LineStyle','--');
 plot(rb,zb,'LineStyle','-','Color','r','LineWidth',2)
 axis image
 colormap 'jet'
 title('Current Density Analytics Vs BEM')
 xlabel('R/Ro') % x-axis label
 ylabel('Z/Zo') % y-axis label
 legend('BEM','Analytics','Location','southwest')
 colormap 'jet'
 
 [iii,jjj] = find(fit_SAR1 == max(max(fit_SAR1)));
 %full colour Psi
 figure(3)
 hold on; box on
 pcolor(XX,YY,fit_SAR1);
 shading interp
 axis image
 colormap 'jet'
 title('Poloidal Flux BEM')
 xlabel('R/Ro') % x-axis label
 ylabel('Z/Zo') % y-axis label
 scatter(XX(iii,jjj),YY(iii,jjj),25,'MarkerEdgeColor','w','MarkerFaceColor','w','LineWidth',1)
 scatter(XX(iii,iii),YY(iii,iii),25,'MarkerEdgeColor','w','MarkerFaceColor','w','LineWidth',1)
 
 %full colour current
 figure(4)
 pcolor(XX,YY,fit_SAR1j);
 shading interp
 axis image
 colormap 'jet'
 title('Current Density BEM')
 xlabel('R/Ro') % x-axis label
 ylabel('Z/Zo') % y-axis label

 % Numerik only contour flux
 figure(5)
 hold on; box on
 contour(XX,YY,fit_SAR1,10);
 plot(rb,zb,'LineStyle','-','Color','r','LineWidth',2)
 axis image
 colormap 'jet'
 %title('Poloidal Flux BEM')
 xlabel('R/Ro') % x-axis label
 ylabel('Z/Zo') % y-axis label
 legend('BEM','Location','southwest')
 colormap 'jet'
 
 % Analytics only contour flux
 figure(6)
 hold on; box on
 contour(XX,YY,fit_SAR2,10,'LineStyle','--');
 plot(rb,zb,'LineStyle','-','Color','r','LineWidth',2)
 axis image
 colormap 'jet'
 %title('Poloidal Flux Analytics')
 xlabel('R/Ro') % x-axis label
 ylabel('Z/Zo') % y-axis label
 legend('Analytic','Location','southwest')
 colormap 'jet'
 