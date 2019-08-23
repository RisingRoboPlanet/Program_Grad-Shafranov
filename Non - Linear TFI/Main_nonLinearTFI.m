clear all;
clc

%% Import boundary file
name = 'inputNONLINEAR.xlsx';
rb = xlsread(name,'A:A');               % R boundary
zb = xlsread(name,'B:B');               % Z boundary
bv = xlsread(name,'C:C');               % Psi in boundary
Lr = 1.32;
Lz = 0.544;

%% Reading matrix file
% adjusted by the number of tab from data on excel
G      = xlsread('G.xlsx','A1:CB80');			
Qtotb  = xlsread('Qtot.xlsx','A1:AB80'); % Q in boundary
Qtot   = xlsread('Qtot.xlsx','A1:AB480'); % Q in domain and boundary
GG     = xlsread('GG.xlsx','A1:CB400');
fprintf('Complete import boundary file\n');

%% initiation
errorpsi = 1;  %initiation error psi
iteration = 0; %initiation iteration
Level = 6;
errorlambda   = 1; % initiation lambda error
lambda1 = 1; %initiation lambda

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
mTFI = 20;
nTFI = 20;
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

while (errorpsi>(1e-5) && errorlambda >(1e-5))
    iteration = iteration + 1;
    fprintf('iteration -:%d \n',iteration);
    if (iteration == 1)
        % Review the boundary area
        for m=1:n
             IT(1,m)=findq(nz(m),lm(m),rb(m),zb(m+1));
             O=0; %initiation for element matrix columns of the Q matrix by level
             for L=0:Level              %L
                 for M=0:Level          %M
                     if((L+M)<=Level)
                         O=O+1;
                         rlzm(m,O)=((rm(m)/Lr)^L)*(((zm(m))/Lz)^M);
                     end
                 end
             end
             psi1(m,1)=0; % psi initial
             psi2(m,1)=0; % psi after 
             ji(m,1)=f1(rm(m),psi1(m,1),1,0.7,4.78,-0.48);
        end
       O=0;
       for L=0:Level
           for M=0:Level
               if((L+M)<=Level)
                   O=O+1;
                   PPP=0;
                   for j=1:n
                       Ptott=findddeln(rb(j),zb(j),rb(j+1),zb(j+1),lm(j),nr(j),nz(j),Lr,Lz,0,L,M);
                       PPP=PPP+Ptott;
                   end
                   V(1,O)=PPP;
               end
           end
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
                         rlzm(i+n,O)=((x(i)/Lr)^L)*(((y(i))/Lz)^M);
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
    for i=1:wTFI
            u(i,1)=GG(i,:)*q+Qtot(i+n,:)*alpha;
    end
    uMAX = max(max(abs(u)));
    u=u./uMAX;
    for i=1:wTFI
            psi2(i+n,1)=u(i,1);
            ji(i+n,1)=f1(x(i),psi2(i+n,1),1,0.7,4.78,-0.48);
            jd(i,1)=ji(i+n,1)/(x(i));
    end
	%determine new alpha
    alpha2=rlzm\(lambda1*ji);
    %find new lambda
    lambda2=lambda1*((V*alpha)/(V*alpha2))/uMAX;
    %%check convergence
    errorpsi=abs((psi2'-psi1')/psi1');
    errorlambda=abs((lambda2-lambda1)/lambda1);
    errorspsi(iteration,1)=errorpsi;
    errorslambda(iteration,1)=errorlambda;
    fprintf('loop: %d, lambda error: %8.7f, psi error: %8.7f\n', iteration,errorlambda,errorpsi);
        if(errorpsi>1e-5)
            psi1=psi2;
            lambda1=lambda2;
            alpha=alpha2;
        end
    if (iteration>100)
        break;
    end
end
xlin=linspace(min(x),max(x),300); 
ylin=linspace(min(y),max(y),300); 
[XX,YY]=meshgrid(xlin,ylin); 
fit_SAR1=griddata(x,y,u,XX,YY);  

% %% Plotting
% %psi
%  figure(1)
%  hold on; box on
%  contour(XX,YY,fit_SAR1,10);
%  axis image
%  colormap 'jet'
%  title('Poloidal Flux BEM')
%  xlabel('R/Ro') % x-axis label
%  ylabel('Z/Zo') % y-axis label
%  colormap 'jet'
%  colorbar

[i,jj]=find(fit_SAR1==max(max(fit_SAR1)));
for JJJ = 1:300
    XZ(JJJ,1) = XX(150,JJJ);
    psiR(JJJ,1) = fit_SAR1(150,JJJ);
end
 
fit_SAR1j=griddata(x,y,jd,XX,YY);  
%  %current 
%  figure(2)
%  hold on; box on
%  contour(XX,YY,fit_SAR1j,10);
%  axis image
%  colormap 'jet'
%  title('Current Density BEM')
%  xlabel('R/Ro') % x-axis label
%  ylabel('Z/Zo') % y-axis label
%  colormap 'jet'
%  colorbar
for JJJ = 1:300
    currentR(JJJ,1) = fit_SAR1j(150,JJJ);
end
%  
%  %full colour Psi
%   figure(3)
%   pcolor(XX,YY,fit_SAR1);
%   shading interp
%   axis image
%   colormap 'jet'
%   title('Poloidal Flux BEM')
%   xlabel('R/Ro') % x-axis label
%   ylabel('Z/Zo') % y-axis label
%   colorbar
%   
% %  %full colour current
%   figure(4)
%   pcolor(XX,YY,fit_SAR1j);
%   shading interp
%   axis image
%   colormap 'jet'
%   title('Current Density BEM')
%   xlabel('R/Ro') % x-axis label
%   ylabel('Z/Zo') % y-axis label
%   colorbar
 
 max(max(jd))
 [iii,jjj] = find(fit_SAR1 == max(max(fit_SAR1)))
 XX(iii,jjj)
 [iiii,jjjj] = find(fit_SAR1j == max(max(fit_SAR1j)))
 XX(iiii,jjjj)