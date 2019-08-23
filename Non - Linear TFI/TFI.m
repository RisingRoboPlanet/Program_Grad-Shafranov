% To demonstarte grid generation using Transfinite Interpolation (TFI)
%{
 Author : Siva Srinivas Kolukula                                
          Senior Research Fellow                                
          Structural Mechanics Laboratory                       
          Indira Gandhi Center for Atomic Research              
          India                                                 
 E-mail : allwayzitzme@gmail.com                                         
          http://sites.google.com/site/kolukulasivasrinivas/                 
%}

% Reference: Fundametnals of Grid Generation - Knupp, Steinberg

clc 
clear all ;

% number of discretizations along xi and eta axis

m = 7 ;
n = 7 ;

% discretize along xi and eta axis
xi = linspace(0.01,0.99,m) ;
eta = linspace(0.01,0.99,n) ;

% Initialize matrices in x and y axis
% X = zeros(m,n) ;
% Y = zeros(m,n) ;

w=0;
for i = 1:m
    Xi = xi(i) ;
    for j = 1:n
        Eta = eta(j) ;
        w=w+1;
        % Transfinite Interpolation 
        XY = (1-Eta)*Rb(Xi)+Eta*Rt(Xi)+(1-Xi)*Rl(Eta)+Xi*Rr(Eta)......
            -(Xi*Eta*Rt(1)+Xi*(1-Eta)*Rb(1)+Eta*(1-Xi)*Rt(0)+(1-Xi)*(1-Eta)*Rb(0)) ;
    
        X(w,1) = XY(1) ;
        Y(w,1) = XY(2) ;
        
    end
end

%plotgrid(X,Y) ;
