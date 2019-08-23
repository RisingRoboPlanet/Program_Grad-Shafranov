% Program for creating matrix G and Q
% Using boundary input file in excel
clear all;
clc

%% Import boundary file
name = 'inputNONLINEAR.xlsx';
rb = xlsread(name,'A:A');               % R boundary
zb = xlsread(name,'B:B');               % Z boundary
bv = xlsread(name,'C:C');               % Psi in boundary
fprintf('Complete import boundary file\n');
Lr = 1.32;
Lz = 0.544;
%% initiation
Level = 2;     %level of expansion current profile

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
mTFI = 30 ;
nTFI = 30 ;
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

for m=1:n
	 for j=1:n
		 %calculate element of matrix G for boundary
		 Gmj=findg(rm(m),zm(m),nr(j),nz(j),lm(j),rb(j),zb(j),rb(j+1),zb(j+1));
		 G(m,j)=-Gmj;
	 end
	 O=0; %initiation for element matrix columns of the Q matrix by level
	 for L=0:Level              %L
		 for M=0:Level          %M
			 if((L+M)<=Level)
				 O=O+1;
				 rlzm(m,O)=((rm(m)/Lr)^L)*(((zm(m))/Lz)^M);
				 Ptot=0;     %calculates the sum of the total integral boundary on the 2nd part of the right hand side of the Q parameter
				 for j=1:n   %element boundary
					 P=findQ3(rm(m),zm(m),nr(j),nz(j),lm(j),rb(j),zb(j),rb(j+1),zb(j+1),L,M,0,Lz,Lr);
					 Ptot=Ptot+P;
				 end
				 Q(m,O)=Q2(rm(m),zm(m),L,M,0,Lz,Lr,0.5)-Ptot; %form matriks Q for boundary area
				 Qtot(m,O)=Q(m,O);
			 end
		 end
	 end
end

%% Review the domain area
% form matriks G and Q for domain area,
w=0;
E = wTFI; %number of discretization in domain, becomes loop too
fprintf('%d loop\n',E);
for i=1:E
         for k=1:n
             GG(i,k)=findg(x(i),y(i),nr(k),nz(k),lm(k),rb(k),zb(k),rb(k+1),zb(k+1)); %% membentuk matriks G untuk bagian dalam
         end
         O=0; %inisiasi untuk matriks kolom elemen dari matrikx part = FF (part dari Q bagian dalam)
         for L=0:Level
             for M=0:Level
                 if((L+M)<=Level)
                     O=O+1;
                     Ptot=0; %menghitung jumlah dari total integral batas pada bagian ke 2 dari bagian FF, 
                     for j=1:n
                         P=findQ3(x(i),y(i),nr(j),nz(j),lm(j),rb(j),zb(j),rb(j+1),zb(j+1),L,M,0,Lz,Lr);
                         Ptot=Ptot+P;
                     end
                     Qtot(i+n,O)=Q2(x(i),y(i),L,M,0,Lz,Lr,1)-Ptot;
                 end
             end
         end 
		 fprintf('loop : %d\n',i);
end
fprintf('Finish to make matrix element G and Q \n');