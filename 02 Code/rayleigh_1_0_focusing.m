clear all vars
clc
close all force

%% Pressure at focal plane (0); source: spherical (focused) transducer (T), phase plate: l=3 (1)
% see lines 61, 65 for curvature of transducer; see supplemental notes for more info

W_bar = waitbar(0,'Please wait...');

%% Part (a) define a geometry and mesh the geometry to perform Rayleigh intergral...
% Define the plane from the transducer to evaluate the rayleigh intergral.

F = 100e-3;         %Focal length
f = 1.092e6;        %Frequency of Ultrasonic wave
rho0 = 1000;        %Density of water
c0 = 1481;          %sound speed in water
K = 2*pi*f/c0;      %wavenumber
a = 50e-3;          %Radius of the circular transducer
u0 = 1;             %Particle velocity on the surface of the transducer in m/s assume uniform across entire surface of tranducer
%size_T=101;        %Line along the transducer's diameter
size_T=151;          %Line along the transducer's diameter
P_T=zeros(size_T,size_T); %Initialize the coordinate mesh on the surface of the transducer->Then fill in the position of each mesh grid
Count=0;
for i=1:length(P_T(:,1))%Index Y
    tic %start the timer
    for j=1:length(P_T(1,:))%Index X
        if ((i-(length(P_T(:,1))+1)/2)^2+(j-(length(P_T(:,1))+1)/2)^2)<=(length(P_T(:,1))/2)^2% Draw a circle on the square matrix
            P_T(i,j)=1;
            Count=Count+1;
        end
    end
    time=toc %stop the timer
        waitbar_show(i,P_T(:,1),time,W_bar)  
end
%imshow(P_T)     %show an image of the transducer
S_T = P_T*pi*a^2/Count; %Assign each mesh square with a differential area dS_T

%introduce phase plate
D_1 = 75e-3;        %Radius of the phase plate 150/2
z1 = 23.4e-3;        %Distance between the transducer and the phase plate
%size_1 = 101;      %Number of partition points along the entire line
size_1 = 151;        %use temporarily for quick computation
P_1 = zeros(size_1,size_1); %Initialize the pressure pattern on phase plate
Count = 0;          %Initialize counter

for i = 1:length(P_1(:,1))%Index Y
    tic
    for j = 1:length(P_1(1,:))%Index X
        Sum = 0;
        if ((i-(length(P_1(:,1))+1)/2)^2+(j-(length(P_1(1,:))+1)/2)^2)<=(length(P_1(:,1))/2)^2%Draw a circle on the square matrix
            for k = 1:length(P_T(:,1))
                for l = 1:length(P_T(1,:))
                    if P_T(k,l) ~= 0
                        Y_T = (length(P_T(:,1))/2)-k;%Positon in number of squares not converted to actual distance
                        X_T = -1*((length(P_T(1,:))/2)-l);%Positon in number of squares not converted to actual distance
                        
                        Y_1 = (length(P_1(:,1))/2)-i;%Positon in number of squares not converted to actual distance
                        X_1 = -1*((length(P_1(1,:))/2)-j);%Positon in number of squares not converted to actual distance
                        
                        Y_T = Y_T*a/(length(P_T(:,1))/2);%Covert to Real Y coordinate of the point
                        X_T = X_T*a/(length(P_T(:,1))/2);%Covert to Real X coordinate of the point
                        Z_T = 0.5*(2*F - sqrt(4*F^2 - 4*(X_T^2+Y_T^2))); %see supplemental notes for derivation of this geometry

                        Y_1 = Y_1*D_1/(length(P_1(:,1))/2);%Covert to Real Y coordinate of the point
                        X_1 = X_1*D_1/(length(P_1(:,1))/2);%Covert to Real X coordinate of the point
                        R_1 = sqrt((X_1-X_T)^2+(Y_1-Y_T)^2+(z1-Z_T)^2);%Distance between the phase plate and the circular transducer
                        
                        Sum = Sum-sqrt(-1)*rho0*K*c0*u0*S_T(k,l)*exp(sqrt(-1)*K*R_1)/(2*pi*R_1);%Perform Rayleigh integral on each mesh point on the transducer
                    end
                end
            end
            P_1(i,j)=Sum;%Assign rayleigh intergral result onto each point on the phase plate
            Count=Count+1;%Record how many partitions there are for the phase plate
       end    
    end
        time=toc
        waitbar_show(i,P_1(:,1),time,W_bar)
end

G = K*a^2/(2*F);

%these seem reasonable
%plots of pressure incident on phase plate
figure(1)
h = surf(abs(P_1)/(rho0*c0*u0*G));
h.EdgeAlpha=0; %hides mesh
title('Pressure magnitude incident on phase plate','Interpreter','latex');
axis square
set(gca,'xtick',[],'ytick',[])
xlabel('$x_0$','interpreter','latex');
ylabel('$y_0$','interpreter','latex');
view(2)

figure(2)
g = surf(abs(angle(P_1)));
g.EdgeAlpha=0; %hides mesh
axis square
%axis off
title('Pressure phase incident on phase plate','Interpreter','latex')
set(gca,'xtick',[],'ytick',[])
xlabel('$x_0$','interpreter','latex');
ylabel('$y_0$','interpreter','latex');
view(2)


%% Part(b) Now that the rayleigh intergral has been performed on the phase plate  
%Now define the necessary phase on each of the point on the
%phase plate to yield a acoustics vortex beam

S_1 = zeros(size_1,size_1); %Assign each mesh square with a differential area dS_1
Count = 0;
for i=1:length(S_1(:,1))%Index Y
    for j=1:length(S_1(1,:))%Index X
           S_1(i,j)=1; %we are integrating over a square, so every differential surface element is assigned the same value
           Count=Count+1;
    end
end
S_1 = S_1*(2*D_1)^2/Count; %Assign each mesh square with a differential area dS_1 (area of S_1 is (2D_1)^2)

%introduce focal plane
D_0 = 5e-3;         %Radius of the focal plane--change to 5e-3 once accounting for focusing 
z0 = 100e-3;        %Distance between the transducer and the focal plane
%size_0 = 101;      %Number of partition points along the entire line
size_0 = 151;        %use temporarily for quick computation
P_0 = zeros(size_0,size_0); %Initialize the pressure pattern on focal plane

for i = 1:length(P_0(:,1))%Index Y on focal plane
    tic
for j = 1:length(P_0(1,:))%Index X on focal plane
        Sum = 0;
            for k = 1:length(P_1(:,1)) %Index Y on phase plate
            for l = 1:length(P_1(1,:)) %Index X on phase plate 
                       
                Y_1 = (length(P_1(:,1))/2)-k;%Positon in number of squares not converted to actual distance
                X_1 = -1*((length(P_1(1,:))/2)-l);%Positon in number of squares not converted to actual distance
                        
                Y_0 = (length(P_0(:,1))/2)-i;%Positon in number of squares not converted to actual distance
                X_0 = -1*((length(P_0(1,:))/2)-j);%Positon in number of squares not converted to actual distance
                        
                Y_1 = Y_1*D_1/(length(P_1(:,1))/2);%Covert to Real Y coordinate of the point
                X_1 = X_1*D_1/(length(P_1(:,1))/2);%Covert to Real X coordinate of the point
                       
                Y_0 = Y_0*D_0/(length(P_0(:,1))/2);%Covert to Real Y coordinate of the point
                X_0 = X_0*D_0/(length(P_0(:,1))/2);%Covert to Real X coordinate of the point
                R0 = sqrt((X_0-X_1)^2+(Y_0-Y_1)^2+(z0-z1)^2);%Distance between the phase plate and the focal plane
                
                if X_1>=0
                   phi = atan(Y_1/X_1);
                else
                   phi=atan(Y_1/X_1)+pi;
                end
                Sum = Sum + 1/(2*pi)* S_1(k,l)* P_1(k,l)*exp(3*sqrt(-1)*phi)*(z0-z1)/R0* (-sqrt(-1)*K/R0+1/R0^2)* exp(sqrt(-1)*K*R0);%Perform 2nd Rayleigh integral on each mesh point on the phase plate                end
            end
            end
      P_0(i,j)=Sum;%Assign 2nd Rayleigh intergral result onto each point on the focal plan 
      Count=Count+1;%Record how many partitions there are for the phase plate
end
      time=toc
      waitbar_show(i,P_0(:,1),time,W_bar)
end
       
subplot(1,2,1)
h = surf(abs(P_0)/(rho0*c0*u0*G));
h.EdgeAlpha=0; %hides mesh
title('magnitude','Interpreter','latex');
axis square
set(gca,'xtick',[],'ytick',[])
xlabel('$x_0$','interpreter','latex');
ylabel('$y_0$','interpreter','latex');
view(2)

subplot(1,2,2)
g = surf(angle(P_0));
g.EdgeAlpha=0; %hides mesh
axis square
%axis off
title('phase','Interpreter','latex')
set(gca,'xtick',[],'ytick',[])
xlabel('$x_0$','interpreter','latex');
ylabel('$y_0$','interpreter','latex');
view(2)

