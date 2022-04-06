clear all vars
clc
close all
%% This script is meant to simulate acoustics wave from a baffled piston incident on a phase plate.

%% Part (a) define a geometry and mesh the geometry to perform Rayleigh intergral...
% And also to define the plane from the piston to evaluate the rayleigh
% intergral.
f=2e5;%Frequency of Ultrasonic wave
a=5e-2;%Radius of the circular piston
co=343;%sound speed in air
uo=1;% Particle velocity on the surface of the piston in m/s assume uniform across entire surface of piston
M_size=1001;% Meaning a linear line along the piston's diameter is M_size
P_piston=zeros(M_size,M_size);%Initialize the coordinate mesh on the surface of the piston->Then fill in the position of each mesh grid
Count=0;
for i=1:length(P_piston(:,1))%Index Y
    for j=1:length(P_piston(1,:))%Index X
        if ((i-(length(P_piston(:,1))+1)/2)^2+(j-(length(P_piston(:,1))+1)/2)^2)<=(length(P_piston(:,1))/2)^2
            P_piston(i,j)=1;
            Count=Count+1;
        end
    end
end
imshow(P_piston)% show an imag of the piston
S_piston=P_piston*pi*a^2/Count; %Assigna each mesh square with an area


        
