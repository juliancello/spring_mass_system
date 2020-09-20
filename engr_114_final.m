%{ 
Author: Julian Kosanovic
Date: 8/26/18
Title: Final Presentation
Version: 1
Notes: Works currently on matlab. See line 44 for octave compatibility 
%}
clear all, clc
%Needed functions: exp(x),

%This small program can calculate the solution system of position functions
%for a simple system of two spring coupled masses. Later version will take
%bigger than 2x2 matrices for the mass and spring constants

%Variable definitions
m1=10; %later user defined
m2=20;
m=[m1,0;0,m2]; %Mass matrix, m1, m2
k1=100;
k2=150;
k3=125;
k=[-(k1+k2),k2;k2,-(k2+k3)]; %obtained by drawing diagram of system
iden=[1,0;0,1];
%n=[m,iden]
%rref(n)
A=inv(m)*k;

%Find eigenvalues. Write separate function for this eventually
a=A(1,1); %The following 4 will be input args 
b=A(1,2);
c=A(2,1);
f=A(2,2);
%Function will return two values: 
e_val_1_A=(((a+f)+sqrt((-a-f)^2-(4*(a*f-b*c))))*.5); %Eigenvalue 1
e_val_2_A=(((a+f)-sqrt((-a-f)^2-(4*(a*f-b*c))))*.5); %Eigenvalue 2
om1=sqrt(abs(e_val_1_A)); %omega
om2=sqrt(abs(e_val_2_A)); %is this needed?
%v= %Find eigenvector of e-val 1
t=0:1:10; %example. Possible user input 
B1=rref(A-(e_val_1_A*iden));
B2=rref(A-(e_val_2_A*iden)); 
if B1(1,2)==0 || B2(1,2)==0
  error('Something is wrong');
end %change to endif if running on octave
%ev1=[B1(1,2)*B1(2,2);B1(2,2)]; %eigenvectors
%ev2=[B2(1,2)*B2(2,2);B2(2,2)];
ev1=[B1(1,2);1];
ev2=[B2(1,2);1];
c1=1;
c2=c1;
c3=c1;
c4=c1; %I do not know about multivariable calc and only a bit about partial derivatives
x=ev1*(c1*cos(om1*t)+c2*sin(om1*t))+ ev2*(c3*cos(om2*t)+c4*sin(om2*t));  %(14) c2, c4 have i in them


