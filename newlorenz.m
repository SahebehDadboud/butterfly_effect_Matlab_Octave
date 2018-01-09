%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sahebeh Dadboud : 1569395
% Assignment 2 - exe 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all;
clear all;
clc;


sigma=10.0;
 b=8.0/3.0;
 r=28.0;    

%%%%%%%%%%%
% Set initial condition, step length 
% and number of steps
%%%%%%%%%%%


y1(1)=1.000;  y2(1)=0.0; y3(1)=0.0; 

h=0.01;  

xM=100.0; %number of steps

h2=h/2;

%%%%%%%%%
% Advance solution
%%%%%%%%%

t(1) = 0;
x = t(1);

i=1;
while x<xM

i=i+1;

%%%
% Step from t(i-1) to t(i) using 
% the provisional points yt1,yt2,yt3 inbetween
%%%

% 1. Compute f(t,y)

     y1m = y1(i-1);  %x
     y2m = y2(i-1);  %y
     y3m = y3(i-1);  %z
     
     
	k1 = sigma*(y2m-y1m);
	k2 = r*y1m - y2m - y1m*y3m;
  k3 = -b*y3m + y1m*y2m;
    
% 2. Set provisional values yt1

	yt11 = y1(i-1) + h2*k1;
	yt12 = y2(i-1) + h2*k2;
  yt13 = y3(i-1) + h2*k3;
    
% 3. Compute phase speed k1 at provisional points

	k11 = sigma*(yt12-yt11);
	k12 = r*yt11 - yt12 - yt11*yt13;
  k13 = -b*yt13 + yt11*yt12;
    
% 4. Set provisional values yt2

	yt21 = y1(i-1) + h2*k11;
	yt22 = y2(i-1) + h2*k12;
  yt23 = y3(i-1) + h2*k13;

% 5. Compute phase speed ft2 at provisional points

	k21 = sigma*(yt22-yt21);
	k22 = r*yt21 - yt22 - yt21*yt23;
  k23 = -b*yt23 + yt21*yt22;
    
% 6. Set provisional values yt3

	yt31 = y1(i-1) + h*k21;
	yt32 = y2(i-1) + h*k22;
  yt33 = y3(i-1) + h*k23;
    
% 7. Compute phase speed ft3 at provisional points

	k31 = sigma*(yt32-yt31);
	k32 = r*yt31 - yt32 - yt31*yt33;
  k33 = -b*yt33 + yt31*yt32;
    
% 8. Compute final phase speeds ff

	kk1 = (1/6)*( k1 + 2*k11 + 2*k21 + k31 );
	kk2 = (1/6)*( k2 + 2*k12 + 2*k22 + k32 );
  kk3 = (1/6)*( k3 + 2*k13 + 2*k23 + k33 );
    
% 9. Set new solution values

	y1(i) = y1(i-1) + h*kk1;
	y2(i) = y2(i-1) + h*kk2;
  y3(i) = y3(i-1) + h*kk3;

% Update independent variable

	t(i) = t(i-1) + h;
        x = t(i-1);	

%x

end


%%%%%%
% Plot results
%%%%%%

plot(t,y1,'r--');
plot(t,y1,'k');
hold on;
plot(t,y2,'r:');
plot(t,y3,'b:');
xlabel('t','fontsize',20);
ylabel('x(t)','fontsize',20);


figure

plot3(y1,y2,y3);
%plot(y1,y3);
xlabel('x'); ylabel('y'); zlabel('z')

 