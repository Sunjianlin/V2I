function [x,y,theta,r,x_comms,y_comms,phi,rcomms] = scatterer(K,xc,yc,L,W,Dx,Dy,vx,vy);

%% 通信Rx 信息
rho = atand(vy/vx);
x_comms = xc + Dx*cosd(rho)-Dy*sind(rho); % comm Rx 初始坐标x
y_comms = yc + Dx*sind(rho)+Dy*cosd(rho); % comm Rx 初始坐标y
% x_comms = xc + sqrt(Dx^2+Dy^2)*cosd( atand(Dy/Dx)+atand(vy/vx) ); % comm Rx 初始坐标x
% y_comms = yc + sqrt(Dx^2+Dy^2)*sind( atand(Dy/Dx)+atand(vy/vx) ); % comm Rx 初始坐标y
phi = arctand(x_comms,y_comms); % comms RX的初始角度
% phi = atand(y_comms/x_comms); % comms RX的初始角度
rcomms = sqrt(x_comms^2+y_comms^2);

%% 散射点信息
% a = xc - L/2;
% b = xc + L/2;
% c = yc - W/2;
% d = yc + W/2;
%  
% % r = a + (b-a).*rand(N,1);
% 
% x = rand(K,1)*(b-a)+a;  %横坐标在(-2500，2500)范围内，生成1*10的随机矩阵
% y = rand(K,1)*(d-c)+c;  %纵坐标在(-2500，2500)范围内，生成1*10的随机矩阵
% theta = arctand(x,y);

%% 四个角散射点 信息
x0_theta1 = xc + (-L/2)*cosd(rho)-(W/2)*sind(rho); % 左上点 初始坐标
y0_theta1 = yc + (-L/2)*sind(rho)+(W/2)*cosd(rho);
% x0_theta1 = xc - sqrt(L^2/4+W^2/4)*cosd(atand(W/L)-atand(vy/vx)); % 左上点 初始坐标
% y0_theta1 = yc + sqrt(L^2/4+W^2/4)*sind(atand(W/L)-atand(vy/vx));
theta0_1 = arctand(x0_theta1,y0_theta1);

x0_theta2 = xc + (-L/2)*cosd(rho)-(-W/2)*sind(rho); % 左下点 初始坐标
y0_theta2 = yc + (-L/2)*sind(rho)+(-W/2)*cosd(rho);
% x0_theta2 = xc - sqrt(L^2/4+W^2/4)*cosd(atand(W/L)+atand(vy/vx)); % 左下点 初始坐标
% y0_theta2 = yc - sqrt(L^2/4+W^2/4)*sind(atand(W/L)+atand(vy/vx));
theta0_2 = arctand(x0_theta2,y0_theta2);

x0_theta3 = xc + (L/2)*cosd(rho)-(-W/2)*sind(rho); % 右下点 初始坐标
y0_theta3 = yc + (L/2)*sind(rho)+(-W/2)*cosd(rho);
% x0_theta3 = xc + sqrt(L^2/4+W^2/4)*cosd(atand(W/L)-atand(vy/vx)); % 右下点 初始坐标
% y0_theta3 = yc - sqrt(L^2/4+W^2/4)*sind(atand(W/L)-atand(vy/vx));
theta0_3 = arctand(x0_theta3,y0_theta3);

x0_theta4 = xc + (L/2)*cosd(rho)-(W/2)*sind(rho); % 右上点 初始坐标
y0_theta4 = yc + (L/2)*sind(rho)+(W/2)*cosd(rho);
% x0_theta4 = xc + sqrt(L^2/4+W^2/4)*cosd(atand(W/L)+atand(vy/vx)); % 右上点 初始坐标
% y0_theta4 = yc + sqrt(L^2/4+W^2/4)*sind(atand(W/L)+atand(vy/vx));
theta0_4 = arctand(x0_theta4,y0_theta4);


x0_theta5 = (xc + x0_theta1)/2; % 左上点 初始坐标
y0_theta5 = (yc + y0_theta1)/2;
% x0_theta1 = xc - sqrt(L^2/4+W^2/4)*cosd(atand(W/L)-atand(vy/vx)); % 左上点 初始坐标
% y0_theta1 = yc + sqrt(L^2/4+W^2/4)*sind(atand(W/L)-atand(vy/vx));
theta0_5 = arctand(x0_theta5,y0_theta5);

x0_theta6 = (xc + x0_theta2)/2; % 左上点 初始坐标
y0_theta6 = (yc + y0_theta2)/2;
% x0_theta1 = xc - sqrt(L^2/4+W^2/4)*cosd(atand(W/L)-atand(vy/vx)); % 左上点 初始坐标
% y0_theta1 = yc + sqrt(L^2/4+W^2/4)*sind(atand(W/L)-atand(vy/vx));
theta0_6 = arctand(x0_theta6,y0_theta6);

x0_theta7 = (xc + x0_theta3)/2; % 左上点 初始坐标
y0_theta7 = (yc + y0_theta3)/2;
% x0_theta1 = xc - sqrt(L^2/4+W^2/4)*cosd(atand(W/L)-atand(vy/vx)); % 左上点 初始坐标
% y0_theta1 = yc + sqrt(L^2/4+W^2/4)*sind(atand(W/L)-atand(vy/vx));
theta0_7 = arctand(x0_theta7,y0_theta7);

x0_theta8 = (xc + x0_theta4)/2; % 左上点 初始坐标
y0_theta8 = (yc + y0_theta4)/2;
% x0_theta1 = xc - sqrt(L^2/4+W^2/4)*cosd(atand(W/L)-atand(vy/vx)); % 左上点 初始坐标
% y0_theta1 = yc + sqrt(L^2/4+W^2/4)*sind(atand(W/L)-atand(vy/vx));
theta0_8 = arctand(x0_theta8,y0_theta8);

x = [x0_theta1;x0_theta2;x0_theta3;x0_theta4;x0_theta5;x0_theta6;x0_theta7;x0_theta8]; 
y = [y0_theta1;y0_theta2;y0_theta3;y0_theta4;y0_theta5;y0_theta6;y0_theta7;y0_theta8];

theta = [theta0_1;theta0_2;theta0_3;theta0_4;theta0_5;theta0_6;theta0_7;theta0_8];

r = sqrt(x.^2+y.^2); % 各个散射点的径向距离

