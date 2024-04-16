clear;clc

tic

%% 参数区
c=3e8;

dt = 0.05; % 跟踪时间间隔
fc = 30e9; % 载频
lambda = c/fc;
p = 1; % 发射功率
G = 10; % 匹配滤波增益

B = 500e6;

L = 5; % 汽车长度 
W = 2; % 汽车宽度
Dx = 1.5; % comms Rx相对质心在x-axes的距离
Dy = 0.5; % comms Rx相对质心在y-axes的距离

K = 8;

%% 目标轨迹 初始状态
vx = -20;

T = 160/abs(vx); % 跟踪时长2.5s

vy = 0;
R = sqrt(60^2+20^2);  % 质心 初始距离
theta0 = asind(20/R);  % 质心 初始角度
xc0 = R*cosd(theta0); % 质心 初始坐标x
yc0 = R*sind(theta0); % 质心 初始坐标y

[x(:,1),y(:,1),theta(:,1),r(:,1),x_comms(1),y_comms(1),phi(1),rcomms(1)] ...
    = scatterer(K,xc0,yc0,L,W,Dx,Dy,vx,vy);  % 四个角散射点坐标和角度

%% 散射点目标状态
xc(1)=xc0;
yc(1)=yc0;
for i = 2:T/dt+1
  
    xc(i) = xc(i-1) + vx*dt;  % 质心 true 轨迹
    yc(i) = yc(i-1) + vy*dt;

    [x(:,i),y(:,i),theta(:,i),r(:,i),x_comms(i),y_comms(i),phi(i),rcomms(i)]...
        = scatterer(K,xc(i),yc(i),L,W,Dx,Dy,vx,vy);  % 四个角散射点 true 轨迹
    
end

%% 预测算法
phi_pre(1:3) = phi(1:3); % comms Rx初始角度预测
rcomms_pre(1:3) = rcomms(1:3);
xc_pre(1:3) = xc(1:3); yc_pre(1:3) = yc(1:3);
xc_hat(1:2) = xc(1:2); yc_hat(1:2) = yc(1:2);

Nr = 128; % RSU接收机天线数
sigma2 = 0.15; % echo噪声方差
% a2 = 3.5e-2; % 可由a0控制距离SNR
% a1 = 1.05e-2; % 可由a1控制角度SNR
% a3 = 1e-2;
a2 = 2.0e-2; % 可由a0控制距离SNR
a1 = 5.05e-3; % 可由a1控制角度SNR
a3 = 1e-2;
KK = 250;

vx = abs(vx);

for kk =1:KK
    
    for i = 1:T/dt+1
        
        Nt(i) = antennanumber(5.5,rcomms(i),phi(i));  % 天线数 计算
        Nt(i) = min(Nt(i),128);
      
%         if Nt(i)>128
%             Nt(i) = 128;
%         end
        
        a0 = 1/sqrt(Nt(i))*exp(-j*pi*[0:Nt(i)-1]'*cosd(phi(i)));
        Kr = sqrt(Nt(i)*Nr);  % 阵列增益factor
        
        beta(i) = 1/(2*rcomms(i))^2;
        
        for ii =1:K
            a(:,ii) = 1/sqrt(Nt(i))*exp(-j*pi*[0:Nt(i)-1]'*cosd(theta(ii,i)));
            SNR_d(ii) =  p*G*abs(beta(i))^2*abs(a(:,ii)'*a0)^2/sigma2;
        end
        
        % 各散射点估计的方差
        var_d(:,i) = a2*3/( 2*pi^2*Nt(i)*Nr)./SNR_d'; % 4x1列向量 
        var_a(:,i) = a1*1/( 2*pi^2*Nt(i)*Nr)./SNR_d'; % 4x1列向量
        var_u(:,i) = a3*1/( 2*pi^2*Nt(i)*Nr)./SNR_d'; % 4x1列向量
        
        %% %%%%%%%%%%%%%%%%%% 各散射点的高分辨 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        r_est(:,i) = r(:,i) + sqrt(var_d(:,i)).*randn(K,1);      % 各散射点距离估计
        angle(:,i) = theta(:,i) + sqrt(var_a(:,i)).*randn(K,1);  % 各散射点角度估计
        
%         vx = 20;
        u(:,i) = 2*vx*cosd(theta(:,i))*fc./c;
        u_est(:,i) = u(:,i) + sqrt(var_u(:,i)).*randn(K,1);  % 各散射点角度估计
        
        clear a0 b a yy X   % 必须clear 因为每个时刻的Nt不一样，导致数据尺寸不同
        
        %% %%%%%%%%%%%%%%%%%% 测量方程：通信Rx %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        xc_hat(i) = sum(r_est(:,i).*cosd(angle(:,i)))/K; % 质心坐标
        yc_hat(i) = sum(r_est(:,i).*sind(angle(:,i)))/K;
        x_comms_hat(i) = xc_hat(i) + sqrt(Dx^2+Dy^2)*cosd( atand(Dy/Dx)+atand(vy/vx) ); % comm Rx 初始坐标x
        y_comms_hat(i) = yc_hat(i) + sqrt(Dx^2+Dy^2)*sind( atand(Dy/Dx)+atand(vy/vx) ); % comm Rx 初始坐标y
        
        phi_hat(kk,i) = arctand(x_comms_hat(i),y_comms_hat(i));  % 估计的通信Rx角度
        tt(i) =  arctand(x_comms_hat(i),y_comms_hat(i)) - arctand(x_comms(i),y_comms(i));
        sigma_n(kk,i) = abs( tt(i) );  % 通信Rx 角度标准差
        
        range_hat(kk,i) = sqrt(x_comms_hat(i)^2+y_comms_hat(i)^2);
        sigma_d(kk,i) = abs( sqrt(x_comms_hat(i)^2+y_comms_hat(i)^2) - sqrt(x_comms(i)^2+y_comms(i)^2)); % 通信Rx 距离标准差

        v_hat(kk,i) = c/(2*fc)*sum(u_est(:,i).*cosd(angle(:,i))./var_u(:,i))./sum((cosd(angle(:,i))).^2./var_u(:,i)); % MLE 速度的最优估计
        sigma_v(kk,i) = abs(v_hat(kk,i)-vx);  % 车辆速度标准差
        
        y_measurement(kk,:,i) = [phi_hat(kk,i),range_hat(kk,i),v_hat(kk,i)].';  % 测量方程

        %% %%%%%%%%%%%%%%%%%% Taylor公式近似 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         K=4; % 4个散射点
        Delta_x = sum(r_est(:,i).*cosd(angle(:,i))) + K*Dx;
        Delta_y = sum(r_est(:,i).*sind(angle(:,i))) + K*Dy;
        % 角度
        dg_r(:,i) = ( sind(angle(:,i)).*Delta_x - cosd(angle(:,i)).*Delta_y )...
            ./( Delta_x.^2 + Delta_y.^2 );
        dg_a(:,i) = ( r_est(:,i).*cosd(angle(:,i)).*Delta_x + r_est(:,i).*sind(angle(:,i)).*Delta_y )...
            ./( Delta_x.^2 + Delta_y.^2 );
        J = [dg_r(:,i);dg_a(:,i)].'; % Jacobian矩阵
        Sigma = diag([var_d(:,i);var_a(:,i)]);
        sigma_n_taylor(kk,i) = sqrt(J*Sigma*J.'); % 通信Rx 角度标准差 Taylor近似
        
        % 距离
        df_r(:,i) = ( cosd(angle(:,i)).*Delta_x + sind(angle(:,i)).*Delta_y )...
            ./(K* sqrt(Delta_x.^2 + Delta_y.^2) );
        df_a(:,i) = ( -r_est(:,i).*sind(angle(:,i)).*Delta_x + r_est(:,i).*cosd(angle(:,i)).*Delta_y )...
            ./(K* sqrt(Delta_x.^2 + Delta_y.^2) );
        F = [df_r(:,i);df_a(:,i)].'; % Jacobian矩阵
        sigma_n_taylor_d(kk,i) = sqrt(F*Sigma*F.'); % 通信Rx 角度标准差 Taylor近似

        A = 2*fc/c.*cosd(angle(:,i));
        Qv = diag(var_u(:,i));
        sigma_n_taylor_v(kk,i) = sqrt(inv(A'*inv(Qv)*A));
    end
    
end

% sigma2_n = sum(sigma_n.^2)/KK;
% sigma2_d = sum(sigma_d.^2)/KK;
% sigma2_v = sum(sigma_v.^2)/KK;
% 
%  plot(sigma2_n,'r');hold on

sigma2_n = sum(sigma_n_taylor.^2)/KK;
sigma2_d = sum(sigma_n_taylor_d.^2)/KK;
sigma2_v = sum(sigma_n_taylor_v.^2)/KK;

% sigma2_n_taylor = sum(sigma_n_taylor.^2)/KK;
% sigma2_taylor_d = sum(sigma_n_taylor_d.^2)/KK;
% sigma2_taylor_v = sum(sigma_n_taylor_v.^2)/KK;

% plot(sigma2_n,'b')

clear var_d var_a



%% 方案2
% KK = 100;
Qw = diag([0.01^2,0.1^2,0.25^2]);
M0 = 128; % 窄波束的发射天线数
Kc_narrow = sqrt(M0*1);  % 窄波束固定32根天线
rho = linspace(0,1,10000); % grid search区间

for kk =1:KK  %% 计算多次，求achievable rate的平均结果
    
%     Qz = diag([sigma2_n(i),sigma2_d(i),sigma2_v(i)]);
    Mn = zeros(3,3);
    phi_hat(kk,1) = phi(1);
    range_hat(kk,1) = rcomms(1);
    v_hat(kk,1) = abs(vx);

    %% 预测算法 方案2
        
    rho0(1) = 0; % 初始点 取rho = 1
    Nr = 128;  % 接收天线数 固定 最大值
    for i = 1:T/dt
        
        Qz = diag([sigma2_n(i),sigma2_d(i),sigma2_v(i)]);

        % 跟踪算法 EKF
        %%  %%%%%%%%%%%%%%%%%%%%%% 1. 状态预测  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        z_hat(:,i) = [phi_hat(kk,i),range_hat(kk,i),v_hat(kk,i)].';
        
        phi_pre(i+1) = z_hat(1,i) + 1/z_hat(2,i)*z_hat(3,i)*dt*sind(z_hat(1,i)) + sqrt(Qw(1,1))*randn(1);
        rcomms_pre(i+1) = z_hat(2,i) - z_hat(3,i)*dt*cosd(z_hat(1,i)) + sqrt(Qw(2,2))*randn(1);
        v_pre(i+1) = z_hat(3,i) + sqrt(Qw(3,3))*randn(1);
        z_pre(:,i+1) = [phi_pre(i+1),rcomms_pre(i+1),v_pre(i+1)].';  

        % 测量方程
        Nt(i+1) = antennanumber(5.5,rcomms_pre(i+1),phi_pre(i+1));  % 天线数 计算
        Nt(i) = min(Nt(i),128);
      
%         if Nt(i)>128
%             Nt(i) = 128;
%         end
        
       %% 优化问题求解 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         p = 1;
        aref = 1;
        an(i+1) = aref./rcomms_pre(i).*exp(j*2*pi*fc/c.*rcomms_pre(i+1)); % LoS路径损失，求优化时未知，用预测值代替
        Kc = sqrt(Nt(i+1)*1); % 通信Rx单天线

        sigma2C = 1;
        h_wide(i+1) = norm(an(i+1)*Kc).^2/sigma2C;
        R_wide(i+1) = log2(1+p*h_wide(i+1));
        
        %     a_narrow = 1/sqrt(M0)*exp(-j*pi*[0:M0-1]'*cosd(phi(i)));
        %     a0_narrow = 1/sqrt(M0)*exp(-j*pi*[0:M0-1]'*cosd(phi_est(i)));  % 窄波束用估计角
        h_narrow(i+1) = norm(an(i+1)*Kc_narrow).^2/sigma2C;
        %     h_narrow(i) = norm(an(i)*Kc_narrow*a_narrow'*a0_narrow).^2/sigma2C;
        R_narrow(i+1) = log2(1+p*h_narrow(i+1));
        
%         delta(i) = 101/M0/sind(phi_pre(i)); % 半波束宽度 窄波束 固定M0根天线
        delta(i+1) = 101/M0/sind(phi_pre(i+1)); % 半波束宽度 窄波束 固定M0根天线  受限于phi_hat(i)而不是phi_pre(i)

        delta_wide(i+1) = 101/Nt(i+1)/sind(phi_pre(i+1)); % 半波束宽度 窄波束 固定M0根天线
%         f = rho.*erf(delta_wide(i+1).*sqrt(rho)./sqrt(2)./sqrt(sigma2_n(i+1))).*R_wide(i+1) + (1-rho).*erf(delta(i+1).*sqrt(rho)./sqrt(2)./sqrt(sigma2_n(i+1))).*R_narrow(i+1);  % 注意这里不应该用质心方差，用comms Rx方差
        f = rho.*R_wide(i+1) + (1-rho).*erf(delta(i+1).*sqrt(rho)./sqrt(2)./sqrt(sigma2_n(i+1))).*R_narrow(i+1);  % 注意这里不应该用质心方差，用comms Rx方差
               
        
        [R_whole(i+1),b0(i+1)] = max(f);
                        
        rho0(i+1) = rho(b0(i+1));  %找到最大值的位置
%         rho0(i+1) = rho(b0(i+1))+0.01;  %找到最大值的位置

        PA(i) = erf(delta(i+1).*sqrt(rho0(i+1))./sqrt(2)./sqrt(sigma2_n(i+1)));

        a0 = 1/sqrt(Nt(i+1))*exp(-j*pi*[0:Nt(i+1)-1]'*cosd(phi_pre(i+1)));
        Kr = sqrt(Nt(i+1)*Nr);  % 阵列增益factor
        
        beta(i+1) = 1/(2*rcomms(i+1))^2;  % 路径损耗 与真实的rcomms有关，与预测值无关
%         beta(i) = 1/(2*rcomms_pre(i))^2;

        for ii =1:K
            a(:,ii) = 1/sqrt(Nt(i+1))*exp(-j*pi*[0:Nt(i+1)-1]'*cosd(theta(ii,i+1)));
            %         SNR_d(ii) =  p*G*Kr^2*abs(beta(i))^2*abs(a(:,ii)'*a0)^2/sigma2;
            SNR_d(ii) =  rho0(i+1)*p*G*abs(beta(i+1))^2*abs(a(:,ii)'*a0)^2/sigma2;
        end
        
        % 各散射点估计的方差
        var_d(:,i+1) = a2*3/( 2*pi^2*Nt(i+1)*Nr)./SNR_d'; % 4x1列向量
        var_a(:,i+1) = a1*1/( 2*pi^2*Nt(i+1)*Nr)./SNR_d'; % 4x1列向量
        var_u(:,i+1) = a3*1/( 2*pi^2*Nt(i+1)*Nr)./SNR_d'; % 4x1列向量

        % 各散射点的高分辨 
        r_est(:,i+1) = r(:,i+1) + sqrt(var_d(:,i+1)).*randn(K,1);      % 各散射点距离估计
        angle(:,i+1) = theta(:,i+1) + sqrt(var_a(:,i+1)).*randn(K,1);  % 各散射点角度估计
        
%         vx = 30;
        u(:,i+1) = 2*vx*cosd(theta(:,i+1))*fc./c;
        u_est(:,i+1) = u(:,i+1) + sqrt(var_u(:,i+1)).*randn(K,1);  % 各散射点角度估计
        
        clear a0 b a yy X   % 必须clear 因为每个时刻的Nt不一样，导致数据尺寸不同
        
        % 测量方程：通信Rx 
        xc_hat(i+1) = sum(r_est(:,i+1).*cosd(angle(:,i+1)))/K; % 质心坐标
        yc_hat(i+1) = sum(r_est(:,i+1).*sind(angle(:,i+1)))/K;
        x_comms_hat(i+1) = xc_hat(i+1) + sqrt(Dx^2+Dy^2)*cosd( atand(Dy/Dx)+atand(vy/vx) ); % comm Rx 初始坐标x
        y_comms_hat(i+1) = yc_hat(i+1) + sqrt(Dx^2+Dy^2)*sind( atand(Dy/Dx)+atand(vy/vx) ); % comm Rx 初始坐标y
        
        phi_hat(kk,i+1) = arctand(x_comms_hat(i+1),y_comms_hat(i+1));  % 估计的通信Rx角度        
        range_hat(kk,i+1) = sqrt(x_comms_hat(i+1)^2+y_comms_hat(i+1)^2);
        v_hat(kk,i+1) = c/(2*fc)*sum(u_est(:,i+1).*cosd(angle(:,i+1))./var_u(:,i+1))./sum((cosd(angle(:,i+1))).^2./var_u(:,i+1)); % MLE 速度的最优估计
        y_m(:,i+1) = [phi_hat(kk,i+1),range_hat(kk,i+1),v_hat(kk,i+1)].';  % 测量方程

%         y_m(:,i+1) = y_measurement(kk,:,i+1);  % 测量方程


        %%  %%%%%%%%%%%%%%%%%%%%%% 2. 线性化 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Hn_1 = [ 1+ (z_hat(3,i)*dt*cosd(z_hat(1,i)))/(z_hat(2,i)), -(z_hat(3,i)*dt*sind(z_hat(1,i)))/(z_hat(2,i))^2, 0 ;
                 z_hat(3,i)*dt*sind(z_hat(1,i)), 1, -dt*cosd(z_hat(1,i));
                 0, 0, 1];
        %%  %%%%%%%%%%%%%%%%%%%%%% 3. MSE矩阵预测 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Mnln_1 = Hn_1*Mn*Hn_1'+Qw;
        %%  %%%%%%%%%%%%%%%%%%%%%% 4. 卡尔曼增益矩阵计算 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Kn = Mnln_1*inv(Mnln_1+Qz);
        %%  %%%%%%%%%%%%%%%%%%%%%% 5. 状态跟踪 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
        z_hat(:,i+1) = z_pre(:,i+1) + Kn*(y_m(:,i+1) - z_pre(:,i+1) ); % y_m: measurement
        %%  %%%%%%%%%%%%%%%%%%%%%% 6. MSE矩阵更新 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Mn = (eye(3)-Kn)*Mnln_1;

        % 计算实际的achievable ate
        an(i+1) = aref./rcomms(i).*exp(j*2*pi*fc/c.*rcomms(i+1)); % LoS路径损失，计算真实值时，用真实的rcomms
        a_wide = 1/sqrt(Nt(i+1))*exp(-j*pi*[0:Nt(i+1)-1]'*cosd(phi(i+1)));
        a0_wide = 1/sqrt(Nt(i+1))*exp(-j*pi*[0:Nt(i+1)-1]'*cosd(phi_pre(i+1))); % 宽波束用预测角
        h_wide(i+1) = norm(an(i+1)*Kc*a_wide'*a0_wide).^2/sigma2C;
%         h_wide(i+1) = norm(an(i+1)*Kc).^2/sigma2C;
        R_wide(i+1) = log2(1+p*h_wide(i+1));
        
%         M0 = 64;
        Kc_narrow = sqrt(M0*1);  % 窄波束固定32根天线
        a_narrow = 1/sqrt(M0)*exp(-j*pi*[0:M0-1]'*cosd(phi(i+1)));
        a0_narrow = 1/sqrt(M0)*exp(-j*pi*[0:M0-1]'*cosd(phi_hat(kk,i+1)));  % 窄波束用估计角
%         h_narrow(i+1) = norm(an(i+1)*Kc_narrow).^2/sigma2C;
        h_narrow(i+1) = norm(an(i+1)*Kc_narrow*a_narrow'*a0_narrow).^2/sigma2C;
        R_narrow(i+1) = log2(1+p*h_narrow(i+1));
        
        R_sum0(i+1) = rho0(i+1)*R_wide(i+1)+(1-rho0(i+1))*R_narrow(i+1);
        
        clear a_c a0_c a_wide a0_wide a_narrow a0_narrow % !!!!! 必须清空

    end
    
    Phi_pre(kk,:) = phi_pre;
    Rcomms_pre(kk,:) = rcomms_pre;
    V_pre(kk,:) = v_pre;
%     Rn(kk,:) = Rn0;

    R_sum(kk,:) = R_sum0;
    R_obj(kk,:) = R_whole; % 优化的目标函数
    rho_whole(kk,:) = rho0;
    
    %% Predicted RMSE计算
    rmse_phi(kk,:)= (phi_pre-phi).^2;
    rmse_r(kk,:)= (rcomms_pre-rcomms).^2;
    rmse_v(kk,:)= (v_pre-abs(vx)).^2;

    rmse_phi_hat(kk,:)= (phi_hat(kk,:)-phi).^2;

end

% KK=50;
Phi_ave = sum(Phi_pre)/KK;
Rcomms_ave = sum(Rcomms_pre)/KK;
V_ave  = sum(V_pre)/KK;
Rn_ave = sum(R_sum)/KK;
Rn_obj_ave = sum(R_obj)/KK;
rho_whole_ave = sum(rho_whole)/KK;

Rmse_phi = sqrt(sum(rmse_phi)/KK);
Rmse_r = sqrt(sum(rmse_r)/KK);
Rmse_v = sqrt(sum(rmse_v)/KK);

Rmse_phi_hat = sqrt(sum(rmse_phi_hat)/KK);


%% 将方案2的预测角和预测距离记录下来
phi_scheme2 = Phi_ave;
rcomms_scheme2 = Rcomms_ave;
phi_hat_scheme2 = sum(phi_hat)/KK;
Rmse_phi_scheme2 = Rmse_phi;
    
Rmse_phi_scheme2_2 = Rmse_phi_hat;

%% 方案1
for kk =1:KK;
    
    % 0. 参数设置
    Qw = diag([0.01^2,0.1^2,0.25^2]);
%     Qz = diag([sigma2_n(i),sigma2_d(i),sigma2_v(i)]);
    Mn = zeros(3,3);
    phi_hat(kk,1) = phi(1);
    range_hat(kk,1) = rcomms(1);
    v_hat(kk,1) = abs(vx);
    for i = 1:T/dt
        
        Qz = diag([sigma2_n(i),sigma2_d(i),sigma2_v(i)]);

        % 跟踪算法 EKF
        %%  %%%%%%%%%%%%%%%%%%%%%%1. 状态预测  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        z_hat(:,i) = [phi_hat(kk,i),range_hat(kk,i),v_hat(kk,i)].';
        
        phi_pre(i+1) = z_hat(1,i) + 1/z_hat(2,i)*z_hat(3,i)*dt*sind(z_hat(1,i)) + sqrt(Qw(1,1))*randn(1);
        rcomms_pre(i+1) = z_hat(2,i) - z_hat(3,i)*dt*cosd(z_hat(1,i)) + sqrt(Qw(2,2))*randn(1);
        v_pre(i+1) = z_hat(3,i) + sqrt(Qw(3,3))*randn(1);
        z_pre(:,i+1) = [phi_pre(i+1),rcomms_pre(i+1),v_pre(i+1)].';  

        % 测量方程
        Nt(i+1) = antennanumber(5.5,rcomms_pre(i+1),phi_pre(i+1));  % 天线数 计算
        Nt(i) = min(Nt(i),128);
      
%         if Nt(i)>128
%             Nt(i) = 128;
%         end
      
        a0 = 1/sqrt(Nt(i+1))*exp(-j*pi*[0:Nt(i+1)-1]'*cosd(phi_pre(i+1)));
        Kr = sqrt(Nt(i+1)*Nr);  % 阵列增益factor
        
        beta(i+1) = 1/(2*rcomms(i+1))^2;  % 路径损耗 与真实的rcomms有关，与预测值无关
%         beta(i) = 1/(2*rcomms_pre(i))^2;

        for ii =1:K
            a(:,ii) = 1/sqrt(Nt(i+1))*exp(-j*pi*[0:Nt(i+1)-1]'*cosd(theta(ii,i+1)));
            %         SNR_d(ii) =  p*G*Kr^2*abs(beta(i))^2*abs(a(:,ii)'*a0)^2/sigma2;
            SNR_d(ii) =  p*G*abs(beta(i+1))^2*abs(a(:,ii)'*a0)^2/sigma2;
        end
        
        % 各散射点估计的方差
        var_d(:,i+1) = a2*3/( 2*pi^2*Nt(i+1)*Nr)./SNR_d'; % 4x1列向量
        var_a(:,i+1) = a1*1/( 2*pi^2*Nt(i+1)*Nr)./SNR_d'; % 4x1列向量
        var_u(:,i+1) = a3*1/( 2*pi^2*Nt(i+1)*Nr)./SNR_d'; % 4x1列向量

        % 各散射点的高分辨 
        r_est(:,i+1) = r(:,i+1) + sqrt(var_d(:,i+1)).*randn(K,1);      % 各散射点距离估计
        angle(:,i+1) = theta(:,i+1) + sqrt(var_a(:,i+1)).*randn(K,1);  % 各散射点角度估计
        
%         vx = 30;
        u(:,i+1) = 2*vx*cosd(theta(:,i+1))*fc./c;
        u_est(:,i+1) = u(:,i+1) + sqrt(var_u(:,i+1)).*randn(K,1);  % 各散射点角度估计
        
        clear a0 b a yy X   % 必须clear 因为每个时刻的Nt不一样，导致数据尺寸不同
        
        % 测量方程：通信Rx 
        xc_hat(i+1) = sum(r_est(:,i+1).*cosd(angle(:,i+1)))/K; % 质心坐标
        yc_hat(i+1) = sum(r_est(:,i+1).*sind(angle(:,i+1)))/K;
        x_comms_hat(i+1) = xc_hat(i+1) + sqrt(Dx^2+Dy^2)*cosd( atand(Dy/Dx)+atand(vy/vx) ); % comm Rx 初始坐标x
        y_comms_hat(i+1) = yc_hat(i+1) + sqrt(Dx^2+Dy^2)*sind( atand(Dy/Dx)+atand(vy/vx) ); % comm Rx 初始坐标y
        
        phi_hat(kk,i+1) = arctand(x_comms_hat(i+1),y_comms_hat(i+1));  % 估计的通信Rx角度        
        range_hat(kk,i+1) = sqrt(x_comms_hat(i+1)^2+y_comms_hat(i+1)^2);
        v_hat(kk,i+1) = c/(2*fc)*sum(u_est(:,i+1).*cosd(angle(:,i+1))./var_u(:,i+1))./sum((cosd(angle(:,i+1))).^2./var_u(:,i+1)); % MLE 速度的最优估计
        y_m(:,i+1) = [phi_hat(kk,i+1),range_hat(kk,i+1),v_hat(kk,i+1)].';  % 测量方程

%         y_m(:,i+1) = y_measurement(kk,:,i+1);  % 测量方程

        %%  %%%%%%%%%%%%%%%%%%%%%% 2. 线性化 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Hn_1 = [ 1+ (z_hat(3,i)*dt*cosd(z_hat(1,i)))/(z_hat(2,i)), -(z_hat(3,i)*dt*sind(z_hat(1,i)))/(z_hat(2,i))^2, 0 ;
                 z_hat(3,i)*dt*sind(z_hat(1,i)), 1, -dt*cosd(z_hat(1,i));
                 0, 0, 1];
        %%  %%%%%%%%%%%%%%%%%%%%%% 3. MSE矩阵预测 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Mnln_1 = Hn_1*Mn*Hn_1'+Qw;
        %%  %%%%%%%%%%%%%%%%%%%%%% 4. 卡尔曼增益矩阵计算 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Kn = Mnln_1*inv(Mnln_1+Qz);
        %%  %%%%%%%%%%%%%%%%%%%%%% 5. 状态跟踪 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
        z_hat(:,i+1) = z_pre(:,i+1) + Kn*(y_m(:,i+1) - z_pre(:,i+1) ); % y_m: measurement
        %%  %%%%%%%%%%%%%%%%%%%%%% 6. MSE矩阵更新 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Mn = (eye(3)-Kn)*Mnln_1;

        % 计算achievable rate 
        pc = 1;
        aref = 1;
        an = aref./rcomms.*exp(j*2*pi*fc/c.*rcomms); % 信号衰减系数，与真实的rcomms有关，与预测值无关
        M = 1;  % 通信Rx天线数
        sigma2C = 1;  % 通信噪声功率        
        Kc = sqrt(Nt(i+1)*M);
        
        a_c = 1/sqrt(Nt(i+1))*exp(-j*pi*[0:Nt(i+1)-1]'*cosd(phi(i+1)));
        a0_c = 1/sqrt(Nt(i+1))*exp(-1i*pi*[0:Nt(i+1)-1]'*cosd(phi_pre(i+1)));
        hn(i+1) = norm(an(i+1)*Kc*a_c'*a0_c).^2/sigma2C;
        Rn0(i+1) = log2(1+pc*hn(i+1));
        
        clear a_c a0_c % 必须清空，很重要
        
    end
    Phi_pre(kk,:) = phi_pre;
    Rcomms_pre(kk,:) = rcomms_pre;
    V_pre(kk,:) = v_pre;
    Rn(kk,:) = Rn0;
    
    %% Predicted RMSE计算
    rmse_phi(kk,:)= (phi_pre-phi).^2;
    rmse_r(kk,:)= (rcomms_pre-rcomms).^2;
    rmse_v(kk,:)= (v_pre-abs(vx)).^2;

end
Phi_ave = sum(Phi_pre)/KK;
Rcomms_ave = sum(Rcomms_pre)/KK;
V_ave  = sum(V_pre)/KK;
Rn_ave1 = sum(Rn)/KK;

Rmse_phi = sqrt(sum(rmse_phi)/KK);
Rmse_r = sqrt(sum(rmse_r)/KK);
Rmse_v = sqrt(sum(rmse_v)/KK);

Rn_ave_scheme1 = Rn_ave1; % 记录方案1的数据


%% 将方案1的预测角和预测距离记录下来
phi_scheme1 = Phi_ave;
rcomms_scheme1 = Rcomms_ave;
Rmse_phi_scheme1 = Rmse_phi;


figure(2)
subplot(131)
plot([1:T/dt+1]*dt,Phi_ave,'b','linewidth',1);
hold on;
plot([1:T/dt+1]*dt,phi,'m--','linewidth',1);
xlabel('Time (s)','Fontname', 'Times New Roman','FontSize',14);
ylabel('Angle (degree)','Fontname', 'Times New Roman','FontSize',14);
h = legend('Predicted Angle','Real Angle')
set(h,'FontName','Times New Roman','FontSize',14,'FontWeight','normal')
grid on

subplot(132)
plot([1:T/dt+1]*dt,Rcomms_ave,'b','linewidth',1);
hold on;
plot([1:T/dt+1]*dt,rcomms,'m--','linewidth',1);
xlabel('Time (s)','Fontname', 'Times New Roman','FontSize',14);
ylabel('Distance (m)','Fontname', 'Times New Roman','FontSize',14);
h = legend('Predicted Distance','Real Distance')
set(h,'FontName','Times New Roman','FontSize',14,'FontWeight','normal')
grid on

subplot(133)
plot([3:T/dt+1]*dt,V_ave(3:end),'b','linewidth',1);
hold on;
plot([1:T/dt+1]*dt,vx,'m--','linewidth',1);
xlabel('Time (s)','Fontname', 'Times New Roman','FontSize',14);
ylabel('Velocity (m/s)','Fontname', 'Times New Roman','FontSize',14);
h = legend('Predicted Velocity','Real Velocity')
set(h,'FontName','Times New Roman','FontSize',14,'FontWeight','normal')
grid on
% xlim([0,8])
% ylim([-90,180])



figure(3)
plot(r_est(1,1:1:end).*cosd(angle(1,1:1:end)),r_est(1,1:1:end).*sind(angle(1,1:1:end)),'ko');hold on;
plot(r_est(2,1:1:end).*cosd(angle(2,1:1:end)),r_est(2,1:1:end).*sind(angle(2,1:1:end)),'ko')
plot(r_est(3,1:1:end).*cosd(angle(3,1:1:end)),r_est(3,1:1:end).*sind(angle(3,1:1:end)),'ko')
plot(r_est(4,1:1:end).*cosd(angle(4,1:1:end)),r_est(4,1:1:end).*sind(angle(4,1:1:end)),'ko')
plot(r_est(5,1:1:end).*cosd(angle(5,1:1:end)),r_est(5,1:1:end).*sind(angle(4,1:1:end)),'ko')
plot(r_est(6,1:1:end).*cosd(angle(6,1:1:end)),r_est(6,1:1:end).*sind(angle(4,1:1:end)),'ko')
plot(r_est(7,1:1:end).*cosd(angle(7,1:1:end)),r_est(7,1:1:end).*sind(angle(4,1:1:end)),'ko')
f8=plot(r_est(8,1:1:end).*cosd(angle(8,1:1:end)),r_est(8,1:1:end).*sind(angle(4,1:1:end)),'ko')
xc_pre = Rcomms_ave.*cosd(Phi_ave);
yc_pre = Rcomms_ave.*sind(Phi_ave);
f9=plot(xc_pre,yc_pre,'b','linewidth',4);hold on
f10=plot(x_comms,y_comms,'y','linewidth',2);
% plot(0,0,'k^','linewidth',3);hold on;
% ylim([0,60])
grid on
xlabel('x-axis (m)','Fontname', 'Times New Roman','FontSize',14);
ylabel('y-axis (m)','Fontname', 'Times New Roman','FontSize',14);
% h = legend('Estimated Scatterer 1: Upper Left','Estimated Scatterer 2: Lower Left',...
%     'Estimated Scatterer 3: Lower Right','Estimated Scatterer 4: Upper Right','Estimated Scatterer 1: Upper Left','Estimated Scatterer 2: Lower Left',...
%     'Estimated Scatterer 3: Lower Right','Estimated Scatterer 4: Upper Right','Predicted Comms Rx','True Comms Rx')
h = legend([f8,f9,f10],'Scatterers','Predicted Comms Rx','True Comms Rx')
set(h,'FontName','Times New Roman','FontSize',14,'FontWeight','normal')


clear Phi_pre

%% ABP算法
for kk =1:KK  %% 计算多次，求achievable rate的平均结果
    M0 =128;
    Nr = 128;  % 接收天线数 固定 最大值
    phi0 = phi/180*pi;
    theta_hat(1) = phi0(1);
    for i = 2:T/dt

        D = pi/32;
        a_abp = 1/sqrt(M0)*exp(-j*pi*[0:M0-1]'*cos(phi0(i)));  
        fl_abp = 1/sqrt(M0)*exp(j*[0:M0-1]'*( -pi*cos(theta_hat(i-1))-D ) );  
        fr_abp = 1/sqrt(M0)*exp(j*[0:M0-1]'*( -pi*cos(theta_hat(i-1))+D ) );  
        
        % 通信接收信号
        aref = 1;
        sigma2C = 1;
        an(i) = aref./rcomms(i).*exp(j*2*pi*fc/c.*rcomms(i)); % LoS路径损失，计算真实值时，用真实的rcomms
        y1 = an(i)*a_abp'*fl_abp + sqrt(sigma2C)*randn(1,1);
        y2 = an(i)*a_abp'*fr_abp + sqrt(sigma2C)*randn(1,1);
        
%         lambda = (norm(y1)^2-norm(y2)^2)/ (norm(y1)^2+norm(y2)^2);
        lambda = (norm(an(i)*a_abp'*fl_abp)^2-norm(an(i)*a_abp'*fr_abp)^2)/ (norm(an(i)*a_abp'*fl_abp)^2+norm(an(i)*a_abp'*fr_abp)^2);
        
        % 算法
        theta_hat(i) = acos( -( -pi*cos(theta_hat(i-1)) - asin( (lambda*sin(D)-lambda*sqrt(1-lambda^2)*sin(D)*cos(D) )...
                                ./(sin(D)^2+lambda^2*cos(D)^2) ) )/pi );
        
        
        % 计算实际的achievable ate
        Kc_narrow = sqrt(M0*1);  % 窄波束固定32根天线
%         a_narrow = 1/sqrt(M0)*exp(-j*pi*[0:M0-1]'*cosd(phi(i+1)));
        a0_narrow = 1/sqrt(M0)*exp(-j*pi*[0:M0-1]'*cos(theta_hat(i)));  % 窄波束用估计角
%         h_narrow(i+1) = norm(an(i+1)*Kc_narrow).^2/sigma2C;
        h_narrow(i) = norm(an(i)*Kc_narrow*a_abp'*a0_narrow).^2/sigma2C;
        R0_abp(kk,i) = log2(1+p*h_narrow(i));

        clear a_c a0_c a_wide a0_wide a_narrow a0_narrow % !!!!! 必须清空

    end
    
    Phi_pre(kk,:) = theta_hat;

    rmse_phi_abp(kk,:) = (theta_hat/pi*180-phi(2:end)).^2;

end

Rmse_phi_abp = sqrt(sum(rmse_phi_abp)/KK);

R_abp = sum(R0_abp)/KK;

% plot([2:T/dt]*dt,R_abp(2:end),'k-.','linewidth',2)

Phi_ave = sum(Phi_pre)/KK;

%% 将ABP的“估计角”记录下来
phi_abp = Phi_ave;
Rmse_phi_abp = Rmse_phi_abp;


clear Phi_pre

%% Liu .TWC 点目标EKF
for kk =1:KK;
    
    % 0. 参数设置
    Qw = diag([0.01^2,0.1^2,0.25^2]);
%     Qz = diag([sigma2_n(i),sigma2_d(i),sigma2_v(i)]);
    Mn = zeros(3,3);
    phi_hat(kk,1) = theta(1,1);
    range_hat(kk,1) = rcomms(1);
    v_hat(kk,1) = abs(vx);
    for i = 1:T/dt
                
        clear a var_d var_a var_u r_est angle u u_est

        Qz = diag([sigma2_n(i),sigma2_d(i),sigma2_v(i)]);

        % 跟踪算法 EKF
        %%  %%%%%%%%%%%%%%%%%%%%%%1. 状态预测  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        z_hat(:,i) = [phi_hat(kk,i),range_hat(kk,i),v_hat(kk,i)].';
        
        phi_pre(i+1) = z_hat(1,i) + 1/z_hat(2,i)*z_hat(3,i)*dt*sind(z_hat(1,i)) + sqrt(Qw(1,1))*randn(1);
        rcomms_pre(i+1) = z_hat(2,i) - z_hat(3,i)*dt*cosd(z_hat(1,i)) + sqrt(Qw(2,2))*randn(1);
        v_pre(i+1) = z_hat(3,i) + sqrt(Qw(3,3))*randn(1);
        z_pre(:,i+1) = [phi_pre(i+1),rcomms_pre(i+1),v_pre(i+1)].';  

        % 测量方程
%         Nt(i+1) = antennanumber(5.5,rcomms_pre(i+1),phi_pre(i+1));  % 天线数 计算
        Nt(i+1) = 128;
      
%         if Nt(i)>128
%             Nt(i) = 128;
%         end
      
        a0 = 1/sqrt(Nt(i+1))*exp(-j*pi*[0:Nt(i+1)-1]'*cosd(phi_pre(i+1)));
        Kr = sqrt(Nt(i+1)*Nr);  % 阵列增益factor
        
        beta(i+1) = 1/(2*rcomms(i+1))^2;  % 路径损耗 与真实的rcomms有关，与预测值无关
%         beta(i) = 1/(2*rcomms_pre(i))^2;

        K = 1;
        for ii =1:K
            a(:,ii) = 1/sqrt(Nt(i+1))*exp(-j*pi*[0:Nt(i+1)-1]'*cosd(theta(ii,i+1)));
            %         SNR_d(ii) =  p*G*Kr^2*abs(beta(i))^2*abs(a(:,ii)'*a0)^2/sigma2;
            SNR_d(ii) =  p*G*abs(beta(i+1))^2*abs(a(:,ii)'*a0)^2/sigma2;
        end
        
        % 各散射点估计的方差
        var_d_point(:,i+1) = a2*3/( 2*pi^2*Nt(i+1)*Nr)./SNR_d'; % 4x1列向量
        var_a_point(:,i+1) = a1*1/( 2*pi^2*Nt(i+1)*Nr)./SNR_d'; % 4x1列向量
        var_u_point(:,i+1) = a3*1/( 2*pi^2*Nt(i+1)*Nr)./SNR_d'; % 4x1列向量

        % 各散射点的高分辨 
        r_est_point(:,i+1) = r(K,i+1) + sqrt(var_d_point(:,i+1)).*randn(K,1);      % 各散射点距离估计
        angle_point(:,i+1) = theta(K,i+1) + sqrt(var_a_point(:,i+1)).*randn(K,1);  % 各散射点角度估计
        
%         vx = 30;
        u_point(:,i+1) = 2*vx*cosd(theta(K,i+1))*fc./c;
        u_est_point(:,i+1) = u_point(:,i+1) + sqrt(var_u_point(:,i+1)).*randn(K,1);  % 各散射点角度估计
        
        clear a0 b a yy X   % 必须clear 因为每个时刻的Nt不一样，导致数据尺寸不同
                

        % 测量方程：通信Rx 
        xc_hat(i+1) = sum(r_est_point(:,i+1).*cosd(angle_point(:,i+1)))/K; % 质心坐标
        yc_hat(i+1) = sum(r_est_point(:,i+1).*sind(angle_point(:,i+1)))/K;
        x_comms_hat(i+1) = xc_hat(i+1); % comm Rx 初始坐标x
        y_comms_hat(i+1) = yc_hat(i+1); % comm Rx 初始坐标y
        
        phi_hat(kk,i+1) = arctand(x_comms_hat(i+1),y_comms_hat(i+1));  % 估计的通信Rx角度        
        range_hat(kk,i+1) = sqrt(x_comms_hat(i+1)^2+y_comms_hat(i+1)^2);
        v_hat(kk,i+1) = c/(2*fc)*sum(u_est_point(:,i+1).*cosd(angle_point(:,i+1))./var_u_point(:,i+1))./sum((cosd(angle_point(:,i+1))).^2./var_u_point(:,i+1)); % MLE 速度的最优估计
        y_m(:,i+1) = [phi_hat(kk,i+1),range_hat(kk,i+1),v_hat(kk,i+1)].';  % 测量方程

%         y_m(:,i+1) = y_measurement(kk,:,i+1);  % 测量方程

        %%  %%%%%%%%%%%%%%%%%%%%%% 2. 线性化 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Hn_1 = [ 1+ (z_hat(3,i)*dt*cosd(z_hat(1,i)))/(z_hat(2,i)), -(z_hat(3,i)*dt*sind(z_hat(1,i)))/(z_hat(2,i))^2, 0 ;
                 z_hat(3,i)*dt*sind(z_hat(1,i)), 1, -dt*cosd(z_hat(1,i));
                 0, 0, 1];
        %%  %%%%%%%%%%%%%%%%%%%%%% 3. MSE矩阵预测 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Mnln_1 = Hn_1*Mn*Hn_1'+Qw;
        %%  %%%%%%%%%%%%%%%%%%%%%% 4. 卡尔曼增益矩阵计算 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Kn = Mnln_1*inv(Mnln_1+Qz);
        %%  %%%%%%%%%%%%%%%%%%%%%% 5. 状态跟踪 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
        z_hat(:,i+1) = z_pre(:,i+1) + Kn*(y_m(:,i+1) - z_pre(:,i+1) ); % y_m: measurement
        %%  %%%%%%%%%%%%%%%%%%%%%% 6. MSE矩阵更新 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Mn = (eye(3)-Kn)*Mnln_1;

        % 计算achievable rate 
        pc = 1;
        aref = 1;
        an = aref./rcomms.*exp(j*2*pi*fc/c.*rcomms); % 信号衰减系数，与真实的rcomms有关，与预测值无关
        M = 1;  % 通信Rx天线数
        sigma2C = 1;  % 通信噪声功率        
        Kc = sqrt(Nt(i+1)*M);
        
        a_c = 1/sqrt(Nt(i+1))*exp(-j*pi*[0:Nt(i+1)-1]'*cosd(phi(i+1)));
        a0_c = 1/sqrt(Nt(i+1))*exp(-j*pi*[0:Nt(i+1)-1]'*cosd(phi_pre(i+1)));
        hn(i+1) = norm(an(i+1)*Kc*a_c'*a0_c).^2/sigma2C;
        Rn0_point(kk,i+1) = log2(1+pc*hn(i+1));
        
        clear a_c a0_c % 必须清空，很重要
               
    end
    
    Phi_pre(kk,:) = phi_pre;
    Rcomms_pre(kk,:) = rcomms_pre;
    
    rmse_phi(kk,:)= (phi_pre-phi).^2;
end

Rn_point = sum(Rn0_point)/KK;

Phi_ave = sum(Phi_pre)/KK;
Rcomms_ave = sum(Rcomms_pre)/KK;

Rmse_phi = sqrt(sum(rmse_phi)/KK);

%% 将点目标EKF的预测角和预测距离记录下来
phi_point = Phi_ave;
rcomms_point = Rcomms_ave;

Rmse_phi_point = Rmse_phi;



figure(1)
Rn_ave(1)=Rn_ave(2);
Rn_obj_ave(1)=Rn_obj_ave(3);Rn_obj_ave(2)=Rn_obj_ave(3);
Rn_point(1:2) = Rn_point(1:2);
% Rn_ave_scheme1(1)=Rn_ave_scheme1(3);Rn_ave_scheme1(2)=Rn_ave_scheme1(3);
plot([1:T/dt+1]*dt,Rn_obj_ave,'b','linewidth',2)
hold on
plot([1:T/dt+1]*dt,Rn_ave,'r--','linewidth',2)
hold on
Rn_ave1(1) = Rn_ave1(2);
plot([1:T/dt+1]*dt,Rn_ave1,'k-.','linewidth',2)

plot([2:T/dt]*dt,R_abp(2:end),'c-.','linewidth',2)
plot([1:T/dt]*dt,Rn_point(2:end),'m-.','linewidth',2)

xlabel('Time (s)','Fontname', 'Times New Roman','FontSize',14);
ylabel('Achievable rate','Fontname', 'Times New Roman','FontSize',14);xlim([0,T])
grid on
h = legend('Optimized Objective Function of Scheme 2','True Achievable Rate of Scheme 2','True Achievable Rate of Scheme 1','ABP-Comms-Only [11]','EKF-Point Target [5]')
set(h,'FontName','Times New Roman','FontSize',14,'FontWeight','normal')

xlim([0,T])


toc
