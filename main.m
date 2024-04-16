clear;clc

tic

%% ������
c=3e8;

dt = 0.05; % ����ʱ����
fc = 30e9; % ��Ƶ
lambda = c/fc;
p = 1; % ���书��
G = 10; % ƥ���˲�����

B = 500e6;

L = 5; % �������� 
W = 2; % �������
Dx = 1.5; % comms Rx���������x-axes�ľ���
Dy = 0.5; % comms Rx���������y-axes�ľ���

K = 8;

%% Ŀ��켣 ��ʼ״̬
vx = -20;

T = 160/abs(vx); % ����ʱ��2.5s

vy = 0;
R = sqrt(60^2+20^2);  % ���� ��ʼ����
theta0 = asind(20/R);  % ���� ��ʼ�Ƕ�
xc0 = R*cosd(theta0); % ���� ��ʼ����x
yc0 = R*sind(theta0); % ���� ��ʼ����y

[x(:,1),y(:,1),theta(:,1),r(:,1),x_comms(1),y_comms(1),phi(1),rcomms(1)] ...
    = scatterer(K,xc0,yc0,L,W,Dx,Dy,vx,vy);  % �ĸ���ɢ�������ͽǶ�

%% ɢ���Ŀ��״̬
xc(1)=xc0;
yc(1)=yc0;
for i = 2:T/dt+1
  
    xc(i) = xc(i-1) + vx*dt;  % ���� true �켣
    yc(i) = yc(i-1) + vy*dt;

    [x(:,i),y(:,i),theta(:,i),r(:,i),x_comms(i),y_comms(i),phi(i),rcomms(i)]...
        = scatterer(K,xc(i),yc(i),L,W,Dx,Dy,vx,vy);  % �ĸ���ɢ��� true �켣
    
end

%% Ԥ���㷨
phi_pre(1:3) = phi(1:3); % comms Rx��ʼ�Ƕ�Ԥ��
rcomms_pre(1:3) = rcomms(1:3);
xc_pre(1:3) = xc(1:3); yc_pre(1:3) = yc(1:3);
xc_hat(1:2) = xc(1:2); yc_hat(1:2) = yc(1:2);

Nr = 128; % RSU���ջ�������
sigma2 = 0.15; % echo��������
% a2 = 3.5e-2; % ����a0���ƾ���SNR
% a1 = 1.05e-2; % ����a1���ƽǶ�SNR
% a3 = 1e-2;
a2 = 2.0e-2; % ����a0���ƾ���SNR
a1 = 5.05e-3; % ����a1���ƽǶ�SNR
a3 = 1e-2;
KK = 250;

vx = abs(vx);

for kk =1:KK
    
    for i = 1:T/dt+1
        
        Nt(i) = antennanumber(5.5,rcomms(i),phi(i));  % ������ ����
        Nt(i) = min(Nt(i),128);
      
%         if Nt(i)>128
%             Nt(i) = 128;
%         end
        
        a0 = 1/sqrt(Nt(i))*exp(-j*pi*[0:Nt(i)-1]'*cosd(phi(i)));
        Kr = sqrt(Nt(i)*Nr);  % ��������factor
        
        beta(i) = 1/(2*rcomms(i))^2;
        
        for ii =1:K
            a(:,ii) = 1/sqrt(Nt(i))*exp(-j*pi*[0:Nt(i)-1]'*cosd(theta(ii,i)));
            SNR_d(ii) =  p*G*abs(beta(i))^2*abs(a(:,ii)'*a0)^2/sigma2;
        end
        
        % ��ɢ�����Ƶķ���
        var_d(:,i) = a2*3/( 2*pi^2*Nt(i)*Nr)./SNR_d'; % 4x1������ 
        var_a(:,i) = a1*1/( 2*pi^2*Nt(i)*Nr)./SNR_d'; % 4x1������
        var_u(:,i) = a3*1/( 2*pi^2*Nt(i)*Nr)./SNR_d'; % 4x1������
        
        %% %%%%%%%%%%%%%%%%%% ��ɢ���ĸ߷ֱ� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        r_est(:,i) = r(:,i) + sqrt(var_d(:,i)).*randn(K,1);      % ��ɢ���������
        angle(:,i) = theta(:,i) + sqrt(var_a(:,i)).*randn(K,1);  % ��ɢ���Ƕȹ���
        
%         vx = 20;
        u(:,i) = 2*vx*cosd(theta(:,i))*fc./c;
        u_est(:,i) = u(:,i) + sqrt(var_u(:,i)).*randn(K,1);  % ��ɢ���Ƕȹ���
        
        clear a0 b a yy X   % ����clear ��Ϊÿ��ʱ�̵�Nt��һ�����������ݳߴ粻ͬ
        
        %% %%%%%%%%%%%%%%%%%% �������̣�ͨ��Rx %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        xc_hat(i) = sum(r_est(:,i).*cosd(angle(:,i)))/K; % ��������
        yc_hat(i) = sum(r_est(:,i).*sind(angle(:,i)))/K;
        x_comms_hat(i) = xc_hat(i) + sqrt(Dx^2+Dy^2)*cosd( atand(Dy/Dx)+atand(vy/vx) ); % comm Rx ��ʼ����x
        y_comms_hat(i) = yc_hat(i) + sqrt(Dx^2+Dy^2)*sind( atand(Dy/Dx)+atand(vy/vx) ); % comm Rx ��ʼ����y
        
        phi_hat(kk,i) = arctand(x_comms_hat(i),y_comms_hat(i));  % ���Ƶ�ͨ��Rx�Ƕ�
        tt(i) =  arctand(x_comms_hat(i),y_comms_hat(i)) - arctand(x_comms(i),y_comms(i));
        sigma_n(kk,i) = abs( tt(i) );  % ͨ��Rx �Ƕȱ�׼��
        
        range_hat(kk,i) = sqrt(x_comms_hat(i)^2+y_comms_hat(i)^2);
        sigma_d(kk,i) = abs( sqrt(x_comms_hat(i)^2+y_comms_hat(i)^2) - sqrt(x_comms(i)^2+y_comms(i)^2)); % ͨ��Rx �����׼��

        v_hat(kk,i) = c/(2*fc)*sum(u_est(:,i).*cosd(angle(:,i))./var_u(:,i))./sum((cosd(angle(:,i))).^2./var_u(:,i)); % MLE �ٶȵ����Ź���
        sigma_v(kk,i) = abs(v_hat(kk,i)-vx);  % �����ٶȱ�׼��
        
        y_measurement(kk,:,i) = [phi_hat(kk,i),range_hat(kk,i),v_hat(kk,i)].';  % ��������

        %% %%%%%%%%%%%%%%%%%% Taylor��ʽ���� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         K=4; % 4��ɢ���
        Delta_x = sum(r_est(:,i).*cosd(angle(:,i))) + K*Dx;
        Delta_y = sum(r_est(:,i).*sind(angle(:,i))) + K*Dy;
        % �Ƕ�
        dg_r(:,i) = ( sind(angle(:,i)).*Delta_x - cosd(angle(:,i)).*Delta_y )...
            ./( Delta_x.^2 + Delta_y.^2 );
        dg_a(:,i) = ( r_est(:,i).*cosd(angle(:,i)).*Delta_x + r_est(:,i).*sind(angle(:,i)).*Delta_y )...
            ./( Delta_x.^2 + Delta_y.^2 );
        J = [dg_r(:,i);dg_a(:,i)].'; % Jacobian����
        Sigma = diag([var_d(:,i);var_a(:,i)]);
        sigma_n_taylor(kk,i) = sqrt(J*Sigma*J.'); % ͨ��Rx �Ƕȱ�׼�� Taylor����
        
        % ����
        df_r(:,i) = ( cosd(angle(:,i)).*Delta_x + sind(angle(:,i)).*Delta_y )...
            ./(K* sqrt(Delta_x.^2 + Delta_y.^2) );
        df_a(:,i) = ( -r_est(:,i).*sind(angle(:,i)).*Delta_x + r_est(:,i).*cosd(angle(:,i)).*Delta_y )...
            ./(K* sqrt(Delta_x.^2 + Delta_y.^2) );
        F = [df_r(:,i);df_a(:,i)].'; % Jacobian����
        sigma_n_taylor_d(kk,i) = sqrt(F*Sigma*F.'); % ͨ��Rx �Ƕȱ�׼�� Taylor����

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



%% ����2
% KK = 100;
Qw = diag([0.01^2,0.1^2,0.25^2]);
M0 = 128; % խ�����ķ���������
Kc_narrow = sqrt(M0*1);  % խ�����̶�32������
rho = linspace(0,1,10000); % grid search����

for kk =1:KK  %% �����Σ���achievable rate��ƽ�����
    
%     Qz = diag([sigma2_n(i),sigma2_d(i),sigma2_v(i)]);
    Mn = zeros(3,3);
    phi_hat(kk,1) = phi(1);
    range_hat(kk,1) = rcomms(1);
    v_hat(kk,1) = abs(vx);

    %% Ԥ���㷨 ����2
        
    rho0(1) = 0; % ��ʼ�� ȡrho = 1
    Nr = 128;  % ���������� �̶� ���ֵ
    for i = 1:T/dt
        
        Qz = diag([sigma2_n(i),sigma2_d(i),sigma2_v(i)]);

        % �����㷨 EKF
        %%  %%%%%%%%%%%%%%%%%%%%%% 1. ״̬Ԥ��  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        z_hat(:,i) = [phi_hat(kk,i),range_hat(kk,i),v_hat(kk,i)].';
        
        phi_pre(i+1) = z_hat(1,i) + 1/z_hat(2,i)*z_hat(3,i)*dt*sind(z_hat(1,i)) + sqrt(Qw(1,1))*randn(1);
        rcomms_pre(i+1) = z_hat(2,i) - z_hat(3,i)*dt*cosd(z_hat(1,i)) + sqrt(Qw(2,2))*randn(1);
        v_pre(i+1) = z_hat(3,i) + sqrt(Qw(3,3))*randn(1);
        z_pre(:,i+1) = [phi_pre(i+1),rcomms_pre(i+1),v_pre(i+1)].';  

        % ��������
        Nt(i+1) = antennanumber(5.5,rcomms_pre(i+1),phi_pre(i+1));  % ������ ����
        Nt(i) = min(Nt(i),128);
      
%         if Nt(i)>128
%             Nt(i) = 128;
%         end
        
       %% �Ż�������� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         p = 1;
        aref = 1;
        an(i+1) = aref./rcomms_pre(i).*exp(j*2*pi*fc/c.*rcomms_pre(i+1)); % LoS·����ʧ�����Ż�ʱδ֪����Ԥ��ֵ����
        Kc = sqrt(Nt(i+1)*1); % ͨ��Rx������

        sigma2C = 1;
        h_wide(i+1) = norm(an(i+1)*Kc).^2/sigma2C;
        R_wide(i+1) = log2(1+p*h_wide(i+1));
        
        %     a_narrow = 1/sqrt(M0)*exp(-j*pi*[0:M0-1]'*cosd(phi(i)));
        %     a0_narrow = 1/sqrt(M0)*exp(-j*pi*[0:M0-1]'*cosd(phi_est(i)));  % խ�����ù��ƽ�
        h_narrow(i+1) = norm(an(i+1)*Kc_narrow).^2/sigma2C;
        %     h_narrow(i) = norm(an(i)*Kc_narrow*a_narrow'*a0_narrow).^2/sigma2C;
        R_narrow(i+1) = log2(1+p*h_narrow(i+1));
        
%         delta(i) = 101/M0/sind(phi_pre(i)); % �벨����� խ���� �̶�M0������
        delta(i+1) = 101/M0/sind(phi_pre(i+1)); % �벨����� խ���� �̶�M0������  ������phi_hat(i)������phi_pre(i)

        delta_wide(i+1) = 101/Nt(i+1)/sind(phi_pre(i+1)); % �벨����� խ���� �̶�M0������
%         f = rho.*erf(delta_wide(i+1).*sqrt(rho)./sqrt(2)./sqrt(sigma2_n(i+1))).*R_wide(i+1) + (1-rho).*erf(delta(i+1).*sqrt(rho)./sqrt(2)./sqrt(sigma2_n(i+1))).*R_narrow(i+1);  % ע�����ﲻӦ�������ķ����comms Rx����
        f = rho.*R_wide(i+1) + (1-rho).*erf(delta(i+1).*sqrt(rho)./sqrt(2)./sqrt(sigma2_n(i+1))).*R_narrow(i+1);  % ע�����ﲻӦ�������ķ����comms Rx����
               
        
        [R_whole(i+1),b0(i+1)] = max(f);
                        
        rho0(i+1) = rho(b0(i+1));  %�ҵ����ֵ��λ��
%         rho0(i+1) = rho(b0(i+1))+0.01;  %�ҵ����ֵ��λ��

        PA(i) = erf(delta(i+1).*sqrt(rho0(i+1))./sqrt(2)./sqrt(sigma2_n(i+1)));

        a0 = 1/sqrt(Nt(i+1))*exp(-j*pi*[0:Nt(i+1)-1]'*cosd(phi_pre(i+1)));
        Kr = sqrt(Nt(i+1)*Nr);  % ��������factor
        
        beta(i+1) = 1/(2*rcomms(i+1))^2;  % ·����� ����ʵ��rcomms�йأ���Ԥ��ֵ�޹�
%         beta(i) = 1/(2*rcomms_pre(i))^2;

        for ii =1:K
            a(:,ii) = 1/sqrt(Nt(i+1))*exp(-j*pi*[0:Nt(i+1)-1]'*cosd(theta(ii,i+1)));
            %         SNR_d(ii) =  p*G*Kr^2*abs(beta(i))^2*abs(a(:,ii)'*a0)^2/sigma2;
            SNR_d(ii) =  rho0(i+1)*p*G*abs(beta(i+1))^2*abs(a(:,ii)'*a0)^2/sigma2;
        end
        
        % ��ɢ�����Ƶķ���
        var_d(:,i+1) = a2*3/( 2*pi^2*Nt(i+1)*Nr)./SNR_d'; % 4x1������
        var_a(:,i+1) = a1*1/( 2*pi^2*Nt(i+1)*Nr)./SNR_d'; % 4x1������
        var_u(:,i+1) = a3*1/( 2*pi^2*Nt(i+1)*Nr)./SNR_d'; % 4x1������

        % ��ɢ���ĸ߷ֱ� 
        r_est(:,i+1) = r(:,i+1) + sqrt(var_d(:,i+1)).*randn(K,1);      % ��ɢ���������
        angle(:,i+1) = theta(:,i+1) + sqrt(var_a(:,i+1)).*randn(K,1);  % ��ɢ���Ƕȹ���
        
%         vx = 30;
        u(:,i+1) = 2*vx*cosd(theta(:,i+1))*fc./c;
        u_est(:,i+1) = u(:,i+1) + sqrt(var_u(:,i+1)).*randn(K,1);  % ��ɢ���Ƕȹ���
        
        clear a0 b a yy X   % ����clear ��Ϊÿ��ʱ�̵�Nt��һ�����������ݳߴ粻ͬ
        
        % �������̣�ͨ��Rx 
        xc_hat(i+1) = sum(r_est(:,i+1).*cosd(angle(:,i+1)))/K; % ��������
        yc_hat(i+1) = sum(r_est(:,i+1).*sind(angle(:,i+1)))/K;
        x_comms_hat(i+1) = xc_hat(i+1) + sqrt(Dx^2+Dy^2)*cosd( atand(Dy/Dx)+atand(vy/vx) ); % comm Rx ��ʼ����x
        y_comms_hat(i+1) = yc_hat(i+1) + sqrt(Dx^2+Dy^2)*sind( atand(Dy/Dx)+atand(vy/vx) ); % comm Rx ��ʼ����y
        
        phi_hat(kk,i+1) = arctand(x_comms_hat(i+1),y_comms_hat(i+1));  % ���Ƶ�ͨ��Rx�Ƕ�        
        range_hat(kk,i+1) = sqrt(x_comms_hat(i+1)^2+y_comms_hat(i+1)^2);
        v_hat(kk,i+1) = c/(2*fc)*sum(u_est(:,i+1).*cosd(angle(:,i+1))./var_u(:,i+1))./sum((cosd(angle(:,i+1))).^2./var_u(:,i+1)); % MLE �ٶȵ����Ź���
        y_m(:,i+1) = [phi_hat(kk,i+1),range_hat(kk,i+1),v_hat(kk,i+1)].';  % ��������

%         y_m(:,i+1) = y_measurement(kk,:,i+1);  % ��������


        %%  %%%%%%%%%%%%%%%%%%%%%% 2. ���Ի� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Hn_1 = [ 1+ (z_hat(3,i)*dt*cosd(z_hat(1,i)))/(z_hat(2,i)), -(z_hat(3,i)*dt*sind(z_hat(1,i)))/(z_hat(2,i))^2, 0 ;
                 z_hat(3,i)*dt*sind(z_hat(1,i)), 1, -dt*cosd(z_hat(1,i));
                 0, 0, 1];
        %%  %%%%%%%%%%%%%%%%%%%%%% 3. MSE����Ԥ�� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Mnln_1 = Hn_1*Mn*Hn_1'+Qw;
        %%  %%%%%%%%%%%%%%%%%%%%%% 4. ���������������� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Kn = Mnln_1*inv(Mnln_1+Qz);
        %%  %%%%%%%%%%%%%%%%%%%%%% 5. ״̬���� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
        z_hat(:,i+1) = z_pre(:,i+1) + Kn*(y_m(:,i+1) - z_pre(:,i+1) ); % y_m: measurement
        %%  %%%%%%%%%%%%%%%%%%%%%% 6. MSE������� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Mn = (eye(3)-Kn)*Mnln_1;

        % ����ʵ�ʵ�achievable ate
        an(i+1) = aref./rcomms(i).*exp(j*2*pi*fc/c.*rcomms(i+1)); % LoS·����ʧ��������ʵֵʱ������ʵ��rcomms
        a_wide = 1/sqrt(Nt(i+1))*exp(-j*pi*[0:Nt(i+1)-1]'*cosd(phi(i+1)));
        a0_wide = 1/sqrt(Nt(i+1))*exp(-j*pi*[0:Nt(i+1)-1]'*cosd(phi_pre(i+1))); % ������Ԥ���
        h_wide(i+1) = norm(an(i+1)*Kc*a_wide'*a0_wide).^2/sigma2C;
%         h_wide(i+1) = norm(an(i+1)*Kc).^2/sigma2C;
        R_wide(i+1) = log2(1+p*h_wide(i+1));
        
%         M0 = 64;
        Kc_narrow = sqrt(M0*1);  % խ�����̶�32������
        a_narrow = 1/sqrt(M0)*exp(-j*pi*[0:M0-1]'*cosd(phi(i+1)));
        a0_narrow = 1/sqrt(M0)*exp(-j*pi*[0:M0-1]'*cosd(phi_hat(kk,i+1)));  % խ�����ù��ƽ�
%         h_narrow(i+1) = norm(an(i+1)*Kc_narrow).^2/sigma2C;
        h_narrow(i+1) = norm(an(i+1)*Kc_narrow*a_narrow'*a0_narrow).^2/sigma2C;
        R_narrow(i+1) = log2(1+p*h_narrow(i+1));
        
        R_sum0(i+1) = rho0(i+1)*R_wide(i+1)+(1-rho0(i+1))*R_narrow(i+1);
        
        clear a_c a0_c a_wide a0_wide a_narrow a0_narrow % !!!!! �������

    end
    
    Phi_pre(kk,:) = phi_pre;
    Rcomms_pre(kk,:) = rcomms_pre;
    V_pre(kk,:) = v_pre;
%     Rn(kk,:) = Rn0;

    R_sum(kk,:) = R_sum0;
    R_obj(kk,:) = R_whole; % �Ż���Ŀ�꺯��
    rho_whole(kk,:) = rho0;
    
    %% Predicted RMSE����
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


%% ������2��Ԥ��Ǻ�Ԥ������¼����
phi_scheme2 = Phi_ave;
rcomms_scheme2 = Rcomms_ave;
phi_hat_scheme2 = sum(phi_hat)/KK;
Rmse_phi_scheme2 = Rmse_phi;
    
Rmse_phi_scheme2_2 = Rmse_phi_hat;

%% ����1
for kk =1:KK;
    
    % 0. ��������
    Qw = diag([0.01^2,0.1^2,0.25^2]);
%     Qz = diag([sigma2_n(i),sigma2_d(i),sigma2_v(i)]);
    Mn = zeros(3,3);
    phi_hat(kk,1) = phi(1);
    range_hat(kk,1) = rcomms(1);
    v_hat(kk,1) = abs(vx);
    for i = 1:T/dt
        
        Qz = diag([sigma2_n(i),sigma2_d(i),sigma2_v(i)]);

        % �����㷨 EKF
        %%  %%%%%%%%%%%%%%%%%%%%%%1. ״̬Ԥ��  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        z_hat(:,i) = [phi_hat(kk,i),range_hat(kk,i),v_hat(kk,i)].';
        
        phi_pre(i+1) = z_hat(1,i) + 1/z_hat(2,i)*z_hat(3,i)*dt*sind(z_hat(1,i)) + sqrt(Qw(1,1))*randn(1);
        rcomms_pre(i+1) = z_hat(2,i) - z_hat(3,i)*dt*cosd(z_hat(1,i)) + sqrt(Qw(2,2))*randn(1);
        v_pre(i+1) = z_hat(3,i) + sqrt(Qw(3,3))*randn(1);
        z_pre(:,i+1) = [phi_pre(i+1),rcomms_pre(i+1),v_pre(i+1)].';  

        % ��������
        Nt(i+1) = antennanumber(5.5,rcomms_pre(i+1),phi_pre(i+1));  % ������ ����
        Nt(i) = min(Nt(i),128);
      
%         if Nt(i)>128
%             Nt(i) = 128;
%         end
      
        a0 = 1/sqrt(Nt(i+1))*exp(-j*pi*[0:Nt(i+1)-1]'*cosd(phi_pre(i+1)));
        Kr = sqrt(Nt(i+1)*Nr);  % ��������factor
        
        beta(i+1) = 1/(2*rcomms(i+1))^2;  % ·����� ����ʵ��rcomms�йأ���Ԥ��ֵ�޹�
%         beta(i) = 1/(2*rcomms_pre(i))^2;

        for ii =1:K
            a(:,ii) = 1/sqrt(Nt(i+1))*exp(-j*pi*[0:Nt(i+1)-1]'*cosd(theta(ii,i+1)));
            %         SNR_d(ii) =  p*G*Kr^2*abs(beta(i))^2*abs(a(:,ii)'*a0)^2/sigma2;
            SNR_d(ii) =  p*G*abs(beta(i+1))^2*abs(a(:,ii)'*a0)^2/sigma2;
        end
        
        % ��ɢ�����Ƶķ���
        var_d(:,i+1) = a2*3/( 2*pi^2*Nt(i+1)*Nr)./SNR_d'; % 4x1������
        var_a(:,i+1) = a1*1/( 2*pi^2*Nt(i+1)*Nr)./SNR_d'; % 4x1������
        var_u(:,i+1) = a3*1/( 2*pi^2*Nt(i+1)*Nr)./SNR_d'; % 4x1������

        % ��ɢ���ĸ߷ֱ� 
        r_est(:,i+1) = r(:,i+1) + sqrt(var_d(:,i+1)).*randn(K,1);      % ��ɢ���������
        angle(:,i+1) = theta(:,i+1) + sqrt(var_a(:,i+1)).*randn(K,1);  % ��ɢ���Ƕȹ���
        
%         vx = 30;
        u(:,i+1) = 2*vx*cosd(theta(:,i+1))*fc./c;
        u_est(:,i+1) = u(:,i+1) + sqrt(var_u(:,i+1)).*randn(K,1);  % ��ɢ���Ƕȹ���
        
        clear a0 b a yy X   % ����clear ��Ϊÿ��ʱ�̵�Nt��һ�����������ݳߴ粻ͬ
        
        % �������̣�ͨ��Rx 
        xc_hat(i+1) = sum(r_est(:,i+1).*cosd(angle(:,i+1)))/K; % ��������
        yc_hat(i+1) = sum(r_est(:,i+1).*sind(angle(:,i+1)))/K;
        x_comms_hat(i+1) = xc_hat(i+1) + sqrt(Dx^2+Dy^2)*cosd( atand(Dy/Dx)+atand(vy/vx) ); % comm Rx ��ʼ����x
        y_comms_hat(i+1) = yc_hat(i+1) + sqrt(Dx^2+Dy^2)*sind( atand(Dy/Dx)+atand(vy/vx) ); % comm Rx ��ʼ����y
        
        phi_hat(kk,i+1) = arctand(x_comms_hat(i+1),y_comms_hat(i+1));  % ���Ƶ�ͨ��Rx�Ƕ�        
        range_hat(kk,i+1) = sqrt(x_comms_hat(i+1)^2+y_comms_hat(i+1)^2);
        v_hat(kk,i+1) = c/(2*fc)*sum(u_est(:,i+1).*cosd(angle(:,i+1))./var_u(:,i+1))./sum((cosd(angle(:,i+1))).^2./var_u(:,i+1)); % MLE �ٶȵ����Ź���
        y_m(:,i+1) = [phi_hat(kk,i+1),range_hat(kk,i+1),v_hat(kk,i+1)].';  % ��������

%         y_m(:,i+1) = y_measurement(kk,:,i+1);  % ��������

        %%  %%%%%%%%%%%%%%%%%%%%%% 2. ���Ի� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Hn_1 = [ 1+ (z_hat(3,i)*dt*cosd(z_hat(1,i)))/(z_hat(2,i)), -(z_hat(3,i)*dt*sind(z_hat(1,i)))/(z_hat(2,i))^2, 0 ;
                 z_hat(3,i)*dt*sind(z_hat(1,i)), 1, -dt*cosd(z_hat(1,i));
                 0, 0, 1];
        %%  %%%%%%%%%%%%%%%%%%%%%% 3. MSE����Ԥ�� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Mnln_1 = Hn_1*Mn*Hn_1'+Qw;
        %%  %%%%%%%%%%%%%%%%%%%%%% 4. ���������������� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Kn = Mnln_1*inv(Mnln_1+Qz);
        %%  %%%%%%%%%%%%%%%%%%%%%% 5. ״̬���� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
        z_hat(:,i+1) = z_pre(:,i+1) + Kn*(y_m(:,i+1) - z_pre(:,i+1) ); % y_m: measurement
        %%  %%%%%%%%%%%%%%%%%%%%%% 6. MSE������� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Mn = (eye(3)-Kn)*Mnln_1;

        % ����achievable rate 
        pc = 1;
        aref = 1;
        an = aref./rcomms.*exp(j*2*pi*fc/c.*rcomms); % �ź�˥��ϵ��������ʵ��rcomms�йأ���Ԥ��ֵ�޹�
        M = 1;  % ͨ��Rx������
        sigma2C = 1;  % ͨ����������        
        Kc = sqrt(Nt(i+1)*M);
        
        a_c = 1/sqrt(Nt(i+1))*exp(-j*pi*[0:Nt(i+1)-1]'*cosd(phi(i+1)));
        a0_c = 1/sqrt(Nt(i+1))*exp(-1i*pi*[0:Nt(i+1)-1]'*cosd(phi_pre(i+1)));
        hn(i+1) = norm(an(i+1)*Kc*a_c'*a0_c).^2/sigma2C;
        Rn0(i+1) = log2(1+pc*hn(i+1));
        
        clear a_c a0_c % ������գ�����Ҫ
        
    end
    Phi_pre(kk,:) = phi_pre;
    Rcomms_pre(kk,:) = rcomms_pre;
    V_pre(kk,:) = v_pre;
    Rn(kk,:) = Rn0;
    
    %% Predicted RMSE����
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

Rn_ave_scheme1 = Rn_ave1; % ��¼����1������


%% ������1��Ԥ��Ǻ�Ԥ������¼����
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

%% ABP�㷨
for kk =1:KK  %% �����Σ���achievable rate��ƽ�����
    M0 =128;
    Nr = 128;  % ���������� �̶� ���ֵ
    phi0 = phi/180*pi;
    theta_hat(1) = phi0(1);
    for i = 2:T/dt

        D = pi/32;
        a_abp = 1/sqrt(M0)*exp(-j*pi*[0:M0-1]'*cos(phi0(i)));  
        fl_abp = 1/sqrt(M0)*exp(j*[0:M0-1]'*( -pi*cos(theta_hat(i-1))-D ) );  
        fr_abp = 1/sqrt(M0)*exp(j*[0:M0-1]'*( -pi*cos(theta_hat(i-1))+D ) );  
        
        % ͨ�Ž����ź�
        aref = 1;
        sigma2C = 1;
        an(i) = aref./rcomms(i).*exp(j*2*pi*fc/c.*rcomms(i)); % LoS·����ʧ��������ʵֵʱ������ʵ��rcomms
        y1 = an(i)*a_abp'*fl_abp + sqrt(sigma2C)*randn(1,1);
        y2 = an(i)*a_abp'*fr_abp + sqrt(sigma2C)*randn(1,1);
        
%         lambda = (norm(y1)^2-norm(y2)^2)/ (norm(y1)^2+norm(y2)^2);
        lambda = (norm(an(i)*a_abp'*fl_abp)^2-norm(an(i)*a_abp'*fr_abp)^2)/ (norm(an(i)*a_abp'*fl_abp)^2+norm(an(i)*a_abp'*fr_abp)^2);
        
        % �㷨
        theta_hat(i) = acos( -( -pi*cos(theta_hat(i-1)) - asin( (lambda*sin(D)-lambda*sqrt(1-lambda^2)*sin(D)*cos(D) )...
                                ./(sin(D)^2+lambda^2*cos(D)^2) ) )/pi );
        
        
        % ����ʵ�ʵ�achievable ate
        Kc_narrow = sqrt(M0*1);  % խ�����̶�32������
%         a_narrow = 1/sqrt(M0)*exp(-j*pi*[0:M0-1]'*cosd(phi(i+1)));
        a0_narrow = 1/sqrt(M0)*exp(-j*pi*[0:M0-1]'*cos(theta_hat(i)));  % խ�����ù��ƽ�
%         h_narrow(i+1) = norm(an(i+1)*Kc_narrow).^2/sigma2C;
        h_narrow(i) = norm(an(i)*Kc_narrow*a_abp'*a0_narrow).^2/sigma2C;
        R0_abp(kk,i) = log2(1+p*h_narrow(i));

        clear a_c a0_c a_wide a0_wide a_narrow a0_narrow % !!!!! �������

    end
    
    Phi_pre(kk,:) = theta_hat;

    rmse_phi_abp(kk,:) = (theta_hat/pi*180-phi(2:end)).^2;

end

Rmse_phi_abp = sqrt(sum(rmse_phi_abp)/KK);

R_abp = sum(R0_abp)/KK;

% plot([2:T/dt]*dt,R_abp(2:end),'k-.','linewidth',2)

Phi_ave = sum(Phi_pre)/KK;

%% ��ABP�ġ����ƽǡ���¼����
phi_abp = Phi_ave;
Rmse_phi_abp = Rmse_phi_abp;


clear Phi_pre

%% Liu .TWC ��Ŀ��EKF
for kk =1:KK;
    
    % 0. ��������
    Qw = diag([0.01^2,0.1^2,0.25^2]);
%     Qz = diag([sigma2_n(i),sigma2_d(i),sigma2_v(i)]);
    Mn = zeros(3,3);
    phi_hat(kk,1) = theta(1,1);
    range_hat(kk,1) = rcomms(1);
    v_hat(kk,1) = abs(vx);
    for i = 1:T/dt
                
        clear a var_d var_a var_u r_est angle u u_est

        Qz = diag([sigma2_n(i),sigma2_d(i),sigma2_v(i)]);

        % �����㷨 EKF
        %%  %%%%%%%%%%%%%%%%%%%%%%1. ״̬Ԥ��  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        z_hat(:,i) = [phi_hat(kk,i),range_hat(kk,i),v_hat(kk,i)].';
        
        phi_pre(i+1) = z_hat(1,i) + 1/z_hat(2,i)*z_hat(3,i)*dt*sind(z_hat(1,i)) + sqrt(Qw(1,1))*randn(1);
        rcomms_pre(i+1) = z_hat(2,i) - z_hat(3,i)*dt*cosd(z_hat(1,i)) + sqrt(Qw(2,2))*randn(1);
        v_pre(i+1) = z_hat(3,i) + sqrt(Qw(3,3))*randn(1);
        z_pre(:,i+1) = [phi_pre(i+1),rcomms_pre(i+1),v_pre(i+1)].';  

        % ��������
%         Nt(i+1) = antennanumber(5.5,rcomms_pre(i+1),phi_pre(i+1));  % ������ ����
        Nt(i+1) = 128;
      
%         if Nt(i)>128
%             Nt(i) = 128;
%         end
      
        a0 = 1/sqrt(Nt(i+1))*exp(-j*pi*[0:Nt(i+1)-1]'*cosd(phi_pre(i+1)));
        Kr = sqrt(Nt(i+1)*Nr);  % ��������factor
        
        beta(i+1) = 1/(2*rcomms(i+1))^2;  % ·����� ����ʵ��rcomms�йأ���Ԥ��ֵ�޹�
%         beta(i) = 1/(2*rcomms_pre(i))^2;

        K = 1;
        for ii =1:K
            a(:,ii) = 1/sqrt(Nt(i+1))*exp(-j*pi*[0:Nt(i+1)-1]'*cosd(theta(ii,i+1)));
            %         SNR_d(ii) =  p*G*Kr^2*abs(beta(i))^2*abs(a(:,ii)'*a0)^2/sigma2;
            SNR_d(ii) =  p*G*abs(beta(i+1))^2*abs(a(:,ii)'*a0)^2/sigma2;
        end
        
        % ��ɢ�����Ƶķ���
        var_d_point(:,i+1) = a2*3/( 2*pi^2*Nt(i+1)*Nr)./SNR_d'; % 4x1������
        var_a_point(:,i+1) = a1*1/( 2*pi^2*Nt(i+1)*Nr)./SNR_d'; % 4x1������
        var_u_point(:,i+1) = a3*1/( 2*pi^2*Nt(i+1)*Nr)./SNR_d'; % 4x1������

        % ��ɢ���ĸ߷ֱ� 
        r_est_point(:,i+1) = r(K,i+1) + sqrt(var_d_point(:,i+1)).*randn(K,1);      % ��ɢ���������
        angle_point(:,i+1) = theta(K,i+1) + sqrt(var_a_point(:,i+1)).*randn(K,1);  % ��ɢ���Ƕȹ���
        
%         vx = 30;
        u_point(:,i+1) = 2*vx*cosd(theta(K,i+1))*fc./c;
        u_est_point(:,i+1) = u_point(:,i+1) + sqrt(var_u_point(:,i+1)).*randn(K,1);  % ��ɢ���Ƕȹ���
        
        clear a0 b a yy X   % ����clear ��Ϊÿ��ʱ�̵�Nt��һ�����������ݳߴ粻ͬ
                

        % �������̣�ͨ��Rx 
        xc_hat(i+1) = sum(r_est_point(:,i+1).*cosd(angle_point(:,i+1)))/K; % ��������
        yc_hat(i+1) = sum(r_est_point(:,i+1).*sind(angle_point(:,i+1)))/K;
        x_comms_hat(i+1) = xc_hat(i+1); % comm Rx ��ʼ����x
        y_comms_hat(i+1) = yc_hat(i+1); % comm Rx ��ʼ����y
        
        phi_hat(kk,i+1) = arctand(x_comms_hat(i+1),y_comms_hat(i+1));  % ���Ƶ�ͨ��Rx�Ƕ�        
        range_hat(kk,i+1) = sqrt(x_comms_hat(i+1)^2+y_comms_hat(i+1)^2);
        v_hat(kk,i+1) = c/(2*fc)*sum(u_est_point(:,i+1).*cosd(angle_point(:,i+1))./var_u_point(:,i+1))./sum((cosd(angle_point(:,i+1))).^2./var_u_point(:,i+1)); % MLE �ٶȵ����Ź���
        y_m(:,i+1) = [phi_hat(kk,i+1),range_hat(kk,i+1),v_hat(kk,i+1)].';  % ��������

%         y_m(:,i+1) = y_measurement(kk,:,i+1);  % ��������

        %%  %%%%%%%%%%%%%%%%%%%%%% 2. ���Ի� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Hn_1 = [ 1+ (z_hat(3,i)*dt*cosd(z_hat(1,i)))/(z_hat(2,i)), -(z_hat(3,i)*dt*sind(z_hat(1,i)))/(z_hat(2,i))^2, 0 ;
                 z_hat(3,i)*dt*sind(z_hat(1,i)), 1, -dt*cosd(z_hat(1,i));
                 0, 0, 1];
        %%  %%%%%%%%%%%%%%%%%%%%%% 3. MSE����Ԥ�� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Mnln_1 = Hn_1*Mn*Hn_1'+Qw;
        %%  %%%%%%%%%%%%%%%%%%%%%% 4. ���������������� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Kn = Mnln_1*inv(Mnln_1+Qz);
        %%  %%%%%%%%%%%%%%%%%%%%%% 5. ״̬���� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
        z_hat(:,i+1) = z_pre(:,i+1) + Kn*(y_m(:,i+1) - z_pre(:,i+1) ); % y_m: measurement
        %%  %%%%%%%%%%%%%%%%%%%%%% 6. MSE������� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Mn = (eye(3)-Kn)*Mnln_1;

        % ����achievable rate 
        pc = 1;
        aref = 1;
        an = aref./rcomms.*exp(j*2*pi*fc/c.*rcomms); % �ź�˥��ϵ��������ʵ��rcomms�йأ���Ԥ��ֵ�޹�
        M = 1;  % ͨ��Rx������
        sigma2C = 1;  % ͨ����������        
        Kc = sqrt(Nt(i+1)*M);
        
        a_c = 1/sqrt(Nt(i+1))*exp(-j*pi*[0:Nt(i+1)-1]'*cosd(phi(i+1)));
        a0_c = 1/sqrt(Nt(i+1))*exp(-j*pi*[0:Nt(i+1)-1]'*cosd(phi_pre(i+1)));
        hn(i+1) = norm(an(i+1)*Kc*a_c'*a0_c).^2/sigma2C;
        Rn0_point(kk,i+1) = log2(1+pc*hn(i+1));
        
        clear a_c a0_c % ������գ�����Ҫ
               
    end
    
    Phi_pre(kk,:) = phi_pre;
    Rcomms_pre(kk,:) = rcomms_pre;
    
    rmse_phi(kk,:)= (phi_pre-phi).^2;
end

Rn_point = sum(Rn0_point)/KK;

Phi_ave = sum(Phi_pre)/KK;
Rcomms_ave = sum(Rcomms_pre)/KK;

Rmse_phi = sqrt(sum(rmse_phi)/KK);

%% ����Ŀ��EKF��Ԥ��Ǻ�Ԥ������¼����
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
