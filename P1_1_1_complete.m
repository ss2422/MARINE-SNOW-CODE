%% P1_1_1 完整版 10s 有post-processing
clear all;
close all;
load('C:\Users\shans\Desktop\complete\P1_1_1_complete\P111.mat')
load('C:\Users\shans\Desktop\complete\P1_1_1_complete\P111_force.mat')

%% Initial condition
rou_p = 1083;
a = 0.00499785881924887*2; % diameter
Vp = (1/6)*pi*a^3;
Ap = (1/4)*pi*a^2;
g = 9.81;
h = 0.154880011; % thickness
rou_up = 976;
rou_down = 1025;
alpha = 18.74999866;
nu = 1.18*10^(-6);

%% m_p*(dV/dt) (i.e. F_sum)
P111_UZT = table2array(P111);
zp = -P111_UZT(:,2)+(h/2);
F_sum = -72*P111_force(:,2) + rou_p*Vp*g;
dV_dt = F_sum/(rou_p*Vp);
U = -P111_UZT(:,3);
t = P111_UZT(:,1);
% plot(P111_UZT(:,1),-P111_UZT(:,3))
[M(4),N(4)] = max(abs(U));
z_new = zp(N(4));
zp = zp-z_new;
[M(1),N(1)] = size(zp);
[M(2),N(2)] = min(abs(zp));
[M(3),N(3)] = min(abs(zp-h));
[M(5),N(5)] = max(U(N(3):end));
U(1:N(4),1) = M(4); % V1
% U(1:N(4),1) = 0.216173276364;
U((61789:end),1) = 0.143797635152066; % V2
% U((N(5)+N(3)):end,1) =  0.150398434088; % V2
% U = U-0.0374; % 0.1618 m/s 版本
U = U - 0.048;
U((61789:end),1) = U(61791,1)*sqrt(0.542056/0.425667)
t_in = t(N(2));
t_out = t(N(3));
% 45310开始变成0
dV_dt((1:23053),1) = 0;
dV_dt((61790:end),1) = 0;
dV_dt((46190:61789),1) = dV_dt(46189);
% dV_dt(1:45309,1) = dV_dt(1:45309,1)-0.224057;
% dV_dt(1:45309,1) = dV_dt(1:45309,1)-0.224057;
% dV_dt(45309:end,1) = -0.0636; % 58971 P111
% dV_dt(1:23053,1) = 0;
F_sum = dV_dt*(rou_p*Vp);
figure(1)
subplot(3,1,1)
plot(t,U)
xlabel('t')
ylabel('U')
subplot(3,1,2)
plot(zp,U)
xline(0,'--')
xline(h,'--')
subplot(3,1,3)
plot(zp,F_sum)
% plot(zp,U/M(4))
% plot(zp/h,F_sum/((rou_p-rou_up)*Vp*g))
t_rec = t(61790);

%% rou_f
rou_f(1:(N(2)),1) = rou_up;
rou_f(N(2)+1:(N(3)-1),1) = (rou_up+rou_down)/2 + ((rou_up-rou_down)/2)*erf(-alpha*(zp(N(2)+1:(N(3)-1))-(h/2)));
rou_f(N(3):M(1),1) = rou_down;
figure(2)
plot(zp,rou_f)
xline(0,'--')
xline(h,'--')
ylim([970 1030])
% plot(t,U)
% plot(t,U)
% xlabel('t')
% ylabel('U')
% xline(t_in,'--')
% xline(t_out,'--')

%% F_WB,F_D
F_WB = (rou_p-rou_f)*Vp*g;
% Ri = (rou_p-rou_f)*g*a./(rou_p.*U.^2);
Re = U*a/nu;
C_D = 0.4 + 24./Re + 6./(1 + Re.^(1/2));
% C_D = 0.67*C_D.*sqrt(Ri);
F_D = (1/2)*Ap*C_D.*rou_f.*abs(U).*U;

%% FA,FH
MA = (dV_dt)*a./(U.^2);
CA = 2.1 - (0.132*MA.^2)./(1+0.12*MA.^2);
CH = 0.48 + (0.52*MA.^3)./((1 + MA).^3);
FA = -Vp*(rou_f/2).*dV_dt;
% s(1) = 0;
% FH_old = 0;
% t(M(1)+1) = 0;
% for i = 1:M(1)
%     FH_term1(i) = -(3/2)*(pi*nu)^(0.5)*rou_f(i).*CH(i);
%     delta_t(i+1) = t(i+1) - t(i);
%     delta_t(1) = t(1);
%     FH_term2(i) = dV_dt(i);    
%     for j = 1:100
%         delta_s = delta_t(i)/100;
%         FH_term3_part(j) = (1./((t(i)-s(j)).^(1/2)))*delta_s;
%         s(j+1) = s(j)+delta_s;
%     end
%     s(1) = s(100);
%     FH_term3(i) = sum(FH_term3_part,"all");
%     FH_term3_part = zeros(1,100);
%     FH(i) = FH_term1(i).*FH_term2(i).*FH_term3(i) + FH_old;
%     FH_old = FH(i);
% end
% FH(1:23053,1) = 0;
% FH = 7.7000e-05;
t = t(1:M(1));
figure(3)
subplot(1,2,1)
plot(t,rou_f.*C_D)
xline(t_in,'--')
xline(t_out,'--')
subplot(1,2,2)
plot(t,F_WB,'r',t,-F_D,'k',t,FA,'b')

%% FS
FS = F_WB+FA-F_sum-F_D;
% FS(1:N(2)+1600) = 0;
% FS(56000:end) = 0;
figure(4)
plot(t,-FS) %-0.00036 % 
UZT = [t zp U];

%% Vc
Vc = (-FS)./((rou_up-rou_f)*g); % -0.00036
Vc(1:23053,1)=0;
figure(5)
% subplot(1,1,1)
plot(zp/h,F_WB./((rou_p-rou_up)*Vp*g),'r',zp/h,-F_D./((rou_p-rou_up)*Vp*g),'k',zp/h,(FA)./((rou_p-rou_up)*Vp*g),'b',zp/h,(-FS)./((rou_p-rou_up)*Vp*g),'m',zp/h,(F_sum)./((rou_p-rou_up)*Vp*g),'g')
xlim([-1.5 7])
xline(0,'--')
xline(1,'--')
legend('G-F_{b}','F_{D}','F_{A}','F_{S}','m_{p}dU/dt')
xlabel('z/h')
ylabel('F/(roup-rou1)*Vp*g')
% subplot(1,2,2)
% plot(t,Vc)
% xlabel('t')
% ylabel('Vc')
% xline(t_in,'--')
% xline(t_out,'--')
% xline(t_rec,'--')
% Vc = 2.342e-06*t - 4.812e-06 (0 ~ 0.0318)
% Vc = 3.542e-07*t^2-1.71623e-06*t + 2.4273e-06 (0.0318 ~ h)
% Vc = -7.5626e-08*t + 6.2719e-07 (h ~ )
% rou_f = 56.089*t + 860.3641 (0 ~ h)
% rou_f = 1025 (h ~ )
% rou_f*CD = 99.7854*t + 353.1017 (0 ~ h)
% rou_f*CD = -15.9914*t + 730.4761 (h ~ )
% 
% ode_dzp_dt_2 = H(:,2);

%% ODE new model

% % (1) 0 ~ zp(25122)
% options = odeset('RelTol',1e-5);
% h0 = 0;
% dh0 = 0.1575;
% Tspan = [t_in,t(25122)]; % Give a time stage (0,t_sum), t_sum=2400s.
% [T,H] = ode45(@hLinkedResTest,Tspan,[h0 dh0],options, Vp, Ap, g, rou_up, rou_p);
% t_new_1 = T;
% U_new_1 = H(:,2);
% 
% function dhdt = hLinkedResTest(t, h, Vp, Ap, g, rou_up, rou_p)
% dhdt = [h(2); -0.5*Ap*(99.7854*t + 353.1017)*h(2)*h(2)/(rou_p*Vp) + Vp*g*(rou_p-56.089*t - 860.3641)/(rou_p*Vp) + g*(rou_up-56.089*t - 860.3641).*(2.342e-06*t - 4.812e-06)/(rou_p*Vp)];
% end

% % (2) zp(25122) ~ h
% options = odeset('RelTol',1e-5);
% h0 = 0.0318;
% dh0 =  0.152010689077203; % 0.1550;
% Tspan = [t(25122),t_out]; % Give a time stage (0,t_sum), t_sum=2400s.
% [T,H] = ode45(@hLinkedResTest,Tspan,[h0 dh0],options, Vp, Ap, g, rou_up, rou_p);
% t_new_2 = T;
% U_new_2 = H(:,2);
% function dhdt = hLinkedResTest(t, h, Vp, Ap, g, rou_up, rou_p)
% dhdt = [h(2); -0.5*Ap*(99.7854*t + 353.1017)*h(2)*h(2)/(rou_p*Vp) + Vp*g*(rou_p-56.089*t - 860.3641)/(rou_p*Vp) + g*(rou_up-56.089*t - 860.3641).*(1.0773e-06*t.^2-5.5133e-06*t + 7.2866e-06)/(rou_p*Vp)];
% end
% % /(rou_p*Vp)

% % h ~ zp(34856)
% options = odeset('RelTol',1e-5);
% h0 = 0.1549;
% dh0 = 0.0902;
% Tspan = [t_out,3.1353]; 
% [T,H] = ode45(@hLinkedResTest,Tspan,[h0 dh0],options, Vp, Ap, g, rou_up, rou_p);
% t_new_3 = T;
% U_new_3 = H(:,2);
% function dhdt = hLinkedResTest(t, h, Vp, Ap, g, rou_up, rou_p)
% dhdt = [h(2); -0.5*Ap*(-15.9914*t + 730.4761)*h(2)*abs(h(2))/(rou_p*Vp) + Vp*g*(rou_p-1025)/(rou_p*Vp) + g*(rou_up-1025)*(1.6465-07*t - 1.0518e-07)/(rou_p*Vp)];
% end
% % /(rou_p*Vp)


% % h ~
% options = odeset('RelTol',1e-5);
% h0 = zp(N(3));
% dh0 = U(N(3));
% Tspan = [t_out,t_rec]; % Give a time stage (0,t_sum), t_sum=2400s.
% [T,H] = ode45(@hLinkedResTest,Tspan,[h0 dh0],options, Vp, Ap, g, rou_up, rou_p);
% t_new_4 = T;
% U_new_4 = H(:,2);
% function dhdt = hLinkedResTest(t, h, Vp, Ap, g, rou_up, rou_p)
% dhdt = [h(2); -0.5*Ap*(-15.9914*t + 730.4761)*h(2)*abs(h(2))/(rou_p*Vp) + Vp*g*(rou_p-1025)/(rou_p*Vp) + g*(rou_up-1025)*(-7.5626e-08*t + 6.2719e-07)/(rou_p*Vp)];
% end
% % /(rou_p*Vp)










%% ODE old model

% % 0 ~ h
% options = odeset('RelTol',1e-5);
% h0 = 0;
% dh0 = 0.1575;
% Tspan = [t_in,t_out]; % Give a time stage (0,t_sum), t_sum=2400s.
% [T,H] = ode45(@hLinkedResTest,Tspan,[h0 dh0],options, Vp, Ap, g, rou_up, rou_p);
% t_old_3 = T;
% U_old_3 = H(:,2);
% function dhdt = hLinkedResTest(t, h, Vp, Ap, g, rou_up, rou_p)
% dhdt = [h(2); -0.5*Ap*(99.7854*t + 353.1017)*h(2)*h(2)/(rou_p*Vp) + Vp*g*(rou_p-56.089*t - 860.3641)/(rou_p*Vp) + g*(rou_up-56.089*t - 860.3641).*(3.612438e-07)/(rou_p*Vp)];
% end


% h ~
options = odeset('RelTol',1e-5);
h0 = zp(N(3));
dh0 = U(N(3));
Tspan = [t_out,t_rec];
[T,H] = ode45(@hLinkedResTest,Tspan,[h0 dh0],options, Vp, Ap, g, rou_up, rou_p,t_out,t_rec);
t_old_4 = T;
U_old_4 = H(:,2);
function dhdt = hLinkedResTest(t, h, Vp, Ap, g, rou_up, rou_p,t_out,t_rec)
dhdt = [h(2); -0.5*Ap*(-15.9914*t + 730.4761)*h(2)*abs(h(2))/(rou_p*Vp) + Vp*g*(rou_p-1025)/(rou_p*Vp) + g*(rou_up-1025)*(3.612438e-07*exp((t_out-t)/(t_rec-2.0536)))/(rou_p*Vp)];
end





% load('C:\Users\shans\Desktop\ode_1_old.mat')
% load('C:\Users\shans\Desktop\ode_2_old.mat')
% ode_t_1_old = ode_1_old(:,1);
% ode_t_2_old = ode_2_old(:,1);
% ode_zp_1_old = ode_1_old(:,2);
% ode_zp_2_old = ode_2_old(:,2);
% ode_u_1_old = ode_1_old(:,3);
% ode_u_2_old = ode_2_old(:,3);
% 
% figure(6)
% plot([ode_t_1_old;ode_t_2_old],[ode_u_1_old;ode_u_2_old])
% title('Old model')
% xlabel('t')
% ylabel('U')
% xline(t_in,'--')
% xline(t_out,'--')
% xline(t_rec,'--')
% 
% figure(7)
% plot([ode_zp_1_old;ode_zp_2_old],[ode_u_1_old;ode_u_2_old])
% title('Old model')
% xlabel('zp')
% ylabel('U')
% xline(0,'--')
% xline(h,'--')

% %%
% for i = 1:24
%     test(i) = find(t==((i-1)*0.1+1.9));
% end