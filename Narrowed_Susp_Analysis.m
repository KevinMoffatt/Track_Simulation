% Narrowed Suspension Analysis
clc, clear
% ORIGIN SET AT LOWER ARM PIVOT
% Define Lengths of links
L = 13.37804588;            % Length of lower arm
U = 7.08464449;             % Length of upper arm when projected on xy plane
Uz = 7.18811192;            % distance into screen from lower arm projection on xy plane
B = 4 - 1*0.625;            % Length of bellcrank to pushrod mnt [in]
A = 2.5 - 0*0.1;            % Length of bellcrank to shock mnt [in]
gamma = 90;                 % angle b/w A and B
P = 16.25 + 2*0.25;         % Length of Pushrod [in]
% F = see below;            % Length of Fixed Frame (pivot of lower arms to pivot of bellcrank) [in]
R = 10;                     % Length from lower pivot to pushrod mount [in]
% Define initial position of lower arm
Ly1 = -1.06049506;
Lx1 = sqrt(L^2 - Ly1^2);
theta1 = atand(Ly1/Lx1);
% Define j as incremental change in lower arm outer end in y-dir
tilt = asind(2.7/18.5);
j = [-1.5:0.05:2]/cosd(tilt);
Ly = Ly1 + j;
Lx = sqrt(L^2 - Ly.^2);
Lz = j.*tand(tilt);
theta = atand(Ly./Lx);
% Define fixed pivot of bellcrank and length
Fx = 1.25;
Fy = 15.10305542;
F = sqrt(Fx^2 + Fy^2);
% Define fixed pivot of spring
Sx = Fx - 6.75; % Fx - 9; 
Sy = Fy - 0.625; % Fy - 0.625;
S_base = sqrt((Fx - Sx)^2 + (Fy - Sy)^2);
% Define Postion of pushrod
Ry = R.*sind(theta);
Rx = sqrt(R^2 - Ry.^2);
% Calculate initial Position of Bellcrank to pushrod mount (Worst Part)
% Transformation angle
    trans = atand(-Fx/Fy);
    Rx_trans = Rx.*cosd(trans) + Ry.*sind(trans);
    Ry_trans = -Rx.*sind(trans) + Ry.*cosd(trans);

    G = ((Rx_trans.^2)./(Ry_trans - F).^2) + 1;
    H = (R^2 - P^2 + B^2 - F^2)./(2.*(Ry_trans - F));
    W = ((F - H).^2) - B^2;
    T = ((2.*Rx_trans.*(F - H))./(Ry_trans - F));

    Bx_temp = ((-T + sqrt((T.^2) - 4.*G.*W))./(2*G));
    By_temp = (H - ((Rx_trans.*Bx_temp)./(Ry_trans - F)));

    Bx = Bx_temp.*cosd(-trans) + By_temp.*sind(-trans);
    By = -Bx_temp.*sind(-trans) + By_temp.*cosd(-trans);
% Define Angles
beta = atand((By - Ry)./(Rx - Bx)) + theta;  % b/w pushrod an lower arm
alpha = atand((By - Fy)./(Bx - Fx));         % b/w bellcrank and lower arm
phi = 180 - beta + theta - alpha;            % b/w bellcrank and pushrod
% Calculate Position of other arm on Bellcrank
Ax_temp = A.*cosd(alpha);
Ay_temp = A.*sind(alpha);
Ax = (Ax_temp.*cosd(-gamma) + Ay_temp.*sind(-gamma)) + Fx;
Ay = (-Ax_temp.*sind(-gamma) + Ay_temp.*cosd(-gamma)) + Fy;
% Calculate position of Upper arm
Hx = 4.75;           % X distance from loweer arm pivot to upper pivot
Hy = 7.04678257;     % Y distance from loweer arm pivot to upper pivot
D = 12.17084016;     % Upright distance between upper and lower arm
syms Uy Ux
for i = 1:length(j)
    [U_x(i,:), U_y(i,:)] = solve(D == sqrt((Ux - Lx(i))^2 + (Uy - Ly(i))^2 + (Uz - Lz(i))^2), U == sqrt((Ux - Hx)^2 + (Uy - Hy)^2));
    % U_y(i,:) = double(solve(D == sqrt((Lx(i) - Hx - sqrt(U - (Uy - Hy)^2))^2 + (Uy - Ly(i))^2 + (Uz - Lz(i))^2),Uy));
end
Uy = double(U_y(:,2)); % Uy = abs(U_y(:,2));
Ux = double(U_x(:,2)); % Ux = Hx + sqrt(U^2 - (Uy - Hy).^2);

% Find Spring Length and angle
S_length = sqrt(A.^2 + S_base.^2 - 2.*A.*S_base.*cosd(180 - alpha - gamma));
spring_angle = acosd((S_base^2 - S_length.^2 - A^2)/(-2.*S_base.*A));
Total_spring_travel = S_length(1) - S_length(end);
% Sanity check if Link lengths are all good
B_checker = sqrt((By - Fy).^2 + (Bx - Fx).^2)';
A_checker = sqrt((Ax - Fx).^2 + (Ay - Fy).^2)';
R_checker = sqrt(Rx.^2 +Ry.^2)';
U_checker = sqrt((Ux - Hx).^2 + (Uy - Hy).^2)';
D_checker = abs(sqrt((Lx - Ux').^2 + (Uy' - Ly).^2 + (Uz - Lz).^2));
for i = 1:length(j)
    P_checker(i,:) = sqrt((Rx(i) - Bx(i)).^2 + (Ry(i) - By(i)).^2);
end
%% Plot Results
figure(1)
plot(Lx,Ly,'r',Rx,Ry,'b',Bx,By,'g',Ax,Ay,'cy',Ux,Uy,'y',[0 Fx Sx Hx],[0 Fy Sy Hy],'k*')
grid on
axis([-10 20 -8 20])
figure(2)
for i = 1:length(j)
    grid on
    axis([-10 20 -8 20])
    line([Lx(i) 0],[Ly(i) 0],'color','r')
    line([Rx(i) 0],[Ry(i) 0],'color','b')
    line([Rx(i) Bx(i)],[Ry(i) By(i)],'color','k')
    hold on
    line([Bx(i) Fx],[By(i) Fy],'color','g')
    line([Ax(i) Fx],[Ay(i) Fy],'color','cy')
    line([Ax(i) Sx],[Ay(i) Sy],'color','m')
    line([Ux(i) Hx],[Uy(i) Hy],'color','y')
    
end
Y = j'; Ly = Ly'; Ry = Ry'; beta = beta'; theta = theta'; phi = phi'; 
alpha = alpha'; spring_angle = spring_angle'; S_length = S_length';
table(Y,beta,alpha,phi,spring_angle,S_length)

% Other Stuff
% Motion Ratio
for i = 2:length(j)
    MR(i) = (S_length(i-1) - S_length(i))./(j(i) - j(i-1));
end
figure(3)
subplot(2,1,1)
plot(j(2:end),MR(2:end))
title('Motion Ratio: ((spring disp)/(Wheel Disp))'); xlabel('Wheel Disp [in]'); ylabel('MR');
grid on
subplot(2,1,2)
plot(j, spring_angle, 'r', j, phi, 'b')
grid on
title('Spring and Pushrod Angle vs Wheel Disp'); xlabel('Wheel Disp [in]'); ylabel('angle [deg]');
legend('Spring Angle (deg)','Pushrod angle','location','southeast')
% Force Calcs
Weight = 2000;          % lbs
F_w = Weight/4;         % Static force on each wheel
for i = 1:length(j)
    F_lp(i) = F_w*cosd(theta(i));
    F_r(i) = (R/L)*F_lp(i);
    F_p(i) = F_r(i)*cosd(90 - beta(i));
    F_b(i) = F_p(i)*cosd(phi(i) - 90);
    F_a(i) = (A/B)*F_b(i);
    F_k(:,i) = F_a(i)*cosd(spring_angle(i) - 90);
    delta_x(i) = S_length(1) - S_length(i);
end
j_static = find(j == 0);
% Spring Constant for a given ride height based on j_static
K_stat = F_k(j_static(1))/delta_x(j_static(1));
Spring_force = K_stat.*delta_x';
figure(4)
plot(j,Spring_force, j,F_k)
title('Spring Force and Applied Force vs Wheel Travel'); xlabel('Wheel Travel [in]'); ylabel('Force [lbf]');
grid on
Spring_Length = S_length(1);
table(Spring_Length,Total_spring_travel, K_stat)
% Solve for camber
theta_s = -1;
tau = atand((Ux - Lx')./(Uy - Ly));
camber = (tau - tau(j_static(1)) + theta_s);
Ox = linspace(0,0); % Defining set of zeros for axis
travel = (j.*cosd(tilt))';
figure(5)
plot(camber,travel,'r-',camber,Ox(1:length(camber)),'k-',Ox(1:length(travel)),travel,'k-','Linewidth',2)
grid on
ylabel('Suspension Travel [in]');
xlabel('Camber Angle');
title('Suspension Travel vs Camber Angle');
table(travel,camber)
% Find Roll center
% Set orgigin at center of tire on ground
d1 = 8.11;          % Height of lower arm pivot
d2 = Hy + d1;       % Height of upper pivot
track1 = 55.4;      % track width
ww = 225/25.4;      % wheel width
Pu = 18.75;         % Distance from car centerline to upper arm mount
PL = Pu - 4.75;     % Distance from car centerline to lower arm mount
% Adjust solved for variables to new coords
Lx_a1 = (track1/2) - (Lx(j_static(1)) + PL); Ly_a = Ly + d1;
Ux_a1 = (track1/2) - (Ux(j_static(1)) + Pu); Uy_a = Uy + d1;   
ang = atand(Lx_a1/Ly_a(j_static(1))) + (camber - theta_s);
dist = sqrt(Lx_a1^2 + Ly_a(j_static(1))^2);
track_x = track1 + 2*(dist*sind(ang));
track_y = (Ly_a - dist*cosd(ang));
Lx_a = (track_x/2) - (Lx' + PL);
Ux_a = (track_x/2) - (Ux + PL);

x = linspace(0,1.2*track1,length(j));
%index1 = zeros(length(j));
index_check = ones(length(j),1);
for i = 1:length(j)
    Ueqt(:,i) = ((d2 - Uy_a(i))/((track_x(i)/2) - Pu - Ux_a(i))).*x + (Uy_a(i) - ((d2 - Uy_a(i))/((track_x(i)/2) - Pu - Ux_a(i))).*Ux_a(i));
    Leqt(:,i) = ((d1 - Ly_a(i))/((track_x(i)/2) - PL - Lx_a(i))).*x + (Ly_a(i) - ((d1 - Ly_a(i))/((track_x(i)/2) - PL - Lx_a(i))).*Lx_a(i));
    index1 = find(Ueqt(:,i) < Leqt(:,i));
    index_check(1:length(index1),i) = index1;
    IC_y(i) = Ueqt(index_check(1,i),i);
    RCeqt(:,i) = (IC_y(i)/x(index_check(1,i))).*x;
    index2 = find(x > (track_x(i)/2));
    RC(i) = RCeqt(index2(1),i);
end
figure(6)
plot(x,Ueqt(:,j_static(1)),x,Leqt(:,j_static(1)),x,RCeqt(:,j_static(1)))
axis([-10 70 -5 20])
grid on
hold on
plot(Ux_a(j_static(1)),Uy_a(j_static(1)),'k*',Lx_a(j_static(1)),Ly_a(j_static(1)),'k*',(track1/2) - Pu,d2,'k*',(track1/2) - PL,d1,'k*')
legend('upper','lower')
hold on
plot(track1/2*ones(length(x),1),linspace(-30,30,length(x)), linspace(0,track1,length(x)), zeros(length(x),1))