%% PLOT RESONANT CIRCLES ON THE B-PLANE

% Generate locus of deflections
RE = 1; RE_km = 6378.140; GE = 398600.43543609593;
th = linspace(0,2*pi,100);
xi_defl = 3*RE*cos(th);
zeta_defl = 3*RE*sin(th);

% Read resonance file
resonances_file = ...
    '2017PDC_resonances_all.dat';
startRow = 4;
fileID = fopen(resonances_file,'r');
dataCell = textscan(fileID, '%3f%4f%15f%f%[^\n\r]', 'Delimiter', '', 'WhiteSpace', '',...
    'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
resonances = struct();
resonances.k = dataCell{1}; resonances.h = dataCell{2};
resonances.D = dataCell{3}; resonances.R = dataCell{4};

% Compute gravitational focussing
% xi_enc   = -8.681064454743289E-01;
% zeta_enc = -1.469758601677242E-01;
xi_enc = -0.8682; % Javi's value
zeta_enc = 0.141; % Javi's value
U_kms    = 1.320963421118330E+01;
focus_factor = sqrt(1 + 2*GE/(RE_km*U_kms^2));
RE_focussed = focus_factor;

% Settings for plotting
%max_res = 49;
max_res = size(resonances.k,1);
k_min = 1; k_max = 20;
xi_lim = 5; zeta_lim = 5;
lw = 1;


%% Plot the resonant circles

% Set figure and axes properties
figure('WindowStyle','normal');
xlabel('\xi (R_\oplus)');
ylabel('\zeta (R_\oplus)');
xlim([-xi_lim,xi_lim]); set(gca,'XTick',(-xi_lim:1:xi_lim));
ylim([-zeta_lim,zeta_lim]); set(gca,'YTick',(-zeta_lim:1:zeta_lim));
axis square
grid on
hold on

% Plot deflection circle
plot(xi_defl,zeta_defl,'-.','LineWidth',1,'Color','black')

% co = vangogh(22);
co = winter(22);
co = co(2:end,:);
for i=1:max_res
    % Only plot resonances between kmin and kmax
    if ((resonances.k(i) < k_min) || (resonances.k(i) > k_max))
        continue
    end
    curr_color = co(resonances.k(i),:);
    
    % Generate circle
    xi_circ   = resonances.R(i)*cos(th);
    zeta_circ = resonances.D(i) + resonances.R(i)*sin(th);
    plot(xi_circ,zeta_circ,'-','LineWidth',lw,'Color',curr_color);
    xlim([-xi_lim,xi_lim]);
    ylim([-zeta_lim,zeta_lim]);
    axis square
    
%     resonances.k(i)
%     resonances.h(i)
%     keyboard
    
end
colormap(co);
caxis([resonances.k(1),resonances.k(end)]);
cb = colorbar;
cb.Label.String = 'k';


% Plot Earth (1RE and focussed)
fill(RE_focussed*cos(th),RE_focussed*sin(th),'white');
plot(RE_focussed*cos(th),RE_focussed*sin(th),'LineWidth',lw,'color','black');
plot(cos(th),sin(th),'--','LineWidth',lw,'color','black');
% Plot b-plane crossing point
plot(xi_enc,zeta_enc,'Marker','x','Color','red','MarkerSize',10,'LineWidth',1)
% Plot deflection
load('bplane_3ER.mat')
plot(xi_opt_defl,zeta_opt_defl,'Marker','o','Color','red','MarkerSize',3,...
    'LineStyle','none')

set(gca,'FontSize',14);

%% Keyhole computation
close all
% Nominal Öpik parameters
U = 0.44705634713083897;
theta = 6.573081547432976E+01*pi/180;
phi = 7.141024642502616E+01*pi/180;
m = 3.0034896149157654e-06;
t0 = 0;
DU = 152003856.97586098;
RE_au = RE_km/DU;

% Additional check for the computation of next-impact probability for
% points in the upper left quadrant
zeta_segm_up = zeta_opt_defl(1) + 0.7130; 
zeta_segm_dw = zeta_opt_defl(1) - 0.7130;
xi_segm = xi_opt_defl(1);

lw=1;
figure;
hold on
kh_up_width = zeros(max_res,2);
kh_down_width = zeros(max_res,2);
kh_width_sum = 0;
for i=1:max_res
    
    % Compute upper and lower keyholes for each of the resonant circles
    [kh_up_xi,kh_up_zeta,kh_down_xi,kh_down_zeta] = ...
        two_keyholes(resonances.k(i),resonances.h(i),...
        resonances.D(i),resonances.R(i),U,theta,phi,m,t0);
    
    % Plot keyhole
    kh_up_xi = kh_up_xi./RE_au; kh_down_xi = kh_down_xi./RE_au;
    kh_up_zeta = kh_up_zeta./RE_au; kh_down_zeta = kh_down_zeta./RE_au;
%     plot(kh_down_xi(:,1),kh_down_zeta(:,1),kh_down_xi(:,1),kh_down_zeta(:,2),...
%         'Color','red','LineStyle','none','Marker','*');
%     plot(kh_up_xi(:,1),kh_up_zeta(:,1),kh_up_xi(:,1),kh_up_zeta(:,2),...
%         'Color','red','LineStyle','none','Marker','*');
    curr_color = co(resonances.k(i),:);
    % Plot top and bottom keyholes
    plot(kh_down_xi(:,1),kh_down_zeta(:,1),kh_down_xi(:,2),kh_down_zeta(:,2),...
        '-','LineWidth',lw,'Color',curr_color);
    plot(kh_up_xi(:,1),kh_up_zeta(:,1),kh_up_xi(:,2),kh_up_zeta(:,2),...
         '-','LineWidth',lw,'Color',curr_color);
    xlim([-xi_lim,xi_lim]);
    ylim([-zeta_lim,zeta_lim]);
    axis square
    
    % Compute & save keyhole widths
    kh_up_width(i,:) = [resonances.k(i),max(abs(kh_up_zeta(:,2) - kh_up_zeta(:,1)))];
    kh_down_width(i,:) = [resonances.k(i),max(abs(kh_down_zeta(:,2) - kh_down_zeta(:,1)))];
    
    % Check for intersection with the segment representing the dispersion
    % in zeta
    zeta_kh = sqrt(resonances.R(i)^2 - xi_segm^2) + resonances.D(i);
    dz = [zeta_segm_up,zeta_segm_dw] - zeta_kh;
    if (dz(1)*dz(2) < 0)
        % Add keyhole width to sum
        kh_width_sum = kh_width_sum + kh_up_width(i,2);
    end
    
end
colormap(co);
caxis([resonances.k(1),resonances.k(end)]);
cb = colorbar;
cb.Label.String = 'k';
% Plot Earth (1RE and focussed)
fill(RE_focussed*cos(th),RE_focussed*sin(th),'white');
plot(RE_focussed*cos(th),RE_focussed*sin(th),'LineWidth',1,'color','black');
xlabel('\xi (R_\oplus)');
ylabel('\zeta (R_\oplus)');
xlim([-xi_lim,xi_lim]); set(gca,'XTick',(-xi_lim:1:xi_lim));
ylim([-zeta_lim,zeta_lim]); set(gca,'YTick',(-zeta_lim:1:zeta_lim));
axis square
grid on
set(gca,'FontSize',14);
plot(cos(th),sin(th),'--','LineWidth',1,'color','black');
% Plot b-plane crossing point
plot(xi_enc,zeta_enc,'Marker','x','Color','red','MarkerSize',10,'LineWidth',1)
% Plot deflection
load('bplane_3ER.mat')
plot(xi_opt_defl,zeta_opt_defl,'Marker','o','Color','red','MarkerSize',3,...
    'LineStyle','none')
plot(xi_defl,zeta_defl,'-.','LineWidth',1,'Color','black')

% Display final probability
p_nextimpact = kh_width_sum/(zeta_segm_up - zeta_segm_dw)

%% Plot keyhole max. widths as a function of k
figure;
bluedef = [0, 0.4470, 0.7410];
hold on
plot(kh_up_width(:,1),kh_up_width(:,2)*RE_km,'Color',bluedef,'Marker','o',...
    'MarkerFaceColor',bluedef,'LineStyle','none','MarkerSize',3);
plot(kh_down_width(:,1),kh_down_width(:,2)*RE_km,'Color',bluedef,'Marker','o',...
    'MarkerFaceColor',bluedef,'LineStyle','none','MarkerSize',3);
xlabel('k');
ylabel('Max. width (km)');
grid on
set(gca,'FontSize',14)