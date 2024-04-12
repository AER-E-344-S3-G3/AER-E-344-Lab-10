% Schlieren and shadowgraph analysis script for AER E 344 lab 10.
% AER E 344 Spring 2024 - Section 3 Group 3
clear; clc; close all;

%% Constants
u = symunit;

figure_dir = "../Figures/";
data_filename = "pressure.csv";

num_taps = 15; % []
P_atm = 1.010; % [bar]
P_atm = double(separateUnits(unitConvert(P_atm*u.bar,u.psi))); % [psi]
gamma = 1.4; % []

%{
Each element of state_idx is a row index/picture number which correspond to
the following states:
    state_idx(1) = under-expanded flow
    state_idx(2) = 3rd critical condition
    state_idx(3) = over-expanded flow
    state_idx(4) = 2nd critical condition
    state_idx(5) = normal shock inside the nozzle
    state_idx(6) = 1st critical condition
NOTE: these indices are 1-based indices, i.e., index 1 corresponds to
picture/data point 0 and the first row in the data file.
%}
state_idx = [4,41,96,122,169,220];

% Pressure tap corresponding to the throat.
throat_idx = 5;

% Pressure tap immediately before the shock.
shock_idx = 10;

% Pressure taps are columns [3,17] in the data file.
pressure_tap_columns = 3:17;

% Only create graphs for states 2, 4, 5, and 6.
printed_states = [2,4,5,6];
plot_titles = [...
    "Under-Expanded Flow", ...
    "3rd Critical Condition", ...
    "Over-Expanded Flow", ...
    "2nd Critical Condition", ...
    "Normal Shock Inside the Nozzle", ...
    "1st Critical Condition"];

options = optimset('Display','off');

%% Equations
% Area ratio solved for 0
AbyAstar_eqn = @(A,Astar,M) (5 + M.^2).^3 ./ (6^3 .* M) - A ./ Astar; % []

M_2_normal_shock_eqn = @(M_1) sqrt((1 + (gamma - 1) / 2 * M_1^2) ...
    / (gamma * M_1^2 - (gamma - 1) / 2)); % []

Astar_eqn = @(M, A) A .* 6^3 .* M ./ (5 + M.^2).^3; % [length^2]

P_t_eqn = @(P, M) P ...
    .* (1 + (gamma - 1) ./ 2 ...
    .* M.^2).^(gamma ./ (gamma - 1)); % [pressure]

total_static_eqn = @(P_t, P, M) ...
    (1 + (gamma - 1) ./ 2 .* M.^2).^(gamma ./ (gamma - 1)) - P_t ./ P; % []

% P_2 / P_1 = (7 * M_1^2 - 1) / 6
P_1_normal_shock_eqn = @(P_2, M_1) 6 * P_2 / (7 * M_1^2 - 1); % [pressure]

%% Data Import
[downstream_dist,tunnel_areas] = getnozzleparams; % [in,in^2]

% The first pressure tap is upstream of the nozzle.
nozzle_dist = downstream_dist(2:end); % [in]
nozzle_area = downstream_dist(2:end); % [in^2]

Astar = tunnel_areas(5); % [in^2]

unzip(data_filename + ".zip");
data_file = readtable(data_filename);

%{
The data measured from the pressure taps.
Row 1 corresponds to state_idx(1), row 2 corresponds to state_idx(2), etc.
Column 1 corresponds to pressure tap 1, column 2 corresponds to pressure
tap 2, etc.
%}
P_ms_gauge = table2array( ...
    data_file(state_idx,pressure_tap_columns)); % [psi]

% Pressure tap 4 was broken/disconnected. Let's average the values from
% pressure tap 3 and 5.
P_ms_gauge(:,4) = (P_ms_gauge(:,3) + P_ms_gauge(:,5)) ./ 2; % [psi]

P_ms = P_ms_gauge + P_atm; % [psi]

% The first pressure tap is upstream of the nozzle.
P_ms_nozzle = P_ms(:,2:end); % [psi]

%% Data Processing
%
% Measured Mach Number and Total Pressure Calculations
%
P_mt = ones(length(state_idx),num_taps) .* P_t_eqn( ...
    P_ms(:,throat_idx),1); % [psi]

M_m = fsolve(@(M) total_static_eqn(P_mt, P_ms, M), ...
    ones(length(state_idx),num_taps) * 0.5,options); % []

M_m1 = M_m(5,shock_idx); % []
M_m2 = M_2_normal_shock_eqn(M_m1); % []

% This calculation should use the measured static pressure immediately
% downstream of the shock, but we didn't have a pressure tap exactly there,
% so, the tap at shock_idx + 1 is the best we can do.
P_mt(5,shock_idx + 1:end) = ones(1,length(P_mt(5,:)) - shock_idx) ...
    .* P_t_eqn(P_ms(5,shock_idx + 1),M_m2); % [psi]

P_mt_nozzle = P_mt(:,2:end); % [psi]

% Recalculate M_m now that we have corrected the P_mt to account for the
% normal shock.
M_m = fsolve(@(M) total_static_eqn(P_mt, P_ms, M), ...
    ones(length(state_idx),num_taps) * 0.5,options); % []

M_m_nozzle = M_m(:,2:end); % []

%
% Theoretical Mach Number Calculations
%

% Calculates the Mach number throughout the supersonic wind tunnel. We must
% use different initial values for taps before the throat and after the
% throat, since the AbyAstar_eqn returns two Mach numbers—a subsonic and
% supersonic value.
M_theory = cat(2, ...
    fsolve(@(M) AbyAstar_eqn(tunnel_areas(1:throat_idx),Astar,M), ...
        ones(1,throat_idx)*0.5,options), ...
    fsolve(@(M) AbyAstar_eqn(tunnel_areas(throat_idx + 1:end),Astar,M), ...
        ones(1,length(tunnel_areas) - throat_idx)*2,options)); % []
M_theory = ones(length(state_idx),num_taps) .* M_theory; % []

% Fix the mach number after the throat for state 6.
M_theory(6,:) = fsolve(@(M) AbyAstar_eqn(tunnel_areas,Astar,M), ...
    ones(1,num_taps) * 0.5,options); % []

% Fix the Mach number after the shock for state 5.
M_1_theory = M_theory(1,shock_idx); % []
M_2_theory = M_2_normal_shock_eqn(M_1_theory); % []

A_2star_theory = Astar_eqn(M_2_theory, ...
    tunnel_areas(shock_idx + 1)); % [in^2]

M_theory(5,shock_idx + 1:end) = ...
    fsolve(@(M) AbyAstar_eqn(tunnel_areas(shock_idx + 1:end), ...
        A_2star_theory,M), ...
    ones(1,length(tunnel_areas) - shock_idx)*0.5,options); % []

M_theory_nozzle = M_theory(:,2:end); % []

%
% Theoretical Pressure Calculations
%

% Calculate total pressure distribution for states 2, 5, and 6.
% In these states, we know the exit pressure is the ambient pressure.
P_t_theory = zeros(length(state_idx),num_taps); % [psi]
P_t_theory([2,5,6],:) = ones(3,num_taps) ...
    .* P_t_eqn(P_atm, M_theory([2,5,6],end)); % [psi]

% Correct the total pressure for state 5 upstream of the shock using:
% P_01 * Astar1 = P_02 * Astar2 (from Durbin's notes)
P_t_theory(5,1:shock_idx) = ones(1,shock_idx) ...
    .* P_t_theory(5,shock_idx + 1) * A_2star_theory / Astar; % [psi]

% Calculate total pressure distribution for state 4.
% In this state, there is a normal shock precisely at the exit of the
% nozzle. The pressure downstream of the normal shock is the ambient
% pressure and the Mach number precisely upstream of the normal shock is
% known.
P_e4_theory = P_1_normal_shock_eqn(P_atm,M_theory(4,end)); % [psi]
P_t_theory(4,:) = ones(1,num_taps) ...
    .* P_t_eqn(P_e4_theory,M_theory(4,end)); % [psi]

P_t_theory_nozzle = P_t_theory(:,2:end); % [psi]

% Calculate pressure distribution for states 2, 4, 5, and 6.
P_theory = zeros(length(state_idx),num_taps); % [psi]
P_theory([2,4,5,6],:) = fsolve(@(P) ...
    total_static_eqn(P_t_theory([2,4,5,6],:),P,M_theory([2,4,5,6],:)), ...
    ones(4,num_taps) * 10, options); % [psi]

P_theory_nozzle = P_theory(:,2:end); % [psi]

%% Plotting
% Plot measured pressure (static and total) as a function of distance along
% the nozzle axis for states 2, 4, 5, and 6.
for i = printed_states
    figure;
    plot(nozzle_dist,P_ms_nozzle(i,:),"LineWidth",2);
    hold on;
    plot(nozzle_dist,P_mt_nozzle(i,:),"LineWidth",2);
    hold off;
    fontname("Times New Roman");
    fontsize(12, "points");
    title_str = "Measured Pressure Distribution - " + plot_titles(i);
    title(title_str);
    xlabel("Distance Downstream of the Nozzle Throat [in]");
    ylabel("Pressure [psi]");
    legend("Measured Static Pressure","Measured Total Pressure", ...
        "Location","southwest");
    grid on;
    saveas(gcf, figure_dir + title_str + ".svg");
end

% Plot theoretically predicted pressure (static and total) as a function of
% distance along the nozzle axis for states 2, 4, 5, and 6.
for i = printed_states
    figure;
    plot(nozzle_dist,P_theory_nozzle(i,:),"LineWidth",2);
    hold on;
    plot(nozzle_dist,P_t_theory_nozzle(i,:),"LineWidth",2);
    hold off;
    fontname("Times New Roman");
    fontsize(12, "points");
    title_str = "Theoretical Pressure Distribution - " + plot_titles(i);
    title(title_str);
    xlabel("Distance Downstream of the Nozzle Throat [in]");
    ylabel("Pressure [psi]");
    legend("Theoretical Static Pressure","Theoretical Total Pressure", ...
        "Location","southwest");
    grid on;
    saveas(gcf, figure_dir + title_str + ".svg");
end

% Plot measured and predicted wall pressure distribution for states 2, 4,
% 5, and 6.
for i = printed_states
    figure;
    plot(nozzle_dist,P_ms_nozzle(i,:),"LineWidth",2);
    hold on;
    plot(nozzle_dist,P_mt_nozzle(i,:),"LineWidth",2);
    plot(nozzle_dist,P_theory_nozzle(i,:),"LineWidth",2);
    plot(nozzle_dist,P_t_theory_nozzle(i,:),"LineWidth",2);
    hold off;
    fontname("Times New Roman");
    fontsize(12, "points");
    title_str = "Measured vs. Theoretical Pressure Distribution - " ...
        + plot_titles(i);
    title(title_str);
    xlabel("Distance Downstream of the Nozzle Throat [in]");
    ylabel("Pressure [psi]");
    if i == 6
        legend("Measured Static Pressure","Measured Total Pressure", ...
        "Theoretical Static Pressure","Theoretical Total Pressure", ...
        "Location","southeast");
    else
        legend("Measured Static Pressure","Measured Total Pressure", ...
        "Theoretical Static Pressure","Theoretical Total Pressure", ...
        "Location","southwest");
    end
    grid on;
    saveas(gcf, figure_dir + title_str + ".svg");
end

% Plot measured and predicted Mach number as a function of distance along
% the nozzle axis for states 2, 4, 5, and 6.
for i = printed_states
    figure;
    plot(nozzle_dist,M_theory_nozzle(i,:),"LineWidth",2);
    hold on;
    plot(nozzle_dist,M_m_nozzle(i,:),"LineWidth",2);
    hold off;
    fontname("Times New Roman");
    fontsize(12, "points");
    title_str = "Measured vs. Theoretical Mach Number - " + plot_titles(i);
    title(title_str);
    xlabel("Distance Downstream of the Nozzle Throat [in]")
    ylabel("Mach Number [ ]");
    legend("Theoretical Mach Number","Measured Mach Number", ...
        "Location","northwest");
    grid on;
    saveas(gcf, figure_dir + title_str + ".svg");
end

%% Tables
tables = cell(1,length(state_idx));
for i = 1:length(state_idx)
    tables{i} = table;
    tables{i}.DownstreamDistance = downstream_dist';
    tables{i}.P_ms = P_ms(i,:)';
    tables{i}.P_theory = P_theory(i,:)';
    tables{i}.P_mt = P_mt(i,:)';
    tables{i}.P_t_theory = P_t_theory(i,:)';
    tables{i}.M_m = M_m(i,:)';
    tables{i}.M_theory = M_theory(i,:)';
    path = convertStringsToChars(figure_dir + i + "-" + plot_titles(i) ...
        + ".tex");
    if ismember(i,[1,3])
        table2latex(tables{i},path, ...
        {'$D_\text{throat}$ [\unit{in}]', ...
        '$P_\text{measured}$ [\unit{psi}]', ...
        '$P_\text{theory}$ [\unit{psi}]', ...
        '$P_{0,\text{measured}}$ [\unit{psi}]', ...
        '$P_{0,\text{theory}}$ [\unit{psi}]', ...
        '$M_\text{measured}$','$M_\text{theory}$'}, ...
        [2,3,3,3,3,3,3],[3,5]);
    else
        table2latex(tables{i},path, ...
        {'$D_\text{throat}$ [\unit{in}]', ...
        '$P_\text{measured}$ [\unit{psi}]', ...
        '$P_\text{theory}$ [\unit{psi}]', ...
        '$P_{0,\text{measured}}$ [\unit{psi}]', ...
        '$P_{0,\text{theory}}$ [\unit{psi}]', ...
        '$M_\text{measured}$','$M_\text{theory}$'}, ...
        [2,3,3,3,3,3,3],[]);
    end
end

%% Clean Up
delete(data_filename);
