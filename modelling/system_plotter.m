clear; close all;
addpath('../maglevFunctions');
load('params.mat');

%% General parameters
approximationType = 0; % 0=fast / 1=accurate

min_height = params.magnets.h/2;
max_height = .15;
height_steps = 128;

%% Solenoids - red
params.solenoids.N = 4; % just one...
params.solenoids.ri = 0;
params.solenoids.ro = 0; % ...and don't draw it

%% Permanent magnets - gray
params.magnets.N = 1;
params.magnets.R = 0;
params.magnets.I = 2798.8;

%% Levitating magnet - blue
lev_height = .0898;
params.levitatingmagnet.ri = params.magnets.ri;
params.levitatingmagnet.ro = params.magnets.ro;
params.levitatingmagnet.h = params.magnets.h;  
params.levitatingmagnet.nr = params.magnets.nr;
params.levitatingmagnet.nh = params.magnets.nh;
params.levitatingmagnet.nl = params.magnets.nl;
params.levitatingmagnet.m = params.magnets.m;
params.levitatingmagnet.I = -params.magnets.I;

%% Plot selected maglev system
Zrs = linspace(min_height,max_height,height_steps); % heights tested
Fzs = zeros(size(Zrs)); % array to store corresponding forces
x0 = zeros(12,1); x0(3) = lev_height;
sys = maglevSystem(x0, params, approximationType);

% waitbar settings
max = length(Zrs);
update_rate = 5; %seconds
h = waitbar(0);
t1 = tic; t3 = tic;
eta = 0;
eta_arr = zeros(1,10);

for j = 1:length(Zrs)
    %simulate system evolution with every height
    temp = sys.f([zeros(1,2),Zrs(j),zeros(1,9)]',zeros(params.solenoids.N,1));
    Fzs(j) = temp(9); %zm_dot should = 0

    % update times
    t2 = toc(t1); t4 = toc(t3);
    t1 = tic;
    curr = j;
    eta_arr(mod(curr,length(eta_arr))+1) = round((max-curr)*t2);
    if(t4>update_rate)
        t3=tic;
        eta = mean(eta_arr);
    end

    % visualize
    waitbar(curr/max,h,sprintf("computing ... %.3f%%\n" + ...
        "eta: %d min %02d sec",curr*100/max,floor(eta/60),round(mod(eta,60))))
end

figure(2);
clf; grid minor; plot(Zrs,Fzs); hold on;
yline(0,"k--"); hold off;