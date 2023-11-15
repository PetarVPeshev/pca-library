close all;
clear;
clc;

[laser, photoconductor] = get_default_config();
laser.plot_time_envelope(laser.pulse_half_pwr / 100);
eff = Efficiency(laser, photoconductor, eff_opt = 0.35);
