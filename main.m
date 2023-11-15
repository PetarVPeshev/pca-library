close all;
clear;
clc;

[laser, photoconductor] = get_default_config();
laser.plot_time_envelope(laser.pulse_half_pwr / 100);
