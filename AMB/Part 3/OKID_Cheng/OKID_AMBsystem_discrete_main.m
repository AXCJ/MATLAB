clc; clear all; close all;

% OKID_AMBsystem_discrete_XY % get XY-axis
% OKID_AMBsystem_discrete_Z % get Z-axis

load('AMB_discrete_parameter.mat')
[T, X, U, Y] = sim('Model_Following_Destrete', Sim_t, opts_sim, [ Sim_t', [r' r' r' r' r']]);

