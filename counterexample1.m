#! /usr/bin/octave-cli -q

clear all;
close all;
clc;

if compare_versions(OCTAVE_VERSION(), '4.0.0', '<')
  error('octave version 4.0.0 or higher required; you have %s\n', OCTAVE_VERSION())
end

addpath('./utils');
pkg load interval;

%% Initialize inverse binary entropy function
ib=ibinent();

%% Parameter a
a=infsupdec('0.8');

%% Calculate the rate R = I(U;X) = I(V;Z) associated with a
R=binent(infsupdec('1/2')*a)-infsupdec('1/2')*binent(a);
%% Calculate mu = I(U;V) associated with a
mu=infsupdec('2')*R-a;

%% Calculate alpha = beta corresponding the the BSC with rate R
alpha=ib.binentinv_int(infsupdec('1')-R);
%% Check alpha
assert(subset(R,infsupdec('1')-binent(alpha)));

%% Calculate mu_b corresponding to two BSCs with crossover probability alpha
mu_b=infsupdec('1')-binent(star(alpha,alpha));

printf('Gap is at least %e\n', inf(mu-mu_b));
