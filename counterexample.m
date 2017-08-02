#! /usr/bin/octave-cli -q

clear all;
close all;

if compare_versions(OCTAVE_VERSION(), '4.0.0', '<')
  error('octave version 4.0.0 or higher required; you have %s\n', OCTAVE_VERSION())
end

addpath('./utils');
pkg load interval;
pkg load symbolic;

load('numeric_distribution.mat');

%% Grid parameters
K=800;
Q=3;
bdp=0.1;

%% Initialize inverse binary entropy function
ib=ibinent();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Obtain a symbolic version of P(UV|XZ) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Floating point P(UV|XZ)
P_UVgXZ=P_lgap./P_XZ;

%% Obtain symbolic distributions in cell arrays
Ps_UVgXZ={};
Ps_XZ={};
for x=1:2
  for z=1:2
    for u=1:2
      for v=1:2
	Ps_UVgXZ{u,v,x,z}=sym(num2str(P_UVgXZ(u,v,x,z),16));
	Ps_XZ{u,v,x,z}=sym(num2str(P_XZ(u,v,x,z),5));
      end
    end
  end
end

%% as is the parameter of the BSCs X -> U and Z -> V
as=Ps_UVgXZ{1,1,2,2}+Ps_UVgXZ{1,2,2,2};

%% ps is the paramter of the BSC X -> Z
ps=sym(num2str(p,5));

%% Construct a distribution satisfying the long Markov chain U - X - Z - V
%% where as is the parameter of the BSCs X -> U and Z -> V
Ps_UVgXZ_ref={};
for x=1:2
  for z=1:2
    for u=1:2
      for v=1:2
	Ps_UVgXZ_ref{u,v,x,z}=sym('1');
	if x==u
	  Ps_UVgXZ_ref{u,v,x,z}=Ps_UVgXZ_ref{u,v,x,z}*(sym('1')-as);
	else
	  Ps_UVgXZ_ref{u,v,x,z}=Ps_UVgXZ_ref{u,v,x,z}*as;
	end
	if z==v
	  Ps_UVgXZ_ref{u,v,x,z}=Ps_UVgXZ_ref{u,v,x,z}*(sym('1')-as);
	else
	  Ps_UVgXZ_ref{u,v,x,z}=Ps_UVgXZ_ref{u,v,x,z}*as;
	end
      end
    end
  end
end

%% Turn the cell arrays into vectors
Psv_UVgXZ=cell2sym(reshape(Ps_UVgXZ,[],1));
Psv_UVgXZ_ref=cell2sym(reshape(Ps_UVgXZ_ref,[],1));
Psv_XZ=cell2sym(reshape(Ps_XZ,[],1));

%% Obtain an exact symbolic solution by projecting onto the null space
%% And get the full distribution by multiplying with P(X,Z)
Psv_UVXZ=( Psv_UVgXZ_ref+get_projection_matrix()*(Psv_UVgXZ-Psv_UVgXZ_ref) ) .* Psv_XZ;

%% turn Psv_UVXZ into 2x2x2x2 cell array
Psc_UVXZ=reshape(num2cell(Psv_UVXZ),2,2,2,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct interval representations %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Rate interval
Ri=1-binent(infsupdec(char(as)));

%% BSC parameters
ppi=infsupdec(char(ps));

%% Interval representation of P(U,V,X,Z)
Piv_UVXZ=infsupdec(zeros(16,1));
for i=1:16
  Piv_UVXZ(i)=infsupdec(char(Psv_UVXZ(i)));
end

%% Construct P(U,V) in interval form
Pic_UVXZ={};
for x=1:2
  for z=1:2
    Pic_UVXZ{x,z}=infsupdec(zeros(2,2));
    for u=1:2
      for v=1:2
	Pic_UVXZ{x,z}(u,v)=infsupdec(char(Psc_UVXZ{u,v,x,z}));
      end
    end
  end
end
Pi_UV=Pic_UVXZ{1,1}+Pic_UVXZ{1,2}+Pic_UVXZ{2,1}+Pic_UVXZ{2,2};

%% Calculate the entropy of (XZ) as 1 + h(p)
Hi_XZ=infsupdec('1')+binent(ppi);
%% Calculate the entropy of (UXZV)
Hi_UVXZ=ent(Piv_UVXZ);
%% Calculate the entropy of (UV)
Hi_UV=ent(reshape(Pi_UV,[],1));

%% This yields an interval for mu
mui=Ri+Ri-Hi_UV-Hi_XZ+Hi_UVXZ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% We now have the point (mui, Ri) as a conterexample                  %
%% We need to upper bound the inner bound and show that there is a gap %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Prepare the grid
Rv1=linspace(0,bdp,floor(K/4))';
Rv2=( 1-linspace((1-bdp)^(1/Q),0,floor(3*K/4)).^Q )';
Rv=unique([Rv1; Rv2]);
K=length(Rv);

%% Intervals:
Rvi=infsupdec(Rv);

%% Probabilities:
avi=infsupdec(zeros(K,1));
for i=1:K
  avi(i)=ib.binentinv_int(1-Rvi(i,1));
end

[aami,bbmi]=meshgrid(avi);
[RRm1, RRm2]=meshgrid(Rv);
[RRmi1, RRmi2]=meshgrid(Rvi);

%% Calculate the intervals for mu
MUmi=1-binent(star(star(aami,ppi),bbmi));

%% Shift everything
MUmi(1:K-1,1:K-1)=MUmi(2:K,2:K);
MUmi(K,1:K-1)=MUmi(K,2:K);
MUmi(1:K-1,K)=MUmi(2:K,K);

%% Construct vectors from arrays
RRv1=[reshape(RRm1, [] ,1)];
RRv2=[reshape(RRm2, [] ,1)];
RRvi1=[reshape(RRmi1, [] ,1)];
RRvi2=[reshape(RRmi2, [] ,1)];
MUvi=[reshape(MUmi, [] ,1)];

%% Use the supremum of the intervals for calculating the upper concave envelope
MUv=sup(MUvi);

%% Get the upper concave envelope, triangulated
T=upper_concave_envelope3_tri(RRv1, RRv2, MUv);

i=tsearchn([RRv1, RRv2], T, sup([Ri, Ri]));
t=T(i,:);
q=hyperplane3(RRvi1(t,:),RRvi2(t,:),MUvi(t,:));
assert(verify_hyperplane(q,RRvi1,RRvi2,MUvi));
assert(strictprecedes(q(3),infsupdec('0')));
MUi=z_interp2(q,Ri,Ri);

printf('Distance is at least %e\n', inf(mui - MUi));
