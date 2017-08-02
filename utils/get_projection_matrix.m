%% Construct a symbolic projection matrix
function Nproj=get_projection_matrix()
  %% Construct Matrix A for checking Conditions
  A=[];

  %% Sum to one
  P_UVgXZ = zeros(2,2,2,2);
  P_UVgXZ(:,:,1,1)=1;
  A=[A;reshape(P_UVgXZ,1,[])];

  P_UVgXZ = zeros(2,2,2,2);
  P_UVgXZ(:,:,1,2)=1;
  A=[A;reshape(P_UVgXZ,1,[])];

  P_UVgXZ = zeros(2,2,2,2);
  P_UVgXZ(:,:,2,1)=1;
  A=[A;reshape(P_UVgXZ,1,[])];

  P_UVgXZ = zeros(2,2,2,2);
  P_UVgXZ(:,:,2,2)=1;
  A=[A;reshape(P_UVgXZ,1,[])];

  %% 1st short Markov chain [ U -- X -- Z ] 
  P_UVgXZ = zeros(2,2,2,2);
  P_UVgXZ(1,:,1,1)=1;
  P_UVgXZ(1,:,1,2)=-1;
  A=[A;reshape(P_UVgXZ,1,[])];

  P_UVgXZ = zeros(2,2,2,2);
  P_UVgXZ(1,:,2,1)=1;
  P_UVgXZ(1,:,2,2)=-1;
  A=[A;reshape(P_UVgXZ,1,[])];

  %% 2nd short Markov chain [ X -- Z -- V ] 
  P_UVgXZ = zeros(2,2,2,2);
  P_UVgXZ(:,1,1,1)=1;
  P_UVgXZ(:,1,2,1)=-1;
  A=[A;reshape(P_UVgXZ,1,[])];

  P_UVgXZ = zeros(2,2,2,2);
  P_UVgXZ(:,1,1,2)=1;
  P_UVgXZ(:,1,2,2)=-1;
  A=[A;reshape(P_UVgXZ,1,[])];

  %% Preserve marginal p(u|x)
  P_UVgXZ = zeros(2,2,2,2);
  P_UVgXZ(1,:,1,1)=1;
  A=[A;reshape(P_UVgXZ,1,[])];

  P_UVgXZ = zeros(2,2,2,2);
  P_UVgXZ(1,:,2,1)=1;
  A=[A;reshape(P_UVgXZ,1,[])];

  %% Preserve marginal p(v|z)
  P_UVgXZ = zeros(2,2,2,2);
  P_UVgXZ(:,1,1,1)=1;
  A=[A;reshape(P_UVgXZ,1,[])];

  P_UVgXZ = zeros(2,2,2,2);
  P_UVgXZ(:,1,1,2)=1;
  A=[A;reshape(P_UVgXZ,1,[])];

  %% Obtain a symbolic projection onto the null space of A
  A_sym=sym(A);
  N_sym=null(A_sym);
  North=orth(N_sym);
  Nproj=North*North';
end
