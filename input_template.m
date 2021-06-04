clear all;
close all;
%% dims
input.dims.nt=1000;
input.dims.nx=2000;
input.dims.nz=1000;
input.dims.X=[0,20000]; % min and max true coordinate in x-axis
input.dims.Z=[0,10000]; % min and max true coordinate in z-axis
input.dims.dt=10^-3;
input.dims.dx=10;
input.dims.dz=10;
%% PML
input.PML.lp=20; % PML layers
input.PML.nPML=2; % PML power
input.PML.Rc=.001; % theoretical reflection coefficient
%% material parameter
input.material_parameter.C44=ones(input.dims.nx,input.dims.nz);
input.material_parameter.C46=ones(input.dims.nx,input.dims.nz);
input.material_parameter.C66=ones(input.dims.nx,input.dims.nz);
input.material_parameter.ts2=ones(input.dims.nx,input.dims.nz); %
input.material_parameter.ts4=ones(input.dims.nx,input.dims.nz); %
input.material_parameter.phi2=ones(input.dims.nx,input.dims.nz); %
input.material_parameter.phi4=ones(input.dims.nx,input.dims.nz); %
input.material_parameter.rho=ones(input.dims.nx,input.dims.nz);
input.material_parameter.C44(:)=10^9;
input.material_parameter.C66(:)=10^9;
input.material_parameter.ts2(:)=.1;
input.material_parameter.ts4(:)=.1;
input.material_parameter.phi2(:)=-.01;
input.material_parameter.phi4(:)=-.01;
input.material_parameter.rho(:)=10^3;
%% source
input.source.s1=30; % source grid location 1
input.source.s3=30; % source grid location 3
input.source.s1t=300; % source true location 1
input.source.s3t=300; % source true location 3
freq=5;
M=2;
ricker=wavelet(freq,input.dims.dt,input.dims.nt,M);
input.source.src=ricker; % source signals
%% receiver
input.receiver.r1=30; % receiver grid location 1
input.receiver.r3=30; % receiver grid location 3
input.receiver.r1t=300; % receiver true location 1
input.receiver.r3t=300; % receiver true location 3
%% visualization
input.visualization.plot_interval=100; % plot inverval
p2=mfilename('fullpath');
if ~exist(p2,'dir')
    mkdir(p2)
end
input.visualization.path=p2; % location for output folder
%% pass input to solver
[u2,R]=SH_2D_solver(input);