function [u2,R]=SH_2D_solver(input)
%% initialize variables in time stepping
u2=zeros(input.dims.nx,input.dims.nz,3);

sigma23=zeros(input.dims.nx,input.dims.nz);
sigma12=zeros(input.dims.nx,input.dims.nz);
e23l=zeros(input.dims.nx,input.dims.nz);
e12l=zeros(input.dims.nx,input.dims.nz);
u3=zeros(input.dims.nx,input.dims.nz);
acceleration=zeros(input.dims.nx,input.dims.nz);

f1=zeros(input.dims.nx,input.dims.nz);
f2=zeros(input.dims.nx,input.dims.nz);

e23l=zeros(input.dims.nx,input.dims.nz);
e12l=zeros(input.dims.nx,input.dims.nz);
e23l_last=zeros(input.dims.nx,input.dims.nz);
e12l_last=zeros(input.dims.nx,input.dims.nz);

dsigma23=zeros(input.dims.nx,input.dims.nz);
dsigma12=zeros(input.dims.nx,input.dims.nz);

sigma23=dsigma23;
sigma12=dsigma12;
e23=zeros(input.dims.nx,input.dims.nz);
e12=zeros(input.dims.nx,input.dims.nz);
e23_last=zeros(input.dims.nx,input.dims.nz);
e12_last=zeros(input.dims.nx,input.dims.nz);

int_dsigma23=zeros(input.dims.nx,input.dims.nz);
int_dsigma12=zeros(input.dims.nx,input.dims.nz);

int_sigma23=zeros(input.dims.nx,input.dims.nz);
int_sigma12=zeros(input.dims.nx,input.dims.nz);


R=zeros(input.dims.nt,length(input.receiver.r3));
%% receiver ind
ind_r=sub2ind(size(u2),input.receiver.r1,input.receiver.r3,3);
ind_s=sub2ind(size(u2),input.source.s1,input.source.s3,3);
%% PML
beta0=sqrt(input.material_parameter.C44./input.material_parameter.rho)*(input.PML.nPML+1) ...
    *log(1/input.PML.Rc)/2/input.PML.lp/(.5*input.dims.dx+.5*input.dims.dz);
beta1=zeros(input.dims.nx,input.dims.nz);
beta3=beta1;
tt=(1:input.PML.lp)/input.PML.lp;
tt2=repmat(reshape(tt,[input.PML.lp,1]),[1,input.dims.nz]);
plane_grad1=zeros(input.dims.nx,input.dims.nz);
plane_grad3=plane_grad1;

plane_grad1(2:input.PML.lp+1,:)=flip(tt2,1);
plane_grad1(input.dims.nx-input.PML.lp:end-1,:)=tt2;
plane_grad1(1,:)=plane_grad1(2,:);
plane_grad1(end,:)=plane_grad1(end-1,:);

tt2=repmat(reshape(tt,[1,input.PML.lp]),[input.dims.nx,1]);
plane_grad3(:,2:input.PML.lp+1)=flip(tt2,2);
plane_grad3(:,input.dims.nz-input.PML.lp:end-1)=tt2;
plane_grad3(:,1)=plane_grad3(:,2);
plane_grad3(:,end)=plane_grad3(:,end-1);

beta1=beta0.*plane_grad1.^input.PML.nPML;
beta3=beta0.*plane_grad3.^input.PML.nPML;
%% weight
weightx1=9/(8*input.dims.dx);
weightx2=-1/(24*input.dims.dx);
weighty1=9/(8*input.dims.dz);
weighty2=-1/(24*input.dims.dz);
kx=3:input.dims.nx-2;
kz=3:input.dims.nz-2;
n_pic=1;
%% create folder
if ~exist([input.visualization.path '/pic/'],'dir')
    mkdir([input.visualization.path '/pic/']);
end
%% time stepping
tic;
for l=1:input.dims.nt-1
    %% shift u2 in time
    for l2=1:2
        u2(:,:,l2)=u2(:,:,l2+1);
    end
    %% Strains
    e23_last=e23;
    e12_last=e12;
    
    e23(kx,kz) = weighty1*(u2(kx,kz,2)-u2(kx,kz-1,2))+ weighty2*(u2(kx,kz+1,2)-u2(kx,kz-2,2));
    e12(kx,kz) = weightx1*(u2(kx,kz,2)-u2(kx-1,kz,2))+ weightx2*(u2(kx+1,kz,2)-u2(kx-2,kz,2));
    %% Memory variables, absorption:
    f1(kx,kz)=2*input.material_parameter.ts2(kx,kz)-input.dims.dt;
    f2(kx,kz)=2*input.material_parameter.ts2(kx,kz)+input.dims.dt;
    e23l_last=e23l;
    e23l(kx,kz)=.5*((2*input.dims.dt*input.material_parameter.ts2(kx,kz).*input.material_parameter.phi2(kx,kz).*e23(kx,kz) ...
        +f1(kx,kz).*e23l(kx,kz))./f2(kx,kz)+e23l_last(kx,kz));
    
    f1(kx,kz)=2*input.material_parameter.ts4(kx,kz)-input.dims.dt;
    f2(kx,kz)=2*input.material_parameter.ts4(kx,kz)+input.dims.dt;
    
    e12l_last=e12l;
    e12l(kx,kz)=.5*((2*input.dims.dt*input.material_parameter.ts4(kx,kz).*input.material_parameter.phi4(kx,kz).*e12(kx,kz) ...
        +f1(kx,kz).*e12l(kx,kz))./f2(kx,kz)+e12l_last(kx,kz));
    %% Stresses:
    sigma23(kx,kz)=input.dims.dt*(input.material_parameter.C44(kx,kz).*(e23(kx,kz)-e23_last(kx,kz))/input.dims.dt ...
        +input.material_parameter.C46(kx,kz).*(e12(kx,kz)-e12_last(kx,kz))/input.dims.dt ...
        +input.material_parameter.C44(kx,kz).*(e23l(kx,kz)-e23l_last(kx,kz))/input.dims.dt ...
        +input.material_parameter.C44(kx,kz).*beta1(kx,kz).*e23(kx,kz) ...
        +input.material_parameter.C46(kx,kz).*beta3(kx,kz).*e12(kx,kz) ...
        -(beta1(kx,kz)+beta3(kx,kz)).*sigma23(kx,kz) ...
        -beta1(kx,kz).*beta3(kx,kz).*int_sigma23(kx,kz)) ...
        +sigma23(kx,kz);
    sigma12(kx,kz)=input.dims.dt*(input.material_parameter.C46(kx,kz).*(e23(kx,kz)-e23_last(kx,kz))/input.dims.dt ...
        +input.material_parameter.C66(kx,kz).*(e12(kx,kz)-e12_last(kx,kz))/input.dims.dt...
        +input.material_parameter.C66(kx,kz).*(e12l(kx,kz)-e12l_last(kx,kz))/input.dims.dt ...
        +input.material_parameter.C46(kx,kz).*beta1(kx,kz).*e23(kx,kz) ...
        +input.material_parameter.C66(kx,kz).*beta3(kx,kz).*e12(kx,kz) ...
        -(beta1(kx,kz)+beta3(kx,kz)).*sigma12(kx,kz) ...
        -beta1(kx,kz).*beta3(kx,kz).*int_sigma12(kx,kz)) ...
        +sigma12(kx,kz);
    
    int_sigma23=int_sigma23+sigma23*input.dims.dt;
    int_sigma12=int_sigma12+sigma12*input.dims.dt;
    %% Differential stresses:
    dsigma23(kx,kz) = weighty1*(sigma23(kx,kz+1)-sigma23(kx,kz))+weighty2*(sigma23(kx,kz+2)-sigma23(kx,kz-1));
    dsigma12(kx,kz) = weightx1*(sigma12(kx+1,kz)-sigma12(kx,kz))+weightx2*(sigma12(kx+2,kz)-sigma12(kx-1,kz));
    
    int_dsigma23=int_dsigma23+dsigma23*input.dims.dt;
    int_dsigma12=int_dsigma12+dsigma12*input.dims.dt;
    %% Acceleration
    acceleration(kx,kz)=((dsigma23(kx,kz)+dsigma12(kx,kz))+beta3(kx,kz).*int_dsigma12(kx,kz)+beta1(kx,kz).*int_dsigma23(kx,kz))./input.material_parameter.rho(kx,kz);
    %% Euler equation
    u2(kx,kz,3)=2*u2(kx,kz,2)-u2(kx,kz,1)+input.dims.dt^2*acceleration(kx,kz) ...
        -input.dims.dt^2*beta1(kx,kz).*beta3(kx,kz).*u2(kx,kz,2) ...
        -input.dims.dt*(beta1(kx,kz)+beta3(kx,kz)).*(u2(kx,kz,2)-u2(kx,kz,1));
    %% source term
    u2(ind_s)=u2(ind_s)+1./input.material_parameter.rho(ind_s-input.dims.nx*input.dims.nz*2).*input.source.src(l,:);
    %%
    for lb=[2,1]
        u2(lb,:,3)=u2(lb+1,:,3);
        u2(:,lb,3)=u2(:,lb+1,3);
    end
    %% Seismograms
    R(l+1,:)=u2(ind_r);
    %%
    % Plot
    if mod(l,input.visualization.plot_interval)==0
        figure('visible','off');
        imagesc(input.dims.X,input.dims.Z,u2(:,:,3)');
        xlabel('x');
        ylabel('y');
        colorbar;
        print([input.visualization.path '/pic/' num2str(n_pic)],'-dpng','-r100');
        n_pic=n_pic+1;
    end
    fprintf('\n%d',l);
end
toc;
end