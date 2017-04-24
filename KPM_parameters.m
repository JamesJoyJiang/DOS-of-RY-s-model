
clear;
global id N gamma ic
sq3=sqrt(3d0);ic=sqrt(-1);
N=30; % the lattice sze for each dimension(Nx=Ny=Nz)

gamma=0.5;



M=2000;  %the expansion terms
Mb=2*M;  %number of the data points
N_disorder=1;  %number of disorder configurations :self-avareage is important!

R=100;      % random integration
eplison=1d-3;
W=0; %the strength of disorder

i=0;
Jac0=0;Jac=0;
Jac0=1/(M+1)*((M-i+1)*cos(pi*i/(M+1))+sin(pi*i/(M+1))*cot(pi/(M+1)));
for i=1:M-1
Jac(i)=1/(M+1)*((M-i+1)*cos(pi*i/(M+1))+sin(pi*i/(M+1))*cot(pi/(M+1)));
end;

id=0;
for nz=1:N
    for ny=1:N
        for nx=1:N
            for a=1:2
              id(a,nx,ny,nz)=a+(nx-1)*2+(ny-1)*N*2+(nz-1)*N^2*2;
            end
        end;
    end;
end;


    
