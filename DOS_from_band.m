
N=30;NE=2000;
rho=0;rho(NE)=0;

s_x=[0,1;1,0];s_y=[0,-ic;ic,0];s_z=[1,0;0,-1];eye2=eye(2);
eta=10^-4;

for kx=-pi:2*pi/N:pi-2*pi/N
    for ky=-pi:2*pi/N:pi-2*pi/N
        for kz=-pi:2*pi/N:pi-2*pi/N
          HH=(2+gamma-cos(kx)-cos(ky)-cos(kz))*s_z+sin(kx)*s_x+sin(ky)*s_y;
          r=eig(HH);
          
          for s=1:2
              nE=0;
            for EE=-1:2/NE:1-2/NE;
              nE=nE+1;
              rho(nE)=rho(nE)+eta*((EE-r(s))^2+eta^2)^-1;
            end;
          end;
        end;
    end;
end
figure;plot(-1:2/NE:1-2/NE,rho)