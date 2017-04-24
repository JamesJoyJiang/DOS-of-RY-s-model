function H0=Hamilton_RY_bcc()
global id N gamma ic

s_x=[0,1;1,0];s_y=[0,-ic;ic,0];s_z=[1,0;0,-1];eye2=eye(2);


dim=2*N^3;
hxp=sparse(dim,dim);
hyp=sparse(dim,dim);
hzp=sparse(dim,dim);
h0=sparse(dim,dim);
H0=sparse(dim,dim);

%bulk hopping

for iz=1:N
    for iy=1:N
        for ix=1:N
            h0(id(:,ix,iy,iz),id(:,ix,iy,iz))=(2+gamma)*s_z;
        end;
    end;
end;

for iz=1:N
    for iy=1:N
        for ix=2:N
            hxp(id(:,ix-1,iy,iz),id(:,ix,iy,iz))=(-i*s_x-s_z)/2;
        end;
    end;
end;

for iz=1:N
    for iy=2,N
        for ix=1:N
            hyp(id(:,ix,iy-1,iz),id(:,ix,iy,iz))=(-i*s_y-s_z)/2;
        end;
    end;
end;

for iz=2:N
    for iy=1:N
        for ix=1:N
            hzp(id(:,ix,iy,iz-1),id(:,ix,iy,iz))=(-s_z)/2;
        end;
    end;
end;
%the boundary terms
for im=1:N
    for in=1:N
      hxp(id(:,N,im,in),id(:,1,im,in))=(-i*s_x-s_z)/2;
      hyp(id(:,im,1,in),id(:,im,1,in))=  (-i*s_y-s_z)/2;
      hzp(id(:,im,in,1),id(:,im,in,1))=  (-s_z)/2;
    end;
end


H0=hxp+hyp+hzp;H0=H0+H0';
H0=H0+h0;


