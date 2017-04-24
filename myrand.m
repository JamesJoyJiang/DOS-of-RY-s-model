
KPM_parameters;
H0=Hamilton_RY_bcc();

filename=strcat('Gamma= ',num2str(gamma),',Disorder= ',...
    num2str(W), ',N=',num2str(N));
filename
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  diagonal terms   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fff(N_disorder,Mb)=0;
A_disorder(N_disorder)=0;
B_disorder(N_disorder)=0;

opts.disp=0; 
dim=2*N^3;

parfor n_disorder=1:N_disorder
ss=W*(rand(dim,1)-0.5);
h0=sparse(reshape(id,dim,1),reshape(id,dim,1),ss,dim,dim);  %on-site disorder


n_disorder

H=h0+H0;

  v=eigs(H,2, 'LM',opts);
  v=sort(v);
  if (prod(v)>0)
     return;
  end
a=(v(2)-v(1))/(2-eplison);
b=(v(2)+v(1))/2;
[a,b];
A_disorder(n_disorder)=v(1);
B_disorder(n_disorder)=v(2);
issparse(H);
H=(H-b*speye(dim))/a;
%issparse(H);
%toc;
%tic;
%"now expansion";
%v;return

mu=0;mu(M-1)=0;
mu0=0;

for ir=1:R
    ir
    v0=0;
    for i=1:dim
      v0(i)=((rand -0.5)*sq3*2);
    end;
   v0=sparse(v0');
   v00=v0;
   
   mu0=v0'*v0+mu0;
   v1=H*v0;

   mu(1)=v1'*v0+mu(1);
  
   v2=2*H*v1-v0;
   mu(2)=v2'*v0+mu(2);;
  
   for im=3:M-1;
     v0=v1;v1=v2;
      v2=2*H*v1-v0;
      mu(im)=v2'*v00+mu(im);
  end;
end
mu0=mu0/R; mu=mu/R;

mu0=mu0*Jac0;
mu=mu.*Jac;

x0=cos(pi/Mb*(1d0/2));
x=cos(pi/Mb*(1d0/2+[1:Mb-1]));  arcx=acos(x);

gma0=mu0+2*mu*cos([1:M-1]'*pi/Mb/2);

gma=mu0+2*mu*cos([1:M-1]'*arcx);

f0=gma0/sqrt(1-x0^2)/pi;
f=gma./sqrt(1-x.^2)/pi;
f=[f0,f]/N^3;
x=[x0,x];

x=x*a+b;  %normaliztion: 2017
fff(n_disorder,:)=f/a; %normaliztion: 2017
plot(x,sum(fff,1)/n_disorder);

figure; plot(x,sum(fff,1)/n_disorder); ylabel('DOS');xlabel('E');
title(filename);
savefig(strcat(filename,'DOS','.fig'));

S = struct('E_region',E_region,'DOS',DOS,'fff',real(fff));
filenamemat=strcat(filename,'.mat')
save(filenamemat,'S');

end



%s=0;%for i=1:2000;for j=1:2000;s=s+cos((i+j)*pi*0.01);end;end

%for i=1:M
%Ind(0)
%mu0=0;mu(1:M-1)=rand(1,M-1);
%ramda(1)=mu0;
%ramda(2:M)=2*mu.*exp(ic*pi*[2:M]/(2*Mb));
%ramda(M+1:Mb)=0;
%ramdaBar=fft(ramda,Mb);
%x=cos(pi/Mb*(1d0/2+[0:Mb-1]));
%stop




%s=0;N=100;

%for i=1:N
%s=s+((rand -0.5)*sq3*2);
%end
%s=s/N;
