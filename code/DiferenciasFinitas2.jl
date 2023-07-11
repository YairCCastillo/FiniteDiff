function finite(nx,nz,nt,t,C,x0,xm,z0,zm,fpk,that,xs,zs)
 dx=abs(xm-x0)/(nx-1)
 dz=abs(zm-z0)/(nz-1)
 dt=t/(nt-1)
 PO=zeros(nx,nz)
 P1=zeros(nx,nz)
 jx=floor((xs-x0)/dx)
 jz=floor((zs-z0)/dz)
 for it=2*dt:dt:t
  for iz=2:(nz-1)
   for ix=2:(nx-1)
    H=(P1[ix-1,iz]-2*P1[ix,iz]+P1[ix+1,iz])/(dx^2)+(P1[ix,iz-1]-2*P1[ix,iz]+P1[ix,iz+1])/(dz^2)
    if ix==jx && iz==jz
     PO[ix,iz]=2*P1[ix,iz]-PO[ix,iz]+(C[ix,iz]^2)*(dt^2)*H+exp(-(fpk^2)*(pi^2)*(it-that)^2)*(dt^2)
    else
     PO[ix,iz]=2*P1[ix,iz]-PO[ix,iz]+(C[ix,iz]^2)*(dt^2)*H
    end
   end
  end
  PAUX=deepcopy(P1)
  P1=deepcopy(PO)
  PO=deepcopy(PAUX)
 end
 return PO
end

C=1000*ones(500,500)

using PyPlot

pcolormesh(finite(500,500,800,0.8,C,-500,500,-500,500,30,1/15,0,0))

finite(200,200,5000,.8,C,-500,500,-500,500,30,1/15,0,0)
 pcolormesh(as)
  colorbar()

function finite2(nx,nz,nt,t,C,x0,xm,z0,zm,fpk,that,xs,zs,xs2,zs2)
 dx=abs(xm-x0)/(nx-1)
 dz=abs(zm-z0)/(nz-1)
 dt=t/(nt-1)
 PO=zeros(nx,nz)
 P1=zeros(nx,nz)
 jx=floor((xs-x0)/dx)
 jz=floor((zs-z0)/dz)
 jx2=floor((xs2-x0)/dx)
 jz2=floor((zs2-z0)/dz)
 for it=2*dt:dt:t
  for iz=2:(nz-1)
   for ix=2:(nx-1)
    H=(P1[ix-1,iz]-2*P1[ix,iz]+P1[ix+1,iz])/(dx^2)+(P1[ix,iz-1]-2*P1[ix,iz]+P1[ix,iz+1])/(dz^2)
    if ix==jx && iz==jz
     PO[ix,iz]=2*P1[ix,iz]-PO[ix,iz]+(C[ix,iz]^2)*(dt^2)*H+exp(-(fpk^2)*(pi^2)*(it-that)^2)*(dt^2)
    elseif ix==jx2 && iz==jz2
     PO[ix,iz]=2*P1[ix,iz]-PO[ix,iz]+(C[ix,iz]^2)*(dt^2)*H+exp(-(fpk^2)*(pi^2)*(it-that)^2)*(dt^2)
    else
     PO[ix,iz]=2*P1[ix,iz]-PO[ix,iz]+(C[ix,iz]^2)*(dt^2)*H
    end
   end
  end
  PAUX=P1
  P1=PO
  PO=PAUX
 end
 return PO
end

using PyPlot

as=finite2(400,400,5000,.8,1000*ones(400,400),-500,500,-500,500,30,1/15,0,0,50,50)
 pcolormesh(as)
  colorbar()
