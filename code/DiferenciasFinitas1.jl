function FDMParallel(nx,nz,nt,t,C,x0,xm,z0,zm,fpk,that,xs,zs)
 dx=abs(xm-x0)/(nx-1)
 dz=abs(zm-z0)/(nz-1)
 dt=t/(nt-1)
 PO=SharedArray(zeros(nx,nz))
 P1=SharedArray(zeros(nx,nz))
 PAUX=SharedArray(zeros(nx,nz))
 jx=floor((xs-x0)/dx)
 jz=floor((zs-z0)/dz)
 for it=2*dt:dt:t
  @parallel  for iz=2:(nz-1)
   for ix=2:(nx-1)
    H=(P1[ix-1,iz]-2*P1[ix,iz]+P1[ix+1,iz])/(dx^2)+(P1[ix,iz-1]-2*P1[ix,iz]+P1[ix,iz+1])/(dz^2)
    if ix==jx && iz==jz
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

CC=1000*ones(200,200)

@time pcolormesh(FDMParallel(200,200,500,.8,CC,-500,500,-500,500,30,1/15,0,0))

SharedArray(zeros(2))

workers()
addprocs(3)

SharedArray{Float64}(1000,10)

@time pmap((nx,nz,nt,t,C,x0,xm,z0,zm,fpk,that,xs,zs)->FDMParallel(nx,nz,nt,t,C,x0,xm,z0,zm,fpk,that,xs,zs),
                                    [(200,200,400,.4,CC,-500,500,-500,500,30,1/15,0,0)])

using PyPlot

a = SharedArray{Float64}(1000)
@parallel for i = 1:1000
    @show a[i] = i
end
a


PO0=zeros(2000,2000)
@time for it=0:10
 for iz=1:2000
   for ix=1:2000
    PO0[ix,iz] =it+iz+ix
  end
 end
end

PO1=SharedArray{Float64}(zeros(2000,2000))
@time for it=0:10
 for iz=1:2000
  @parallel for ix=1:2000
    @show PO1[ix,iz] = it+iz+ix
  end
 end
end

PO2=SharedArray{Float64}(zeros(100,100))
@time for it=0:10
 for iz=1:100
  @parallel for ix=1:100
    PO2[ix,iz] = it+iz+ix
  end
 end
end


PO3=SharedArray{Float64}(10000)
@time @parallel for ix=1:10000
  PO3[ix] = ix
end

PO4=zeros(10000)
@time for ix=1:10000
  PO4[ix] = ix
end

PO5=SharedArray{Float64}(zeros(10000))
@time @parallel for ix=1:10000
  @show PO3[ix] = ix
end
