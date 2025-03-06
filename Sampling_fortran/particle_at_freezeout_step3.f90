!******************************************************************************************************
!! Example code of various particle sampling from a freezeout surface!!!
!! gfortran particle_at_freezeout_step3.f90 !!
!! prerequisites pyfort.py (calculates total number of hadrons from each cell) and my_file_mass.txt (various hadron informations) !!
!******************************************************************************************************

real:: aa,bb,cc,dd,ee,ff,gg,hh,ii,jj,kk,ll,mm,nn,oo,pp,qq,rr,ss,tt,uu,vv,ww,xx,yy,zz
real:: efakt,skapro,ska,expon,kurz,mt,max_spec,max_val,y1,y2,u,u1,u2,pre_fac,m,y,y_local,etas
real,parameter:: hbarc=0.19733,y_minus_eta_cut = 4.0
real:: phi,pt,pi=dacos(-1.0d0),min_pt=0.01,max_pt=3.0,ymax=5.0d0,prob,prob2,val,log_val,deg 
integer:: n,sign1,iii,jjj,tot_had,flag,flag_rand=1,flag_pos,lll,mmm,i,j,k,l,pti,phii,sign1_B,ik,ik2
integer, parameter :: npt=16,nphi=40,ny=51
real,dimension(ny+1,npt,nphi)::dndyptdptdphi
real,dimension(npt,nphi)::temp_cell
real,allocatable,dimension(:)::numb,x_f,y_f,tau_f,eta_f,ef_Plus_Pf_by_Tf, pi_tautau, pi_taux, pi_tauy, pi_tauetat
real,allocatable,dimension(:)::u_tau, u_x, u_y, u_etat, ef,Tf, mu_B,pi_yetat,pi_etatetat, pi_xetat, pi_yy, pi_xx, pi_xy
real,allocatable,dimension(:):: d3sigma_tau, d3sigma_x, d3sigma_y, d3sigma_eta
character(len=15)::chr
real,dimension(nphi)::cos_phi,sin_phi
real,dimension(npt)::pt_array
real::cosh_y_local,sinh_y_local,cosh_eta_s,sinh_eta_s,deltay,y3,prap,rapidity,cosh_y,sinh_y

  call system("rm my_file.txt")
  call system("rm spec_had_phipt.dat")
  
!**********************************************************************!
!********************** total number of freezeout elements ********************!
!***********************************************************************
call countlines(n)
write(*,*)n
allocate(numb(n))
allocate(x_f(n))
allocate(y_f(n))
allocate(tau_f(n))
allocate(eta_f(n))
allocate(d3sigma_tau(n))
allocate(d3sigma_x(n))
allocate(d3sigma_y(n))
allocate(d3sigma_eta(n))
allocate(u_tau(n))
allocate(u_x(n))
allocate(u_y(n))
allocate(u_etat(n))
allocate(ef(n))
allocate(Tf(n))
allocate(mu_B(n))
allocate(ef_Plus_Pf_by_Tf(n))
allocate(pi_tautau(n))
allocate(pi_taux(n))
allocate(pi_tauy(n))
allocate(pi_tauetat(n))
allocate(pi_xx(n))
allocate(pi_xy(n))
allocate(pi_xetat(n))
allocate(pi_yy(n))
allocate(pi_yetat(n))
allocate(pi_etatetat(n))
!**********************************************************************!
!*******************   starts to read surface information  ***********************!
!***********************************************************************
open(10,file="surface_exmp3.dat") !!! surface file
open(210,file="spec_had_phipt.dat") !!! diff  spectra stored this file for each hdrons
 iii=1
do while (iii.le.n) 
read(10,*)aa,bb,cc,dd,ee,ff,gg,hh,ii,jj,kk,ll,mm,nn,oo,pp,qq,rr,ss,tt,uu,vv,ww,xx,yy,zz
tau_f(iii)=aa !! spacetime informations
x_f(iii)=bb   !! needs to be recalled later
y_f(iii)=cc   !! needs to be recalled later
eta_f(iii)=dd   !! needs to be recalled later
d3sigma_tau(iii)=ee
d3sigma_x(iii)=ff
d3sigma_y(iii)=gg
d3sigma_eta(iii)=hh
u_tau(iii)=ii
u_x(iii)=jj
u_y(iii)=kk
u_etat(iii)=ll
ef(iii)=mm
Tf(iii)=nn*hbarc
mu_B(iii)=oo
ef_Plus_Pf_by_Tf(iii)=pp
pi_tautau(iii)=qq
pi_taux(iii)=rr
pi_tauy(iii)=ss
pi_tauetat(iii)=tt
pi_xx(iii)=uu
pi_xy(iii)=vv
pi_xetat(iii)=ww
pi_yy(iii)=xx
pi_yetat(iii)=yy
pi_etatetat(iii)=zz
iii=iii+1
end do
 close(10)
!**********************************************************************!
!*******************  openning mass information file    ************************!
!***********************************************************************
open(30,file="my_file_mass.txt")
!**********************************************************************!
open(18,file="sampling_final.txt")
do imass= 1,1  !!! hadron loop
read(30,*)m,sign1_B,spin,chr
!***************check whether baryon or meson************************!
!**********************************************************************!
if (sign1_B.eq.0)then
sign1=-1.0d0
else
sign1=1.0d0
endif
!***************check degeneracy************************!
!**********************************************************************!
deg=2*spin+1.0d0
pre_fac=deg/(((2.*pi)** 3.)*(hbarc** 3.))
!**********************************************************************!
!***********************************************************************!
!****************calculating differential spectra *******************************!
!***********************************************************************!     
call clean_cell(dndyptdptdphi)
call clean_cell1(numb,n)

   do k=1,nphi
        phi=0.0+(k-1)*(2*pi/nphi)   
        cos_phi(k)=cos(phi)
        sin_phi(k)=sin(phi)
   end do     
   
  do j=1,npt
     pt_array(j)=min_pt+(max_pt-min_pt)*(float((j-1)**2)/float((npt-1)**2))
  end do
  
do i=1,ny+1
     deltay=2*abs(ymax)/(ny)
     etay=-ymax+deltay*(i-1)   
       y_local=etay           
       rapidity = y_local
       cosh_y = cosh(y_local)
       sinh_y = sinh(y_local)  
           
  do  iii=1,n
     call clean_cell2(temp_cell)
     etas=eta_f(iii)
     cosh_eta_s=cosh(etas)
     sinh_eta_s=sinh(etas)
     do j=1,npt

      pt = pt_array(j)
      cosh_y_local = cosh_y
      sinh_y_local = sinh_y
      mt = sqrt(m*m + pt*pt)    ! all in GeV
      y = rapidity
          if (abs(y - etas) .lt. y_minus_eta_cut) then
           ptau = mt*(cosh_y_local*cosh_eta_s- sinh_y_local*sinh_eta_s)
           peta = mt/tau_f(iii)*(sinh_y_local*cosh_eta_s- cosh_y_local*sinh_eta_s)
              do k=1,nphi
            px = pt*cos_phi(k)
            py = pt*sin_phi(k)
          
	  psigma_dsigma=tau_f(iii)*(ptau*d3sigma_tau(iii)+px*d3sigma_x(iii)+py*d3sigma_y(iii)+peta*d3sigma_eta(iii))
	  E=(ptau*u_tau(iii)-px*u_x(iii)-py*u_y(iii)-tau_f(iii)*peta*u_etat(iii))
	  
	  w_factor = (ptau*pi_tautau(iii)*ptau-2*ptau*pi_taux(iii)*px-2*ptau*pi_tauy(iii)*py-2*tau_f(iii)*tau_f(iii)* &
		ptau*pi_tauetat(iii)/tau_f(iii)*peta+px*pi_xx(iii)*px+2*px*pi_xy(iii)*py+2*tau_f(iii)*tau_f(iii)* &
		px*pi_xetat(iii)/tau_f(iii)*peta+py*pi_yy(iii)*py+2*tau_f(iii)*tau_f(iii)*py*pi_yetat(iii)/tau_f(iii)*peta+ &
		tau_f(iii)*tau_f(iii)*tau_f(iii)*tau_f(iii)*peta*pi_etatetat(iii)/tau_f(iii)/tau_f(iii)*peta)
	  pre_factor_shear=1./(2.*ef_Plus_Pf_by_Tf(iii)*Tf(iii)*Tf(iii)*Tf(iii))*hbarc

	  f=1./(exp(1./Tf(iii)*(E - mu_B(iii))) + sign1)

	  delta_f=(f*(1.-sign1*f)*pre_factor_shear*w_factor)
	  spcell=(f+delta_f)*psigma_dsigma
         
          dndyptdptdphi(i,j,k)=dndyptdptdphi(i,j,k)+spcell
          temp_cell(j,k)=spcell
       end do   
      end if  
    end do

!***********************************************************************
!********************* diff spectra ready ***********************************!
!***********************************************************************
!***********************************************************************
!****************** max log value for each cell  stored ************************!
!***********************************************************************
max_spec=numb(iii)
do j=1,npt
  do k=1,nphi
    if (temp_cell(j,k).le. 0.0d0)then
    temp_cell(j,k)=10**(-30)
    end if
    xx=temp_cell(j,k)*pre_fac
    if(max_spec.lt.log(xx)) then
       max_spec=log(xx)
    end if
  end do
end do  
numb(iii)=max_spec
!**************************************************************************
!**************************************************************************
   end do

end do
!**************************************************************************
!**************************************************************************
max_val=-10.0**34
do  i=1,ny+1
   do j=1,npt
       do k=1,nphi
       if (dndyptdptdphi(i,j,k) .lt. 0.0d0)then
       dndyptdptdphi(i,j,k)=0.0d0
       end if   
     log_val=dndyptdptdphi(i,j,k)*pre_fac
     if(max_val.lt.log(log_val))then
        max_val=log(log_val)
     end if   
       write(210,*)dndyptdptdphi(i,j,k)*pre_fac,(k-1)*(2*pi/nphi),pt_array(j),-ymax+deltay*(i-1)
       end do
     end do
end do                     
  close(210)
!**************************************************************************
!**************************************************************************  
!do iii=1,n
!write(*,*)numb(iii)
!end do
 !***********************************************************************
 !***************** calculates tot number of a hadron  species***********************!
 !***********************************************************************
  call system("python3 pyfort.py")
 open(20,file="my_file.txt")
 read(20,*)tot_had 
 write(*,*)tot_had
!**********************************************************************!
!***********************************************************************
!*******************      particle sampling intiates            ***********************!
!***********************************************************************
i=0
call change_seed()
do   while(i.lt.tot_had)
flag=0
       do while (flag.eq.0)
				call random_number(y1)
				call random_number(y2)
				call random_number(y3)
				call random_number(u)
				
				phi=y2*(2.0*pi) !!! phi sampled
				prap=-ymax+2*y3*ymax   !!!  scenario y can be any value 
                                y_local=prap

				pt= sqrt(y1*max_pt*max_pt) !!! pt sampled
				
				if (pt.lt.min_pt)then
				exit
				endif
	      call giveprob(phi,pt,y_local,dndyptdptdphi,prob,max_val)  !!! probability of a phase-space value : phi,pt,rap
	      if (u.lt.prob)then  !!! prob checked acceptance rejection method
				i=i+1
              flag_pos=0
!          write(18,*)y_local,pt*cos(phi),pt*sin(phi),i,m,chr !!! writes particle informations x,y,eta,tau,px,py,number,mass,name				


!**************************************************************************
!************ sampling surface cell which can emit the above particle  ***********!
!**************************************************************************

     do while (flag_pos.eq.0)
     phi=y2*(2.0*pi) !!! phi sampled
     prap=-ymax+2*y3*ymax   !!!  scenario y can be any value 
     y_local=prap
     pt= sqrt(y1*max_pt*max_pt) !!! pt sampled	  
	call random_number(u1)
	call random_number(u2)
	      
	j = 1 + floor(n*u1)  !!! We want to choose one from n surface cells
	!write(*,*)j
	mt = sqrt(m*m + pt*pt)    ! all in GeV
	     etas=eta_f(j)
             cosh_eta_s=cosh(etas)
             sinh_eta_s=sinh(etas)
             cosh_y_local = cosh(y_local)
             sinh_y_local = sinh(y_local)
             if (abs(y_local - etas) .lt. y_minus_eta_cut) then
             ptau = mt*(cosh_y_local*cosh_eta_s- sinh_y_local*sinh_eta_s)
             peta = mt/tau_f(j)*(sinh_y_local*cosh_eta_s- cosh_y_local*sinh_eta_s)
             px=pt*cos(phi)
             py=pt*sin(phi)
             psigma_dsigma=tau_f(j)*(ptau*d3sigma_tau(j)+px*d3sigma_x(j)+py*d3sigma_y(j)+peta*d3sigma_eta(j))
			 E=(ptau*u_tau(j)-px*u_x(j)-py*u_y(j)-tau_f(j)*peta*u_etat(j))
				
				w_factor = (ptau*pi_tautau(j)*ptau-2*ptau*pi_taux(j)*px-2*ptau*pi_tauy(j)*py-2*tau_f(j)*tau_f(j)* &
				ptau*pi_tauetat(j)/tau_f(j)*peta+px*pi_xx(j)*px+2*px*pi_xy(j)*py+2*tau_f(j)*tau_f(j)* &
				px*pi_xetat(j)/tau_f(j)*peta+py*pi_yy(j)*py+2*tau_f(j)*tau_f(j)*py*pi_yetat(j)/tau_f(j)*peta+ &
				tau_f(j)*tau_f(j)*tau_f(j)*tau_f(j)*peta*pi_etatetat(j)/tau_f(j)/tau_f(j)*peta)
				pre_factor_shear=1./(2.*ef_Plus_Pf_by_Tf(j)*Tf(j)*Tf(j)*Tf(j))*hbarc

				f=1./(exp(1./Tf(j)*(E - mu_B(j))) + sign1)

				delta_f=(f*(1.-sign1*f)*pre_factor_shear*w_factor)
				spcell=(f+delta_f)*psigma_dsigma
				prob2=exp(log(spcell*pre_fac)-numb(j))
				else
				prob2=0.00
				end if
	if (u2.lt.prob2) then						
	write(18,*)x_f(j),y_f(j),eta_f(j),tau_f(j),y_local,pt*cos(phi),pt*sin(phi),i,m,chr !!! writes particle informations x,y,etas,tau,rap,px,py,number,mass,name
	flag_pos=1
	end if
				   end do
         !***************************************************************************!	
				 	flag=1
	       else
	       flag=0
	       end if

            end do
end do
!***************************************************************************!
!****************************** deallocating arrays******************************!
!***************************************************************************!
 close(20)
 !*************************************************************************!
 !***************************** removing files ********************************!
!*************************************************************************!
  call system("rm my_file.txt")
  call system("rm spec_had_phipt.dat")
 end do
 close(30)
 close(18)
 
end

!**************************************************************************!
!************************ END *********************************************!
!**************************************************************************!


!**************************************************************************
  !!!                                           SUBROUTINES                                                                    !!!
!**************************************************************************
!**********************************************************************!
!******************************************************************
!*
!**counting lines
!*
!**********************************************************************!
subroutine countlines(nlines)
nlines = 0 
open (1, file = 'surface_exmp3.dat')
do
  read(1,*,iostat=io)
  if (io/=0) exit
  nlines = nlines + 1
end do
 close (1)
end subroutine countlines
!******************************************************************
!*
!**clean array for each cell
!*
!**********************************************************************!
subroutine clean_cell(arr)
integer, parameter :: ny=51,npt=16,nphi=40
real,dimension(ny+1,npt,nphi)::arr
do k=1,ny+1
do i=1,npt
do j=1,nphi
arr(k,i,j)=0.00
end do
end do
end do
end subroutine clean_cell
!******************************************************************
!*
!**clean array for each cell
!*
!**********************************************************************!
subroutine clean_cell1(arr,n)
real,dimension(n)::arr
do i=1,n
arr(i)=-10**8
end do
end subroutine clean_cell1
!******************************************************************
!*
!**clean array for each cell
!*
!**********************************************************************!
subroutine clean_cell2(arr)
integer, parameter :: npt=16,nphi=40
real,dimension(npt,nphi)::arr
do i=1,npt
do j=1,nphi
arr(i,j)=0.00
end do
end do
end subroutine clean_cell2
!**********************************************************************!
!******************************************************************
!*
!**random seed
!*
!**********************************************************************! 
 subroutine change_seed()
 integer :: i, n, clock
 integer, dimension(:), allocatable :: seed
 call random_seed(size = n)
 allocate(seed(n))
 call system_clock(count=clock)
 seed = clock + 37 * (/ (i - 1, i = 1, n) /)
 call random_seed(put = seed)
 deallocate(seed)
 end subroutine change_seed
!**********************************************************************!
 !**********************************************************************!
!******************************************************************
!*
!**gives prob
!*
!**********************************************************************!

subroutine giveprob(phi,pt,rap,spectra,prob,max_val)
real::phi,pt,rap,pi=dacos(-1.0d0),min_pt=0.01,max_pt=3.0,prob,val,max_val
integer, parameter :: npt=16,nphi=40,ny=51
integer::pt_index,phi_index,phi_index1,y_index
real,dimension(npt)::pt_spacing
real,dimension(ny+1,npt,nphi)::spectra
real::ymax=5.0
   do i=1,npt-1
    pt_spacing(i)=givept(i+1,min_pt,max_pt,npt)-givept(i,min_pt,max_pt,npt)
   end do
   pt_spacing(npt)=pt_spacing(npt-1)
    y_index=giveyindex(rap,-ymax,ny)
    phi_index=givephiindex(phi,nphi)
    if (phi_index.eq.40)then
        phi_index1=0
    else
       phi_index1=phi_index+1
    end if          

    pt_index=giveptindex(pt,min_pt,max_pt,npt)

    rap = scaleToZeroToOne(rap, giverapidity(y_index,-ymax,ny), 2.*abs(ymax)/(ny))
    phi = scaleToZeroToOne(phi, givephi(phi_index,nphi),2.*pi/nphi)
    pt = scaleToZeroToOne(pt, givept(pt_index,min_pt,max_pt,npt),pt_spacing(pt_index))
    
   val= trilinear(rap,pt,phi,spectra(y_index,pt_index,phi_index),spectra(y_index,pt_index,phi_index1),&
                        spectra(y_index,pt_index+1,phi_index),spectra(y_index,pt_index+1,phi_index1),&
                        spectra(y_index+1,pt_index,phi_index),spectra(y_index+1,pt_index,phi_index1),&
                        spectra(y_index+1,pt_index+1,phi_index),spectra(y_index+1,pt_index+1,phi_index1))

    prob=exp(val-max_val)
end subroutine giveprob

!**********************************************************************!
!******************************************************************
!*
!**linear interpolation
!*
!**********************************************************************!
function linear(x,v0,v1)
real::x,v0,v1
linear= (1.0-x)*v0 + x*v1
end function linear
!**********************************************************************!
!******************************************************************
!*
!**bilinear interpolation
!*
!**********************************************************************!  
function bilinear(x, y, v00, v01,v10, v11)
real::x, y, v00, v01,v10, v11
bilinear= (1.0-x)*(1.0-y)*v00 + x*(1.0-y)*v10 + (1.0-x)*y*v01 + x*y*v11
end function bilinear
!**********************************************************************!
!******************************************************************
!*
!**trilinear interpolation
!*
!**********************************************************************!  
function trilinear( x,  y,  z,  v000, v001,  v010,  v011, v100,  v101, v110, v111)
real::x,  y,  z,  v000, v001,  v010,  v011, v100,  v101, v110, v111
trilinear= (1.0-x)*(1.0-y)*(1.0-z)*v000 + x*(1.0-y)*(1.0-z)*v100 +(1.0-x)*y*(1.0-z)*v010 + &
x*y*(1.0-z)*v110 + (1.0-x)*(1.0-y)*z*v001 +x*(1.0-y)*z*v101 + (1.0-x)*y*z*v011 + x*y*z*v111 
end function trilinear
!**********************************************************************!
!******************************************************************
!*
!**scale 0 to1
!*
!**********************************************************************! 
function scaleToZeroToOne(value, value_at_zero,interval)
real::value, value_at_zero,interval
scaleToZeroToOne=(value-value_at_zero)/interval
end function scaleToZeroToOne
!**********************************************************************!
!******************************************************************
!*
!**gives pT index
!*
!**********************************************************************! 
function giveptindex(pt,min_pt,max_pt,npt)
real::pt,max_pt,min_pt
integer::npt
 giveptindex=int(sqrt((pt-min_pt)/(max_pt-min_pt)*float((npt-1)**2)))+1
end function giveptindex
!**********************************************************************!
!******************************************************************
!*
!**gives phi index
!*
!**********************************************************************! 
function givephiindex(phi,nphi)
real::phi,dphi
integer::nphi
real::pi=dacos(-1.d0)
dphi=2.*pi/nphi
 givephiindex=int(phi/dphi)+1
end function givephiindex
!**********************************************************************!
!******************************************************************
!*
!**gives y index
!*
!**********************************************************************! 
function giveyindex(rap,ymin,ny)
real::rap,deltay,ymin
integer::ny
deltay=2.0*abs(ymin)/ny
giveyindex=int((rap-ymin)/deltay)+1
end function giveyindex
!**********************************************************************!
!******************************************************************
!*
!**gives rapidity
!*
!**********************************************************************! 
function giverapidity(i,ymin,ny)
real::ymin,deltay
integer::i,ny
deltay=2.0*abs(ymin)/ny
giverapidity= ymin+(i-1)*deltay
end function giverapidity
!**********************************************************************!
!******************************************************************
!*
!**gives pt
!*
!**********************************************************************! 
function givept(i,min_pt,max_pt,npt)
real::min_pt,max_pt
integer::i,npt
givept=min_pt+(max_pt-min_pt)*(float((i-1)**2)/float((npt-1)**2))
end function givept
!**********************************************************************!
!******************************************************************
!*
!**gives phi
!*
!**********************************************************************! 
function givephi(i,nphi)
real::dphi
integer::i,nphi
real::pi=dacos(-1.d0)
dphi=2*pi/nphi
givephi= (i-1)*dphi
end function givephi

