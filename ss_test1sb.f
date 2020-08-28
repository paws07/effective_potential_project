	implicit double precision(a-h,o-z),integer(i-n)
	double precision r(10,3),v(10,3),a(10,3),m(10)
	character*70 filout
	dimension rl(5,2),rctr(3),rpt(3)
c
	call getinp(nbod,r,v,m,bsize,npts,filout)
c
 1	format(3d18.10)
c
	  open(75,file=filout)
c Findng L4 and L5 (assume bodies ordered M1,M2 ...)
	cosv=r(2,1)/(sqrt(r(2,1)**2+r(2,2)**2))
	sinv=r(2,2)/(sqrt(r(2,1)**2+r(2,2)**2))
	
	rearth=sqrt(r(2,1)**2+r(2,2)**2+r(2,3)**2)
	
	rlfourx=rearth/2.0*(m(1)-m(2))/(m(1)+m(2))
	rlfoury=rearth*sqrt(3.0)/2.0
	write(*,*) 'L4 Point: ', (rlfourx*cosv-rlfoury*sinv),(rlfourx*sinv+rlfoury*cosv)
	
	rlfivex=rearth/2.0*(m(1)-m(2))/(m(1)+m(2))
	rlfivey=rearth*-sqrt(3.0)/2.0
	write(*,*) 'L5 Point: ', (rlfivex*cosv-rlfivey*sinv),(rlfivex*sinv+rlfivey*cosv)
	
	rctr(1)=(0,0)
	rctr(2)=(0,0)
	rctr(3)=0.0
c
	rpt(3)=0.0
	do 60 ix=1,npts
	  rpt(1)=rctr(1)-0.5*bsize+(1.*ix-1.0)/(1.*npts-1.0)*bsize
	  do 61 iy=1,npts
	    rpt(2)=rctr(2)-0.5*bsize+(1.*iy-1.0)/(1.*npts-1.0)*bsize
	    call gravp(nbod,v,rpt,r,m,potl)
	    write(75,1) rpt(1),rpt(2),potl
 61	  continue
 60	continue
	close(75)
	stop
	end
c-----------------------------------------------------------
	subroutine gravp(nbod,v,rpt,r,m,potl)
	implicit double precision(a-h,o-z),integer(i-n)
	double precision rpt(3),r(10,3),m(10)
c G in AU,day from original ss.f code via wikipedia on G
	g=1.488180714e-34
	
c a= the distance between two bodies
	a=sqrt((r(1,1)-r(2,1))**2+(r(1,2)-r(2,2))**2)
	cosiv=r(2,1)/(sqrt(r(2,1)**2+r(2,2)**2))
	sinev=r(2,2)/(sqrt(r(2,1)**2+r(2,2)**2))
c om=Reduced Mass Ratio
	om=m(1)/(m(1)+m(2))
c compute potential of the body at the center. 
	term2=-om/sqrt((rpt(1)/a+cosiv*(1-om))**2+(rpt(2)/a+sinev*(1-om))**2)
c centrifugal potential due to only earth and sun
	term3=-0.5*(((rpt(1)**2)/a**2)+((rpt(2)**2)/a**2))
	potl=0
	do 110 inb=2, nbod
	  cosiv=r(inb,1)/(sqrt(r(inb,1)**2+r(inb,2)**2))
	  sinev=r(inb,2)/(sqrt(r(inb,1)**2+r(inb,2)**2))
	  a=sqrt((r(1,1)-r(inb,1))**2+(r(1,2)-r(inb,2))**2)
c om=Reduced Mass Ratio
	  om=m(1)/(m(1)+m(inb))
c compute m2/abs(r-r1)
	  term1=(-1+om)/sqrt((rpt(1)/a-cosiv*om)**2+(rpt(2)/a-sinev*om)**2)
	  potl=potl+term1
 110	continue
	potl=potl+term2+term3
	return
	end
c-----------------------------------------------------------
	subroutine getinp(nbod,r,v,m,bsize,npts,filout)
	implicit double precision(a-h,o-z),integer(i-n)
	double precision r(10,3),v(10,3),m(10)
	character*70 junk,fil,filout
c
	read(*,*) junk
	read(*,*) nbod
	do 10 ibod=1,nbod
	  read(*,*) fil
	  open(25,file=fil)
	    read(25,*) m(ibod)
	    read(25,*) r(ibod,1),r(ibod,2),r(ibod,3)
	    read(25,*) v(ibod,1),v(ibod,2),v(ibod,3)
	  close(25)
 10	continue
	read(*,*) junk
	read(*,*) filout
	read(*,*) bsize
	read(*,*) npts
	return
	end
c------------------------------------------------------------
