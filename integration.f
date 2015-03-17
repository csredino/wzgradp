	subroutine integration(n,it1,n1,it2,n2,ierr,n_proc,my_id)
c
c
c
c	n: number of dimensions
c	it1: number of iterations in the grid search
c	n1: number of points for each dimensions in the grid search
c	it2: number of iterations in the precision calculation(<50))
c	n2: number of points for each dimensions in the precision 
c			calculation
c	note: -in the grid search alpha is progressively adjusted
c	      -in the precision calculation alpha=0.
c	      -ihist is used to decide if one wants to plot things
c	      -reno is vegas with ranmar as a random number generator
c
c
	implicit none	
	integer n,it1,n1,it2,n2
	integer ncall1,ncall2
	real*8 alph
	integer ndmx,mds
	common/bveg3/alph,ndmx,mds
	real*8 kern
	external kern
	real*8 stot,dstot,chi2tot,dstotper,stoti(50),sdeli(50)
c
c All communicated varaibles are of same type, so lets pack them into a small vector 
c
	real*8 Sresults(3),results(3)
c
	integer ndim,it,ihist
c
c       new variables for parallel version
c
	integer ierr,n_proc,my_id
c
	common/par_integration/ndim,ihist,it,ncall2
	common/par_res_integration/stot,dstot
	common/par_integration_reno/stoti,sdeli

	real*8 s_avg,ds,sdel_avg,dsdel
c	real*8 a,b,c
	integer i
c
c MPI header to to allow for parallel processing
c
	include 'mpif.h'
c
c the number of calls to reno will be rescaled for parallelization, so that the size of the actual task remains the same when the number of processors increases
c need to be careful to avoid rounding errors, number of calls is typically a  multiple of 100, so 1,2,5,10,20,50,100 processors are safe.	
	
c
c
	ncall1=int(n*n1/n_proc)
	ncall2=int(n*n2/n_proc)
	ndim=n
	ihist=0
	alph = 1.5d0
	if(it1.eq.0)then
	  if (my_id == 0) then
	  print*,' No grid searching!'
	  endif
	else 
	  if (my_id == 0) then
	  print*,' Grid searching:'
	  endif
	  call reno(it1,ncall1,n,kern,stot,dstot,chi2tot,ierr,n_proc,my_id)
	end if
	alph=0.d0
	ihist=1
	it=it2
	if (my_id == 0) then
	print*,' Precision calculation:'
	endif
c the precision calculation will run in parallel, each with their own set of random numbers
c
	call reno1(it2,ncall2,n,kern,stot,dstot,chi2tot,ierr,n_proc,my_id)
c
c
c results from precision calculation are reduced(summed and averaged)stot,dstot,chi2tot
c
	results(1)=stot
	results(2)=dstot
	results(3)=chi2tot
	call  MPI_REDUCE(results,Sresults,3,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD,ierr)

c
	stot=dble(Sresults(1)/n_proc)
	dstot=dble(Sresults(2)/n_proc)
	chi2tot=dble(Sresults(3)/n_proc)
c
c
c
 	dstotper=100*dstot/stot
c
	if(my_id ==0) then
	write(8,20)stot,dstot,dstotper,chi2tot
20	format(/,' Sig from Reno:',E14.8,' +- ',E10.4,'(',E10.4,'%)
	1    chi2:',E8.2)
	endif
c
c       averages from the it2 shots
c
	s_avg=0.d0
	ds=0.d0
	sdel_avg=0.d0
	dsdel=0.d0
	do i=1,it2
	   s_avg=s_avg+stoti(i)
	   ds=ds+stoti(i)**2
	   sdel_avg=sdel_avg+sdeli(i)
	   dsdel=dsdel+sdeli(i)**2
	end do
	s_avg=s_avg/dble(it2)
	sdel_avg=sdel_avg/dble(it2)
	ds=ds/dble(it2)
	dsdel=dsdel/dble(it2)
	ds=sqrt(dble(it2)*(ds-s_avg**2)/(dble(it2)-1.d0))
	dsdel=sqrt(dble(it2)*(dsdel-sdel_avg**2)/(dble(it2)-1.d0))
c
	if (my_id == 0) then
	write(8,21)s_avg,ds/sqrt(dble(it2)), 100*ds/sqrt(dble(it2))/s_avg
 21	format(/,' Sig from averages:',E14.8,' +- ',E10.4,'(',E10.4,'%)')
	write(8,22)sdel_avg,dsdel,ds 
 22	format(/,' Average shot uncertainty:',E10.4,' +-',E10.4,'(',E10.4,')')
	endif

c
c       chi2 fit (gives same result as reno)
c
c	a=0.d0
c	b=0.d0
c	c=0.d0
c	do i=1,it2
c	   a=a+1.d0/sdeli(i)**2
c	   b=b+stoti(i)/sdeli(i)**2
c	   c=c+stoti(i)**2/sdeli(i)**2
c	end do
c	s_avg=b/a
c	ds=sqrt(1.d0/a)
c	chi2tot=(s_avg**2*a-2.d0*s_avg*b+c)/(dble(it2)-1.d0)
c
	if (my_id == 0) then
c	write(8,23)s_avg,ds, 100*ds/s_avg,chi2tot
c 23	format(/,' Sig from chi2 fit:',E14.8,' +- ',E10.4,'(',E10.4,'%)  chi2:',E8.2)
	endif

	return
	end

	subroutine setup_integration(i,my_id)
c
c	can be used to run several jobs with different seeds
c	ranmar is the random number generator used in reno(=vegas)
c
	implicit none
	integer i
	integer*4 my_id
	integer iseed,n1,n2
	if (my_id == 0) then
	if(i.eq.0)then
	      print*, ' Seed for the random number generator RANMAR:'
	      read*,iseed
	else
	  iseed=i
	endif
	endif
	call rmarin(iseed,n1,n2)
	return
	end
