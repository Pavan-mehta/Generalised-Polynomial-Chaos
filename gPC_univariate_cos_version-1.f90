
program gpc

!************************************************************************************************************************************************************!
! This is a Fortran 90 code which can be complied in any standard Fortran complier such a gfortran (Linux - my personal favrate), NAG or Siverfrost
!
! This code is for an univariate non - intrusive Generalised Polynomial Choas Technique, also known as Stoachastic Galerkin Method 
!
! While this is for an univariate case it can be readily extended to multivariate case. 
!
! The test function used here is F(x) = cos (x); x in RADIANS. However, you can change it as per your will, even taking real world data
!
! I'm a Thermo-fluids Engineer not a fulltime programmer. Hence, this code may not be efficient but it does the job :)
!
! My name is Pavan Mehta, I wrote this while my MSc Dissertation at The University of Manchester, UK in summer 2016. Needless to mention validated as well
!
! I hope you find it usefull. Please don't hesitate to contact me, I'll try my best to answer your queries. Email: mehtapavanp@gmail.com
!
! Constructive critisim always welcome
!
! Date of publication: 28 March 2017
!************************************************************************************************************************************************************!


implicit none
integer y,m
parameter (y=6) !
real Qp(y),w(y)


double precision :: x, x1, x2, p1, p2, pp, z, z1, xm, p3
double precision :: Tmin, Tmax, XT(y), CT, MT, F(y)	
integer::i,j, k, Plt

parameter (Plt = 20)
real::P(y+1,y), P_Plt(y+1,Plt+1)

double precision EPS, P_proj2(y), SumP,  S_gPC_Cof(y), Sum_S_gPC, Value, tmp1
parameter (EPS=3.d-14)


!Ya! there are alot of variables, but you need not worry about it :)
!The key variables are y, Plt and x

!y - is the total number of roots (Quad Points). These roots are obtained using Gaussian Quadrature. 
!From literetaures it is known that for an orthogonal polynomial, the most accurate results are obtained at its roots.

!Tensor grid is constructed 
!Here Plt tells you the maximum number of points
!While x is the step size

!********************************************************************************************************************

!Computing the roots - Gauss Quatrature

!***********************************************************************************************************************


x1=-1.d0
x2=1.d0

!OPEN(unit=5,FILE="Gaus-quad.dat")
!write(5,*) 'Range for P/D, Lower limit -',x1, 'upper linmit -',x2

m = y !as the roots are symmetric, its better sometimes to use m = (n+1)/2

xm = 0.5d0*(x2+x1)
x1 = 0.5d0*(x2-x1)
j = 1


do i=1,m
	z = cos(3.141592654d0*(i-0.25d0)/(y+0.5d0))
	
	1 continue
		p1=1.d0
		p2=0.d0
		do j=1,m
			p3=p2
			p2=p1
			p1 = ((2.d0*j-1)*z*p2-(j-1.d0)*p3)/j
		enddo
		
		pp = m*(z*p1-p2)/(z*z-1.d0)
		z1 = z
		z = z1-p1/pp
	
	
	if(abs(z-z1).gt.EPS) goto 1
	Qp(i)=xm-x1*z
	Qp(y+1-i)=xm+x1*z
	w(i)=2.d0*x1/((1.d0-z*z)*pp*pp)
	w(y+1-i)=w(i)
	!	write(*,*) 'root(',i,')= ',Qp(i), 'w(',i,')= ',w(i)
		
	enddo
	
	!***********************************************************************
	!orthogonal polynomial values (Legender) calculation at quad points in hilbert space
	!**********************************************************************

k = -1
!open(unit=1,file='poly-valida.dat')

!open(unit=2,file='x-poly-valida.dat')
do i = 1,y
	!x = -1
	do j=1,y
 
	x = Qp(j) 

		if (i.eq.1) then
			P(1,j)=1
		!	write(*,*) P(1,j), x, i
		!	write(2,*) x
		else if (i.eq.2) then
			P(2,j)=x
		!	write(*,*) P(2,j), x, i
		else

			P(i,j)=(((2*k)+1)*x*(P(i-1,j))-((k)*P(i-2,j)))/(k+1)	!coeffiecents are computed using regression formula. 

		!	write(*,*) P(i,j) , x, i
		end if
	!	write(*,*) x
	!x = x + 0.1
	end do

	k = k + 1

enddo




	!**************************************************
	!Quad points [test - cos(x) ] x in radians
	!****************************************************
	
!OPEN(unit=6,FILE="Quad.dat")



	Tmin=0.D0
	Tmax=2.D0	!Min and Max value of our real interval
	
	CT= 0.5D0*(Tmin+Tmax)
	MT= 0.5D0*(Tmax-Tmin)

	Do i=1,y
		XT(i)= CT + MT*Qp(i) 	!Transformation is done such that there is one - to - one correspondce
	
	!write(*,*) 'Quad(',i,')=,' ,XT(i)
	
	End do !i
	
	
	!**************************************************
	!Function [test - cos(x) ] x in radians
	!****************************************************
	
	Do j=1,y	
			F(j)=cos(XT(j))	!Here cos can be replaced by any other function or real world data


	!write(*,*) 'Cos values at quad point (',j,')=,' ,F(j)
		End do !j



	!**************************************************
	!polynomial projection (Lengenger)
	!****************************************************

	Do k=1,y
		SumP=0.D0
		Do i=1,y
		SumP= SumP + w(i)*P(k,i)*P(k,i)
		End do ! i
		P_proj2(k)=SumP
	
	!	write(*,*) 'Poly projection(',k,')=,' ,P_proj2(k)
	!	
	End do ! k



	!**************************************************
	!GPC coeefienct polynomial (test - sin wave) y = sin(T)
	!****************************************************
	
	Do k=1, y
	
			Sum_S_gPC=0.D0
			
				Do j=1,y
					Sum_S_gPC= Sum_S_gPC + (w(j)*P(k,j)*F(j))
   
				End do !j
			
			S_gPC_Cof(k)=Sum_S_gPC/((P_proj2(k)))
	
!OPEN(unit=7,FILE="Gpc-coef.dat")
		
!	Write(*,*) 'S_gPC_Cof(',k,')=' ,	S_gPC_Cof(k)	
!	Write(7,*) 'S_gPC_Cof(',k,')=' ,	S_gPC_Cof(k)
	
	End do !k
	

!**************************************************************************************
!**	Genrerating and Plotting the Response surface
!**************************************************************************************	

k = -1
!open(unit=1,file='poly-valida.dat')

!open(unit=2,file='x-poly-valida.dat')
do i = 1,y
	x = -1
	do j=1, (Plt + 1)
 
		if (i.eq.1) then
			P_Plt(1,j)=1
		!	write(*,*) P_Plt(1,j), x, i
		!	write(2,*) x
		else if (i.eq.2) then
			P_Plt(2,j)=x
		!	write(*,*) P_Plt(2,j), x, i
		else

			P_Plt(i,j)=(((2*k)+1)*x*(P_Plt(i-1,j))-((k)*P_Plt(i-2,j)))/(k+1)

		!	write(*,*) P_Plt(i,j), x, i
		end if
	!	write(*,*) x
				
	
	x = x + 0.1	!You can here change the step size; do remember to adjust Plt, which is the max. value

			!If you change the value of "0.01" and the output has abursurdly high values. 
			!It is because the tranformation is made with the interval size of 2 for both Real and Hilbert space
			!Create your own transfromation. You would be good to go!
	
	end do

	k = k + 1

enddo

!**************************************************************************************
! Final step
!**************************************************************************************
	

x = -1

write(*,*) "Univariate non - intrusive Generalised Polynomial Choas Technique, also known as Stoachastic Galerkin Method"
write(*,*) " F(x) = cos(x)"
write(*,*) "x values in Radians" !here this does not mean the above x; it refer's to the angles. You will undertand after you exceute it
write(*,*) "Author: Pavan Mehta"
write(*,*)

Do i=1, (Plt+1)
	
			Value=0.D0
								
			Do k=1,y
			
				  Value = Value + S_gPC_Cof(k)*P_Plt(k,i)
			
			End do

!	write(7,*)'(',i,')',del_i,'(',j,')=',del_j,'Value=',Value

	tmp1=CT+MT*x

	
	write(*,*) 'x =',tmp1,'; F(x) =',Value

!	Assigning Values here for Tecplot files
		!	ValueI(i)=Value
		!	LocI_gPC(i)=x
		!	LocI_Phy(i)=tmp1	!I haven't written code here to write these in a file. But it's not cumbersome to write those values. :)

	
	x = x + 0.1	!Caution. Check if you have changed the value in the previous step as well
	
End do
	

end program gpc !Thank You :) Hope you enjoyed and found useful

