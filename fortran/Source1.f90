
    program vario_simples
    
    real :: coordx(300), coordy(300), coordz(300), racio(300), Pop(300)
    real :: dist(300,300), w(300,300),gama1(300,300)
    real :: lag, soma1(60), soma2(60), gama2(60), distot(60), media
    real :: poptotal, varw, varp, variance
    integer :: var2,i,nv, nh(60),nobs
    
 !   Character(LEN=20) :: var1,var3,var4,var5,var6,var7
    
    OPEN (10,FILE='c:/VARIO/DATA/deathRacioPop20200511.out',STATUS='OLD')
    
    OPEN (20,FILE='c:/VARIO/OUTPUT/vadeathRacio.out',STATUS='unknown')
    OPEN (30,FILE='c:/VARIO/OUTPUT/debug.out',STATUS='unknown')

 !   READ (10, *) var1
 !   write (*,*) var1
 !   READ (10, *) var2
 !   write (*,*) var2
 !   READ (10, *) var3
 !   write (*,*) var3
 !   READ (10, *) var4
 !   write (*,*) var4
 !   READ (10, *) var5
 !   write (*,*) var5
 !   READ (10, *) var6
 !   write (*,*) var6
 !   READ (10, *) var7
 !   write (*,*) var7
    
    
    nobs=278
    
    do i=1, nobs
        READ (10, *) coordx(i), coordy(i), coordz(i), racio(i), Pop(i)
!        write (30,*) coordx(i), coordy(i), coordz(i), racio(i), Pop(i)
    end do
    media=sum(racio)/nobs
   
    variance=0
    poptotal=0
    
    do i=1,nobs
        varw=pop(i)*pop(i)/(pop(i)+pop(i))
        varp=varw*((racio(i)-media)**2)
        poptotal=poptotal+varw
        variance=variance+varp
        do j=1, nobs
	       dist(i,j)=sqrt((coordx(i)-coordx(j))**2+(coordy(i)-coordy(j))**2)
           w(i,j)=pop(i)*pop(j)/(pop(i)+pop(j))
            gama1(i,j)=w(i,j)*((racio(i)-racio(j))**2)
           write(30,*) i, j, dist(i,j), w(i,j), gama1(i,j)
        end do
    end do
   
   variance=variance/poptotal
   
   lag=5000
   nv=20
   do k=1,nv
       nh(k)=0
       soma1(k)=0
       soma2(k)=0
       distot(k)=0
   end do
       
   do i=1,nobs
        do j=1,nobs
            do k=1,nv
                 if ((dist(i,j)<=(lag*k)) .and. (dist(i,j)>lag*(k-1))) nh(k)=nh(k)+1
                 if ((dist(i,j)<=(lag*k)) .and. (dist(i,j)>lag*(k-1))) soma1(k)=soma1(k)+gama1(i,j)
                 if ((dist(i,j)<=(lag*k)) .and. (dist(i,j)>lag*(k-1))) distot(k)=distot(k)+dist(i,j)
                 if ((dist(i,j)<=(lag*k)) .and. (dist(i,j)>lag*(k-1))) soma2(k)=soma2(k)+w(i,j)
           end do          
       end do
   end do
   
   write(20,*) variance
   do k=1,nv
       distot(k)=distot(k)/nh(k)
       gama2(k)=(1/(2*soma2(k)))*soma1(k)
       write(20,*) nh(k)/2, distot(k), gama2(k)
   end do
   
   end program vario_simples

