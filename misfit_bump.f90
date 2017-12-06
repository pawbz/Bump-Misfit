
￼
!\end{comment}
!\begin{sourcecode}
subroutine fg_bump(p,      & ! 
                        eps,     & !
                        p0,     & !
                        q0,     & !
                        q0freq, & !
                        J,      & !
                        gp,     & !
                        n1pow2, & !
                        n2pow2, & !
                        n1,     & !
                        n2,     & !
                        n3,     & !    
                        saved,  & !
                        i3save  & !
                        ) !\end{sourcecode}
!\begin{comment}
implicit none
real, intent(in)                        :: p(n1*n2*n3), eps,  p0(n1*n2*n3)
real, intent(in)                        :: q0freq(n1pow2*n2pow2), q0(n1*n2*n3) 
real, intent(out)                       :: gp(n1*n2*n3)
integer, intent(in)                     :: n1, n2, n3, n1pow2, n2pow2
real      , intent(out)                 :: J
integer, intent(in)                     :: i3save
real, intent(inout)                     :: saved(n1*n2*n3*5)

integer                                 :: i2, i3, ap1, ap2, i1, aphat1, aphat2
real                                    :: Jtemp

complex, allocatable                    :: pbump(:,:), p0bump(:,:), dJdpbump(:,:) 
complex, allocatable                    :: temp2dcomplex(:,:)
real, allocatable                       :: temp1D(:) 

! gradient w.r.t p
gp = rzero_de

! functional output to zero
J = rzero_de

! allocate temporary arrays in f-k domain
allocate(p0bump(n1pow2, n2pow2), pbump(n1pow2, n2pow2))

! other temp
allocate(temp1D(n1))

! other temp for gradient
allocate(dJdpbump(n1pow2, n2pow2))

! loop over sources
do i3 = 1, n3
        ! initialize f-k domain
        p0bump = cmplx(rzero_de, rzero_de)
        pbump = cmplx(rzero_de, rzero_de)

        ! zero padding p and p0 
        do i2 = 1, n2
                ap1 = 1 + (i2-1) * n1 + (i3-1) * n1 * n2;   
                ap2 =   (i2) * n1 + (i3-1) * n1 * n2 ;
        
                ! absolute value of the observed and modeled data
                ! after applying weighting q0
                temp1D = q0(ap1:ap2) * sqrt((p0(ap1:ap2))**2 + eps);
                call nlag_npow2_pad_truncate(x = temp1D, &
                         xpow2_complex = p0bump(:,i2), nplags = n1-1, nnlags = 0,&
                        npow2 = n1pow2, flag = 1)

                temp1D = q0(ap1:ap2) * sqrt(( p(ap1:ap2))**2 + eps);
                call nlag_npow2_pad_truncate(x = temp1D, &
                                 xpow2_complex = pbump(:,i2), &
                                nplags = n1-1, nnlags=0, npow2 = n1pow2, flag = 1)
        enddo

        ! cfft2 forward
        call cfft2_fftw(p0bump, n1pow2, n2pow2, 1)
        call cfft2_fftw(pbump, n1pow2, n2pow2, 1)

        ! functional -- q0freq helps to choose low frequencies
        Jtemp = rzero_de;
        ! multiplication with q0 will apply blurring
        do i2 = 1, n2pow2
                do i1 = 1, n1pow2
                        Jtemp = Jtemp + real(q0freq(i1 + (i2-1)*n1pow2)) * &
                                real((abs(pbump(i1, i2) - p0bump(i1, i2)))**(2.0e0))
                enddo
        enddo
        J = J + Jtemp

        forall(i2=1:n2pow2,i1=1:n1pow2)
                        dJdpbump(i1,i2) = &
                        cmplx(2.0e0, rzero_de) * cmplx(q0freq(i1 + (i2-1)*n1pow2), rzero_de) * &
                        (pbump(i1,i2) - p0bump(i1,i2))
        endforall
        
        ! inverse cffts
        call cfft2_fftw(dJdpbump, n1pow2, n2pow2, -1)

        ! scaling
        dJdpbump = dJdpbump / cmplx(real(n1pow2),0.0) / cmplx(real(n2pow2),0.0)

        do i2 = 1, n2
                ap1 = 1 + (i2-1) * n1 + (i3-1) * n1 * n2;   
                ap2 =   (i2) * n1 + (i3-1) * n1 * n2 ;
                call nlag_npow2_pad_truncate(x = gp(ap1:ap2), &
                         xpow2_complex = dJdpbump(:,i2), nplags = n1-1, &
                        nnlags=0, npow2 = n1pow2, flag = -1)
        
                do i1 = ap1, ap2
                        gp(i1) =  gp(i1) * p(i1) * q0(i1) * (sqrt((p(i1))**2 + eps))**(-1.e0)
                enddo
        enddo

        if(i3 .eq. i3save) then
                ! p
                saved(n1*n2+1:2*n1*n2) = p(n1*n2*(i3save-1)+1:n1*n2*i3save)

                ! p0
                saved(1:n1*n2) = p0(n1*n2*(i3save-1)+1:n1*n2*i3save)

                ! pbump
                allocate(temp2dcomplex(n1pow2, n2pow2))
                forall(i2=1:n2pow2,i1=1:n1pow2)
                        temp2dcomplex(i1,i2) =  pbump(i1,i2) * &
                        cmplx(q0freq(i1 + (i2-1)*n1pow2), rzero_de)
                endforall
                call cfft2_fftw(temp2dcomplex, n1pow2, n2pow2, -1)
                saved(3*n1*n2+1:4*n1*n2) = pack(real(temp2dcomplex(1:n1,1:n2)),.true.)
                
                ! p0bump
                forall(i2=1:n2pow2,i1=1:n1pow2)
                        temp2dcomplex(i1,i2) =  p0bump(i1,i2) * &
                        cmplx(q0freq(i1 + (i2-1)*n1pow2), rzero_de)
                endforall
                call cfft2_fftw(temp2dcomplex, n1pow2, n2pow2, -1)
                saved(2*n1*n2+1:3*n1*n2) = pack(real(temp2dcomplex(1:n1,1:n2)),.true.)
                deallocate(temp2dcomplex)

                ! gp
                saved(4*n1*n2+1:5*n1*n2) = gp(n1*n2*(i3save-1)+1:n1*n2*i3save)

        endif
enddo
if(allocated(temp1D)) deallocate(temp1D)
deallocate(pbump, p0bump)
if(allocated(dJdpbump)) deallocate(dJdpbump)

end subroutine fg_bump
