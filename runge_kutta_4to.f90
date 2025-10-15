program rungekutta4
        !
        ! Método de Runge-kutta de cuarto orden para resolver numéricamente la ecuación logística.
        ! Los parámetros se definen de la siguiente forma:
        ! k1 = f(x, t)          k2 = f(x + k1*h/2.0, t + h/2.0) 
        ! k3 = f(x + k2*h/2.0, t + h/2.0)               k4 = f(x + k*3*h, t + h)
        ! Finalmente, la pendiente formada con estos parámetros lo sustituimos 
        !       x = x + h * f((k1 + 2*k2 + 2*k3 + k4)/6.0)
        ! Autor                 Semestre
        ! Esteban J             2026-1
        !

        implicit none

        integer :: i, n
        real    :: t, x, t0, x0, h, xe
        real    :: r, k
        real    :: k1, k2, k3, k4               ! Son los parámetros del metodo

        open(100, file='rk_4to.inp')
                read(100,*)
                read(100,*)n, t0, x0, h, r, k
        close(100)

        ! Condiciones iniciales
        x = x0
        t = t0

        open(200, file='apsol_logistica.dat')
        open(300, file='exsol_logistica.dat')

        write(200,*)t, x
        write(300,*)t, x

        ! Ciclo para asignar iterativamente los parámetros.
        do i = 1, n
                k1 = flog(x, t, r, k)
                k2 = flog(x+k1*h/2.0, t+h/2.0, r, k)
                k3 = flog(x+k2*h/2.0,t+h/2.0, r, k)
                k4 = flog(x+k3*h,t+h, r, k)
                x = x + h * (k1+2.0*k2+2.0*k3+k4)/6.0
                t = t + h
                xe = sol(t, t0, x0, r, k)
                write(200,20)t, x
                write(300,20)t, xe
        end do

        20 format(2F15.8)

        contains

                function flog(xx,tt,rr,kk) result(yy)
                        implicit none
                        real :: xx, tt, rr, kk, yy

                        yy = rr * xx * (1.0 - xx/kk)

                 end function flog

                 function sol(tt, tt0, xx0, rr, kk) result(yy)
                         implicit none
                         real :: tt, tt0, xx0, rr, kk, yy

                         yy = (kk-xx0)*exp(-rr*(tt-tt0))/xx0
                         yy = 1.0 + yy
                         yy = kk/yy
                 end function sol

end program rungekutta4
