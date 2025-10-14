program Taylor
! Los metodos numericos para aproximar la solucion de una EDO de primer orden están relacionados
! con la aproximacion de Taylor de la funcion. 
! Aqui se veran los metodos de Taylor de primer y segundo orden. 
! Autor: Esteban Juarez           Fecha: 2025
        implicit none

        integer :: i
        integer :: n            ! Numero de pasos
        real    :: x            ! Variable dependiente para Taylor de primer orden
        real    :: u            ! Variable dependiente para Taylor de segundo orden
        real    :: t            ! Variable independiente
        real    :: x0, t0       ! Condiciones iniciales
        real    :: h, xe        ! Incremento y solucion exacta
        real    :: r, k         ! Parametros de la ecuacion logistica

        ! Archivo de entrada
        open(100, file='Taylor_data.inp')        
                read(100,*)
                read(100,*)n, t0, x0, h
                read(100,*)
                read(100,*)r, k
        close(100)

        x = x0          ! Posicion para Taylor de primer orden
        u = x0          ! Posicion para Taylor de segundo orden
        t = t0
        
        ! Archivos de almacenamiento
        open(200,file="apsol_Taylor.dat")
        open(300,file="exsol_Taylor.dat")

        write(200,20)t, x, u
        write(300,20)t, x
 
        do i = 1, n
                x = x + funcion(x, r, k) * h                                               ! Formula de Taylor primer orden
                u = u + (funcion(u, r, k) * h) + ((dfunc(u, r, k) / 2.0) * h**2)     ! Formula de Taylor de segundo orden
                t = t + h                                                                  ! Avance en la variable independiente
                xe = sol(t, t0, x0, r, k)                                                  ! Funcion con la solucion exacta
                write(200, 20)t, x, u
                write(300, 20)t, xe
        end do

        close(200)
        close(300)

        ! Formato de salida
        20 format(3F15.8)

        ! Definicion de funciones, cambiara la ecuacion diferencial segun la naturaleza del problema. 

        contains
  
                function funcion(x, r, k)
                        implicit none
                        real :: x, r, k, funcion
                        funcion = r * x * (1.0 - x/k)           ! Ecuacion Logistica
                        return
                end function funcion

                ! Derivada de la funcion de arriba
                function dfunc(x, r, k)
                        implicit none
                        real    :: x, r, k, dfunc
                        dfunc = r**2 * x * (1.0 - 2.0 * x /k) * (1.0 - x / k)
                end function dfunc

                ! Solucion exacta
                function sol(t, t0, x0, r, k)
                        implicit none
                        real :: t, t0, x0, r, k, sol
                        sol = (k - x0) * exp(-r * (t - t0)) / x0
                        sol = 1.0 + sol
                        sol = k / sol
                        return
                end function sol

end program Taylor
