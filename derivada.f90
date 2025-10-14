program derivada

        ! Programa: se presentan los arreglos para aproximar las primeras dos derivadas
        ! de una función real por medio de la Regla Central, tanto la de 3 puntos como 
        ! la de 5. Además, se agrega el valor exacto para su comparación. 

        implicit none

        ! Declaración de variables
        integer :: i           ! Variables contadoras (en Fortran se puede repetir
                                   ! la letra para ciclos diferentes porque lo que hace
                                   ! el lenguaje es REASIGNAR.
        integer :: n               ! Número de pasos
        real :: h               ! Separación
        !dimension es el análogo a tener un vector (renglón, columnas)
        ! allocatable (alojable) guarda un espacio (tamaño desconocido) en la memoria
        real, allocatable, dimension(:)   :: x             ! Elementos del dominio
        real, allocatable, dimension(:)   :: f             ! Función conocida
        real, allocatable, dimension(:)   :: df_exact, ddf_exact      ! Derivadas exactas
        real, allocatable, dimension(:)   :: df_aprox3, df_aprox5, ddf_aprox3, ddf_aprox5      ! Derivadas aproximadas

        ! Valores iniciales
        n = 20
        h = 0.1

        ! allocate (alojar) indica el tamaño de los vectores dado por el usuario
        allocate(x(1:n), f(1:n))
        allocate(df_exact(1:n), ddf_exact(1:n))
        allocate(df_aprox3(2:n-1), df_aprox5(3:n-2))
        allocate(ddf_aprox3(2:n-1), ddf_aprox5(3:n-2))

        ! Ciclo para asignar los valores de la funcion.
        do i = 1, n
                x(i) = real(i) * h      ! Equivale al tiempo
                f(i) = sin(x(i))        ! Evaluacion de la funcion en la varibale independiente.
        end do

        write(*,*) "Primera Derivada."

        write(*,19) "i", "Tiempo", "Exacto", "Aprox 3", "Aprox 5"

        ! Ciclo para la primera derivada.
        do i = 1, n
                ! Derivada exacta.
                df_exact(i) = cos(x(i))                
                if (1 < i .and. i < n)then
                ! Aproximacion por la regla central de 3 puntos.
                        df_aprox3(i) = (f(i+1) - f(i-1)) / (2.0*h) 
                end if      
                if (2 < i .and. i < n-1)then
                       ! Aproximacion por la regla central de 5 puntos.
                       df_aprox5(i) = (8.0 * f(i+1) - 8.0 * f(i-1) + f(i-2) - f(i+2)) / (12.0 * h)
                end if      
                ! Salida a la pantalla. 
                write(*,20) i, x(i), df_exact(i), df_aprox3(i), df_aprox5(i)
        end do

        write(*,19) "i", "Tiempo", "Exacto", "Aprox 3", "Aprox 5"
        
        write(*,*) "Segunda Derivada."

        ! Ciclo para la segunda derivada.
        do i = 1, n
                ! Derivada exacta.
                ddf_exact(i) = -sin(x(i))       
                if (1 < i .and. i < n)then
                        ! Aproximacion por 3 puntos.
                        ddf_aprox3(i) = (f(i+1) - 2.0 * f(i) + f(i-1)) / (h**2)         
                end if
                if (2 < i .and. i < n-1)then
                        ! Aproximacion por 5 puntos.
                        ddf_aprox5(i) = (- f(i+2) + 16.0 * f(i+1) - 30.0 * f(i) + 16.0 * f(i-1) &
                               & - f(i-2)) / (12.0 * h**2) 
                end if
                ! Salida a la pantalla. 
                write(*,20) i, x(i), ddf_exact(i), ddf_aprox3(i), ddf_aprox5(i)
        end do 

        deallocate(x, f, df_exact, ddf_exact)
        deallocate(df_aprox3, df_aprox5, ddf_aprox3, ddf_aprox5)

        ! Escritura con formato dirigido. 
        ! Es la instruccion para indicar el tipo de formato de salida. 
        ! La sintaxis es la siguiente.
        ! rFwd --> r : repeticion; F : formato; w : numero de caracteres; d : numero de digitos.
        19 format(5A15)
        20 format(1I15, 4F15.8)

end program derivada
