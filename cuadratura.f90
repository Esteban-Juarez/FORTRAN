program cuadratura
! EL valor de una integral definida (una variable) es el punto central de este programa. 
! A través del metodo de cuadraturas, regla del Trapecio y la de Simpson, se hace una 
! estimación de la integral comparandola con el vacor exacto. 
! Juarez, E., E.                Agosto 2025

        implicit none 

        integer :: i                    ! Variable contadora
        integer :: n                    ! Numero de pasos
        real    :: h                    ! Longitud de los intervalos.
        real    :: a                    ! Limite inferior de la integral
        real    :: b                    ! Limite superior de la integral
        real    :: trapecio, simpson    ! Cuadraturas. 
        real    :: factor               ! Es el factor de multiplicacion para el metodo de Simpson.
        real, allocatable, dimension(:) :: x    ! Son los puntos donde se evalúa la función
        real, allocatable, dimension(:) :: f    ! Es la función para integrar

        ! open  sirve para indicarle al compilador que tomará los datos de entrada de un archivo externo. 
        ! A esto se le llama: estructura de lectura de datos almacenados.
        ! La estructura tiene la etiqueta 100 y el nombre del archivo de lectura es "input_data.inp".
        open(100, file="input_data_cuadratura.inp")         ! Lectura de los datos de entrada
                read(100, *)n, a, b
        close(100)

        ! Si el número de intervalos es un número impar se sigue con el método de cuadratura. De lo contrario se salta el ciclo.
        if(mod(n, 2)==0)then            ! Paso de control, n debe ser impar
                write(*,*)"Se ha proporcionado n par y debería ser impar"
                write(*,*)"La ejecución se finalizará"
                stop            ! La instrucción stop es similar al break de Python.
        end if

        ! Cuando n es impar se alojan los arreglos vectoriales de la variable independiente y la función. 
        allocate(x(1:n), f(1:n))

        ! Se determina la longitud de los subintervalos.
        h = (b-a)/real(n-1)

        ! Ciclo para asignar los valores a la funcion
        do i = 1, n
                ! Comenzamos en el punto "a" y despues avanzamos una cantidad "h" hasta "b"
                x(i) = a + real(i-1)*h
                ! Integrando.
                f(i) = exp(x(i))
        end do

        ! Se calcula la integral en los extremos para ambos metodos.
        trapecio = (f(1) + f(n))/2.0
        simpson = f(1) + f(n)

        ! Ciclo para calcular las contribuciones restantes.
        do i = 2, n-1
                ! Regla del trapecio.
                trapecio = trapecio + f(i)
                ! Regla de Simpson
                if(mod(i,2) == 0)then
                        factor = 4.0
                else
                        factor = 2.0
                end if
                simpson = simpson + factor * f(i)
        end do

        ! Multiplicacion por la longitud del intervalo.
        trapecio = h * trapecio
        simpson = h * simpson / 3.0

        ! Salida de los resultados
        write(*,20) "Exacto", "Trapecio", "Simpson"
        write(*,21) exp(1.0)-exp(-1.0), trapecio, simpson

        20 format(3A15)         ! La letra A es de carácteres alfanúmericos.
        21 format(3F15.8)

end program cuadratura
