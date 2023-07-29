import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import simpledialog

#estos son los parámetros que el usuario debería de introducir enj pantalla
def obtener_parametros():
    root = tk.Tk()
    root.withdraw()
    m = simpledialog.askfloat("Ingresar valor", "Ingrese el valor de 'm':")
    b = simpledialog.askfloat("Ingresar valor", "Ingrese el valor de 'b':")
    k = simpledialog.askfloat("Ingresar valor", "Ingrese el valor de 'k':")

    return m, b, k

#los valores de los parámetros necesarios
m, b, k = obtener_parametros()
print(m,b,k) #para verificar que se introdujo los parámetros correctos


# Función que representa el sistema de ecuaciones de primer orden
def differential_equations(y, t, m, b, k):
    v, u = y  # y[0] = v, y[1] = u
    dvdt = -b / m * v - k / m * u
    dudt = v
    return [dvdt, dudt]

# Parámetros del sistema (lo pueden introducir aquí también)
#m = 1  # masa (positiva)
#b = 1  # coeficiente de amortiguamiento (positivo)
#k = 1  # constante del resorte (positiva)

# Condiciones iniciales
v0 = 0.0  # velocidad inicial
u0 = 1.0  # posición inicial

# Tiempo
t = np.linspace(0, 40, 4000)  # Intervalo de tiempo de 0 a 100 con 10000 puntos

# Resolver la ecuación diferencial de forma numérica utilizando odeint
sol = odeint(differential_equations, [v0, u0], t, args=(m, b, k))

# Extraer la posición y velocidad de las soluciones
v = sol[:, 0]
u = sol[:, 1]


#la solución análitica del sistema está resuelto para diferentes valores de los parámetros m, b y k (directamente de la literatura)
def solt_analitica(m, b, k): #la solución analítica en función a los parámetros 
    if b**2-4*k*m > 0: #sistema sobreamortiguado (amortiguamiento fuerte o supercrítico)
        lambda1=-(b+np.sqrt(b**2-4*k*m))/(2*m)
        lambda2=(-b+np.sqrt(b**2-4*k*m))/(2*m)
        A1=(u0*lambda2-v0)/(lambda2-lambda1)
        A2=(v0-u0*lambda1)/(lambda2-lambda1)
        y=A1*np.exp(lambda1*t)+A2*np.exp(lambda2*t) #solución analítica
        coment='$b^2-4km > 0$, sist. sobreamort.'
    elif b**2-4*k*m == 0: #sistema con amortiguamiento crítico
        lambda1=-(b+np.sqrt(b**2-4*k*m))/(2*m)
        lambda2=(-b+np.sqrt(b**2-4*k*m))/(2*m)
        A1=u0
        A2=v0+u0*b/2*m
        y=A1*np.exp(lambda1*t)+A2*t*np.exp(lambda2*t) #solución analítica
        coment='$b^2-4km=0$, sist. con amort. crítico'
    else: #sistema con amortiguamiento débil o subcrítico, b**2-4*k*m<0
        omeg=np.sqrt(k/m - (b/2*m)**2) 
        A=u0/(np.cos(np.arctan((-b*u0/(2*m) - v0)/(u0*omeg))))
        y=A*np.exp(-b*t/(2*m))*np.cos(omeg*t+np.arctan((-b*u0/(2*m) - v0)/(u0*omeg))) #solución analítica
        coment='$b^2-4km<0$, sist. con amort. subcrítico'
    return y,coment

solan=solt_analitica(m,b,k) #solución analítica para cada valor dado de los parámetros
#solan[0] corresponde a la solución análitica y mientras que solan[1] es el comentario que se le agrega

#Graficas para comparar las soluciones analíticas y numéricas
plt.figure(1)
plt.plot(t, solan[0], '-.', label='Sol. Analítica')
plt.plot(t, u, label='Sol. Numérica')
plt.xlabel('Tiempo')
plt.ylabel('Posición (x)')
plt.legend(loc='best')
#plt.ylim(0)
plt.xlim(0)
plt.title('$md^2x/dt^2 + b dx/dt + kx = 0$, {}'.format(solan[1]))
#plt.grid(True)
#plt.show()

plt.figure(2)
plt.plot(u, v, label='Sol. Numérica')
plt.ylabel('Velocidad (V)')
plt.xlabel('Posición (x)')
plt.legend(loc='best')
plt.title('Espacio de fase') #para observar el punto de convergencia
plt.grid(True)
