import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d
from scipy.optimize import curve_fit #Use non-linear least squares to fit a function, f, to data.
from scipy.integrate import solve_ivp, quad
from scipy.linalg import eig # Solve an ordinary or generalized eigenvalue problem of a square matrix.


########################################################
def systemBS(x, yV, arg):
    """
    SISTEMA de EDO que describen una estrella de Bosones en el contexto de EKG
    
    L = k R - (CD[-a][PhiC]*CD[a][Phi]/2 + m_0^2 Phi PhiC+\lambda\phi^4)

    donde R: escalar de Ricci, k=c^4/(16 G pi), CD indica la derivada covariante y m la masa del campo.

    El sistema es adimensionalizado de tal forma que:
    rF = r/m, 
    PhiF = Mp/ Phi,
    wF = m w

    las cantidades con F son las físicas (no adimensionales)
    
    COMENTARIOS:
    El sistema es dividido en r=0, r>0, ver notebook de mathematica

    IMPLEMENTACION:

    Sistema:
    A' = f0
    B' = f1
    Phi' = f2
    Phi'' = f3

    Variables:
    [g, g', N, N', Phi, Phi'] -> [g0, g1, n0, n1, p0, p1]

    In:
    r  -> ri 
    yV -> [g(ri), n(ri), p(ri), p1(ri)]
    arg -> [w, aB]

    Out:
    [f0, f1, f2, f3, f4, f5]  ->  [g', N', Phi', Phi'', m', n']
    """

    g0, n0, p0, p1, m0, N0 = yV
    w, Lambda = arg

    if p0>80 or p0<-80:
        #sys.exit('El perfil se indeterminó')
        f0, f1, f2, f3, f4, f5 = 0, 0, 0, 0, 0, 0
    elif x > 0:
        f0 = -(-g0 + g0**2)/x +(g0**2)*x*(Lambda*(p0**4) + (p1**2)/g0 + (p0**2)*(1 + (w**2)/n0))
        f1 = ((-1 + g0)*n0)/x + g0*n0*x*(-Lambda*(p0**4) + (p1**2)/g0 + (p0**2)*(-1 + (w**2)/n0))
        f2 = p1
        g1, n1 = f0, f1
        f3 = g0*(2*Lambda*(p0**3) - p0*(-1 + (w**2)/n0))-p1*(-g1/(2*g0)+n1/(2*n0) + 2/x)
        f4 = ((x**2)/2)*((w**2/n0+1)*p0**2+Lambda*p0**4+(p1**2)/g0)
        f5 = (x**2)*(p0**2)*((g0/n0)**(1/2))*w
    else:
        f0 = 0
        f1 = 0
        f2 = p1
        f3 = (p0*(2-w**2/n0**2)+2*Lambda*(p0**3))/3
        f4 = 0
        f5 = 0
        
    return [f0, f1, f2, f3, f4, f5]

########################################################
# MI implementación
def Freq_solveNodos(in0, wmin, wmax, Lambda, nodos, rmin=0, rmax=50, met='RK45', Rtol=1e-06, Atol=1e-7,
               limw=1e-14, info=False):
    """
    Implementación de un algoritmo de shooting usando un método de bisección para encontrar el valor de la 
    frecuencia w0, dado una amplitud central p0.

    In:
    system -> sitema de ecuaciones: eds.systemBS
    in0 -> condiciones iniciales del problema: [g0, n0, p0, p1] recordar que [g, g', N, N', Phi, Phi'] -> [g0, g1, n0, n1, p0, p1]
    wmin, wmax -> rango de valores en el que se buscará la frecuencia (recomendado añadir un check al código que verifique que wmin<wmax)
    rmax, rmin -> intervalo que se discretizará, por defecto rmin=0, rmax=1e03 (notar que lo q se busca es que alcance el límite en w0)
    met -> metodología usada por solve_ivp, las opciones son: 'DOP853', 'LSODA', 'RK45', por defecto está este último
    Rtol, Atol=1e-10 -> representan la tolerancia relativa y absoluta al usar solve_ivp, por defecto son Rtol=1e-09, Atol=1e-10
    limw -> representa la diferencia limítrofe que aceptaremos: abs((wmax-wmin)/2)<=limw Por defecto es limw=1e-14 
    alphB -> Valor de alphaB, por defecto alphB=1
    nodos -> los nodos de la configuracion
    info -> Imprime información complementaria, por defecto info=False

    Out:
    w0, rTemp, nodosPosit  -> valor de la frecuencia encontrada, el radio máximo de la iteración y la posición de los nodos
    """
    
    print('Finding a profile with ', nodos, 'nodes')
    
    def Sig(r, yV, arg): return yV[2]
    def dSig(r, yV, arg): return yV[3]

    # establece las direcciones
    Sig.direction = 0  # como pueden ser varios nodos se pone 0. Notar que no hay acción que tomar, solo almacenamos 
    dSig.direction = 0
    

    while True:
        w0 = (wmax+wmin)/2
        arg = [w0, Lambda]

        sol = solve_ivp(systemBS, [rmin, rmax], in0, events=(Sig, dSig),
                         args=(arg,), method=met,  rtol=Rtol, atol=Atol)
        
        if sol.t_events[1].size == nodos+1 and sol.t_events[0].size == nodos:
            print('Found', w0)
            rTemp = rmax
            nodosPosit = sol.t_events[0]
            break
        elif sol.t_events[1].size > nodos+1:  # una vez por nodo
            if sol.t_events[0].size > nodos:  # dos veces por nodo
                wmax = w0
                rTemp = sol.t_events[0][-1]
            else:  # si pasa por cero más veces que 2*nodos se aumenta la w, sino se disminuye
                wmin = w0
                rTemp = sol.t_events[1][-1]
        elif sol.t_events[1].size <= nodos+1:
            if sol.t_events[0].size > nodos:  # dos veces por nodo
                wmax = w0
                rTemp = sol.t_events[0][-1]
            else:
                wmin = w0
                rTemp = sol.t_events[1][-1]

        # checking the lim freq.
        if abs((wmax-wmin)/2)<=limw:
            print('Maxima precisión alcanzada', w0, 'radio', rTemp)
            nodosPosit = sol.t_events[0]
            break
        
        if nodos==0:
            nodosPosit = None

    return w0, rTemp, nodosPosit 


########################################################

def ParticleNumber(r, sigma_0, A, B, omega):
    """

    """
    sigF = interp1d(r, sigma_0, kind='quadratic') 
    InVAf = 1/A
    

    Bf = interp1d(r,  B, kind='quadratic') 
    InvAf = interp1d(r, InVAf, kind='quadratic') 
    
    NInt = lambda r: (r**2*sigF(r)**2)*((Bf(r)*InvAf(r))**(1/2))*omega
 
        
    rmin = r[0]
    rfin = r[-1]
    
    N = quad(NInt,rmin, rfin)[0]
    
    return N
########################################################

# [f0, f1, f2, f3, f4]  ->  [g', N', Phi', Phi'', m']

def profiles(datos, in0, rmin, Rtol, Atol, info = False): 

        w = datos[0]
        # w, Lambda = arg
        arg = [w, datos[3]]
        sol2 = solve_ivp(systemBS, [rmin, datos[1]], in0, args=(arg,), method='RK45', rtol=Rtol, atol=Atol)
        

        
        N = ParticleNumber(sol2.t, sol2.y[2], sol2.y[0], sol2.y[1], w)
        
        if info: 
            plt.plot(sol2.t, sol2.y[4])
            plt.plot(sol2.t, sol2.y[5])
            plt.plot(sol2.t, sol2.y[2])
            plt.xlabel(r'r')
            plt.xlim(0,datos[1])
            plt.ylabel(r'$\sigma(r)$')
            plt.axhline(y=0, color='r', linestyle='--')
            plt.show()
        
        return sol2.t, sol2.y[0], sol2.y[1], sol2.y[2], sol2.y[3], sol2.y[4], sol2.y[5], N


########################################################
def progressbar(current_value, total_value, bar_lengh, progress_char): 
    """
    Barra de progreso
    """
    percentage = int((current_value/total_value)*100)                                                # Percent Completed Calculation 
    progress = int((bar_lengh * current_value ) / total_value)                                       # Progress Done Calculation 
    loadbar = "Progress: [{:{len}}]{}%".format(progress*progress_char,percentage, len = bar_lengh)    # Progress Bar String
    print(loadbar, end='\r')

########################################################




########################################################




########################################################




########################################################



########################################################



########################################################



########################################################


########################################################


########################################################

########################################################

########################################################
########################################################