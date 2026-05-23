from main import IsingModel
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid
from scipy.interpolate import CubicSpline
import numpy as np
 #3) - (a)
def equilibrium_ener(N_meas, L, T):
  #Função que devolve a energia de equilíbrio
  Model = IsingModel(L,T)
  Model.iter_monte_carlo(N_meas)
  return (np.mean(Model.energies[200000:]))
#Lista de temperaturas e função a integrar 
T_list = np.arange(0.1, 10.1, 0.1)
F = lambda T: equilibrium_ener(1000000,128,T)/(T**2)
print('3) - (a)\n')
Value_list = []
#Número de pontos a utilizar para integrar
Slices = 100
x = np.linspace(0.1,10.1, Slices+1)
y = np.array([F(t) for t in x])
for i in T_list:
  #A energia livre de Helmholtz
  Value = -i*(trapezoid(y[int(Slices*(i-0.1)/10):], x[int(Slices*(i-0.1)/10):]) + np.log(2))*(128**2)
  Value_list.append(Value)
  print(f'T = {i} --> F = {Value}')
plt.plot(T_list, Value_list)
plt.xlabel("T")
plt.ylabel("F")
plt.show()
#3) - (b)
#Aproximação da função
Spline = CubicSpline(T_list, Value_list)
Not_found = True
l = 1
Zero = 0
lastfound = 0
Precision = 10
#Algoritmo de procura binária
while Not_found:
  if (Spline(Zero+10/(2**l))+4*(128**2))*(Spline(Zero)+4*(128**2)) > 0:
    Zero+=10/(2**l)
    if Zero >= 10:
      Zero = lastfound
      l += 1
  elif (Spline(Zero+10/(2**l))+4*(128**2))*(Spline(Zero)+4*(128**2)) < 0:
    if l >= Precision:
      Zero = Zero + 10/(2**(l+1))
      Not_found = False
      break
    lastfound = Zero
    l += 1
print(Zero)
#3) - (c)
second_derivative = Spline.derivative(2)
T_rep = np.linspace(0.1,10.1,1001)
Value_rep = np.array([second_derivative(T) for T in T_rep])
plt.plot(T_rep, Value_rep)
plt.show()
#A primeira derivada fica aproximadamente constante na temperatura crítica, isto justifica-se com a aproximação de F(T) ~ -T*N*ln(2) visto que a energia a parttir da temperatura crítica é aproximadamente zero o integral de Tc até infinito de E(T)/T**2 vai para zero
#Dado que a derivada de uma constante é zero, a segunda derivada é zero 