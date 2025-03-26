import numpy as np
import matplotlib.pyplot as plt

x_0 = 1.0e11 
y_0 = 0.0     
vx_0 = 3.0e4    
vy_0 = 3.0e4  
Mc = 1.989e30  
m = 5.972e24   
dt = 86400     
N = 1000      
U = [0,1,2]

def sectorial_velocity(R, V):
    return np.cross(R, V) * 1/2 

def force(R, fyc):
    a = 1
    n = 10
    r_norm = np.linalg.norm(R)  
    if fyc == 0:
        return -a * R / r_norm  
    if fyc == 1:
        return -1 / r_norm**3 * R / r_norm  
    if fyc == 2:
        return -n / r_norm**(n+1) * R / r_norm   



def movement2D(x_0, y_0, vx_0, vy_0, dt, N, force):
    fyc = 1
    r = np.zeros((N, 3))  
    t = np.zeros(N)       

    r[0] = [x_0, y_0, 0]  
    v = np.array([vx_0, vy_0, 0])  
    t[0] = 0  

    for i in range(1, N):
        r_norm = np.linalg.norm(r[i - 1])
        a = force(r[i-1],fyc) 

        r[i] = r[i - 1] + v * dt + 0.5 * a * dt**2

        r_norm_new = np.linalg.norm(r[i])
        a_new = force(r[i], fyc)

        v += 0.5 * (a + a_new) * dt

        t[i] = t[i - 1] + dt

    return r, t

r, t = movement2D(x_0, y_0, vx_0, vy_0, dt, N, force)

print (np.array([sectorial_velocity(r[i], np.array([vx_0, vy_0, 0])) for i in range(N)]))
angular_moments = np.array([sectorial_velocity(r[i], np.array([vx_0, vy_0, 0])) for i in range(N)])
print("первые 10 значений секториальной скорости:")
angular_moments_z = angular_moments[:, 2]
print(angular_moments_z[:10])

plt.figure(figsize=(10, 6))
plt.plot(t, angular_moments[:, 2], label="Секториальная скорость вокруг оси z", color='red')
plt.xlabel('Время (с)')
plt.ylabel('Секториальная скорость (м·кг/с)')
plt.title('Закон сохранения секториальной скорости')
plt.grid(True)
plt.legend()
plt.show()