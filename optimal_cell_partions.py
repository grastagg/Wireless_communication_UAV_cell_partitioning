import numpy as np
import params as P


def average_data_service_at_xy(x,y,uav_index):
    #this function finds the average data service at a given location x,y from the given UAV position
    #we wish to maximize the average data service
    lambda_x = P.B/(P.N*P.alpha[uav_index])

def elevation_angle(x,y,uav_index):
    #this function finds the elevation angle of the given UAV at a given location x,y
    h = P.h_UAV[uav_index]
    d = np.sqrt((x-P.x_locations[uav_index])**2 + (y-P.y_locations[uav_index])**2 + h**2)
    theta = np.arcsin(h/d)
    return theta
    

def recieved_power(x,y,uav_index):
    P_i = P.P[uav_index]
    K_o = ((4*np.pi*P.f_c)/(3e9))**2
    theta_i = elevation_angle(x,y,uav_index)
    di_squared = (x-P.x_locations[uav_index])**2 + (y-P.y_locations[uav_index])**2 + P.h_UAV[uav_index]**2
    P_los = P.b1*(180/np.pi * theta_i - 15)**P.b2
    P_nlos = 1 - P_los
    
    den = K_o*di_squared*(P_los * P.mu_los+P_nlos*P.mu_nlos)
    num = P_i/(P.N/P.num_UAVs)
    return num/den

def gamma(x,y,uav_index):
    return recieved_power(x,y,uav_index)/(recieved_interference(x,y,uav_index) + 10**(P.N_O/10)/1000)
    
def recieved_interference(x,y,uav_index):
    interference = 0
    for j in range(P.num_UAVs):
        if j != uav_index:
            interference += recieved_power(x,y,j)
    return P.beta_interference_factor*interference

def effective_data_tranmission_time(uav_index):
    return P.hover_time[uav_index] - 10*(P.N/P.num_UAVs)**2

def lambda_i(uav_index):
    sum = 0
    for i in range(uav_index):
        sum += effective_data_tranmission_time(i)
    return P.B/(P.N*P.alpha[uav_index])*sum



lambda_i_precomputed = np.zeros(P.num_UAVs)
for i in range(P.num_UAVs):
    lambda_i_precomputed[i] = lambda_i(i)

def J(x,y,uav_index):
    return -lambda_i_precomputed[uav_index]*np.log2(1+gamma(x,y,uav_index))

def omega_i(uav_index):
    sum = 0
    for i in range(P.num_UAVs):
        sum += P.alpha[i] * effective_data_tranmission_time(i) 
    return P.alpha[uav_index]*effective_data_tranmission_time(uav_index)/sum


        
precomputed_omega_i = np.zeros(P.num_UAVs)
for i in range(P.num_UAVs):
    precomputed_omega_i[i] = omega_i(i)


def integrate_d_i(psi):
    #d_i is defined in equation 27
    int_vals = np.zeros(P.num_UAVs)
    for x in P.x_int:
        for y in P.y_int:
            index = -1
            cost_val = np.inf
            for i in range(P.num_UAVs):
                if J(x,y,i) - psi[i] < cost_val:
                    index = i
                    cost_val = J(x,y,i) - psi[i]
            int_vals[index] += P.f(x,y)*P.dx*P.dy
    return int_vals
                


def compute_grad_f(psi):
    #defined in equation 27
    grad = np.zeros(P.num_UAVs)
    integral_over_d_i = integrate_d_i(psi)
    for i in range(P.num_UAVs):
        grad[i] = precomputed_omega_i[i] - integral_over_d_i[i] 
    return grad

    
    
    

    
def find_optimal_partiions():
    psi = np.ones(P.num_UAVs)
    # 
    grad_f = compute_grad_f(psi)
    print(grad_f)


if __name__ == '__main__':
    find_optimal_partiions()