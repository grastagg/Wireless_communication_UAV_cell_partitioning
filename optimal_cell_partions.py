import numpy as np
import params as P


def average_data_service_at_xy(x,y,uav_index):
    #this function finds the average data service at a given location x,y from the given UAV position
    #we wish to maximize the average data service
    lambda_x = P.B/(P.N*P.alpha[uav_index])

def elevation_angle(x,y,uav_index):
    #this function finds the elevation angle of the given UAV at a given location x,y
    h = P.h[uav_index]
    d = np.sqrt((x-P.x[uav_index])**2 + (y-P.y[uav_index])**2 + h**2)
    theta = np.arcsin(h/d)
    return theta
    

def recieved_power(x,y,uav_index):
    P_i = P.P[uav_index]
    K_o = ((4*np.pi*P.f_c)/(3e9))**2
    theta_i = elevation_angle(x,y,uav_index)
    di_squared = (x-P.x[uav_index])**2 + (y-P.y[uav_index])**2 + P.h[uav_index]**2
    P_los = P.b1*(180/np.pi * theta_i - 15)**P.b2
    P_nlos = 1 - P_los
    
    den = K_o*di_squared*(P_los * P.mu_los+P_nlos*P.mu_nlos)
    
def recieved_interference(x,y,uav_index):
    interference = 0
    for j in range(P.num_UAVs):
        if j != uav_index:
            interference += recieved_power(x,y,j)
    return P.beta_interference_factor*interference


def gamma_i(x,y,uav_index):
    