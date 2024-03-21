import numpy as np
import params as P


def average_data_service_at_xy(x,y,uav_index):
    #this function finds the average data service at a given location x,y from the given UAV position
    #we wish to maximize the average data service
    lambda_x = P.B/(P.N*P.alpha[uav_index])
    