# -*- coding: utf-8 -*-
"""
Created on Sun Mar  1 09:25:45 2020

@author: lishuqi
"""

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif']=['SimHei'] #用来正常显示中文标签
plt.rcParams['axes.unicode_minus']=False #用来正常显示负号

"""
参数修改区
"""
steps = 1
seir= [9999, #S
         0,    #E
         1,    #I
         0]    #R

#基础SEIR的参数
para = [20,     # I每天接触的人  r
        0.03,   # S和E被感染的概率  beta
        0.1,    # 潜伏者变为感染者的概率 a
        0.1,    # 康复/治愈概率  gamma
        10000]  # 总人数

#潜伏期也具有传染性SEIRS_II
para_II = [20, # 感染者每天接触的人 r 
          20, # 潜伏者每天接触的人 r2 
          0.03, # 感染者传染的概率 beta
          0.03, # 潜伏者传染的概率 beta2 
          0.1, # 潜伏者变为感染者的概率 a 
          0.1, # 康复/治愈概率 gamma 
          10000]


def SEIR(seir,para,steps):
    S,E,I,R = seir
    r,beta,a,gamma,N = para #赋参数
    dS = -(r*beta*I*S)/N #表示微分方程
    dE =  (r*beta*I*S)/N - a*E
    dI = a*E - gamma*I
    dR = gamma*I
    return [S+dS*steps, E+dE*steps, I+dI*steps, R+dR*steps]

def SEIR_II(seir,para,steps): 
    S,E,I,R = seir 
    r,r2,beta,beta2,a,gamma,N = para 
    dS = -(r*beta*I*S)/N - (r2*beta2*E*S)/N 
    dE = (r*beta*I*S)/N - a*E + (r2*beta2*E*S)/N 
    dI = a*E - gamma*I 
    dR = gamma*I 
    return [S+dS*steps, E+dE*steps, I+dI*steps, R+dR*steps]

def plot_graph(np_res):#画图函数
    plt.figure()
    plt.plot(np_res[:,0])
    plt.plot(np_res[:,1])
    plt.plot(np_res[:,2])
    plt.plot(np_res[:,3])
    plt.title("SEIR")
    plt.xlabel("日期/天")
    plt.ylabel("人数/人)")
    plt.legend(['易感者','潜伏者','感染者','康复者'])
    plt.show()
    
def calculate(func,seir,para,intervene_N):  
    steps = 1
    t = np.arange(0,200,steps)
    res=[]
    for itm in t:
        if intervene_N!=0:
            if itm>intervene_N:
                para[0]=4
                para[1]=4
        seir=func(seir,para,steps)
        res.append(seir)
    return np.array(res)

result = calculate(SEIR_II,seir,para_II,0)
plot_graph(result)
