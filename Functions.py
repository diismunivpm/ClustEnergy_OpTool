"""
ClustEnergy OpTool version 1.0.0
Released on February 29, 2024
@author: DIISM UNIVPM
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""
import numpy as np


"""
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 

Funzione che permette di scrivere il sistema per i vincoli lineari 
(1) Per edifici (optForm_build())
richiede: 
    samples -> numero timestep simulazione
    Adopt   -> matrice A del sistema discretizzato scritto nella forma SSM
    Xo      -> condizione iniziale del sistema (variabili di stato)
    Edopt   -> matrice E del sistema discretizzato che isola la variabile decisionale 
    Uopt    -> vettore di input 
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
"""

def optForm_build(samples, Adopt, Xo, Bdopt, Uopt, Edopt, p_vd):
    """The optForm_build() function allows to write the building space state system for the linear constraints of the optimization problem. 
        For the calculation the following are required:
            samples -> simulation timestep number
            Adopt -> matrix A of the discretized system written in SSM form
            Xo -> initial condition of the system (state variables)
            Edopt -> matrix E of the discretized system isolating the decision variable 
            Uopt -> input vector """
    
    Nvs       = Adopt.shape[0]
    Pvv       = p_vd   # variable position (system state) to be bound - Tair
    duration  = range(0,samples)
    
    M_opt_max = np.empty([samples,1])
    M         = np.empty([Nvs,samples])
    Z_opt     = np.zeros([samples,samples])
    
    for k in duration:
        if k == 0:
            m          = np.dot(Adopt, Xo) + np.dot(Bdopt, Uopt[:,k]).reshape(Nvs,1)
            Z_opt[0,0] = Edopt[Pvv,0]
        else: 
            Adis  = np.dot(np.linalg.matrix_power(Adopt,(k+1)), Xo)
            b_dis = np.empty([Nvs,k])
            
            for i in range(0,k):
                b_dis[:,i] = np.dot(np.linalg.matrix_power(Adopt,i),np.dot(Bdopt,Uopt[:,(k-i)]).reshape(Nvs,1)).reshape(Nvs,1)[:,0]
           
            for i in range(0,k+1):
                Z_opt[k,(k-i)] = np.dot(np.linalg.matrix_power(Adopt,i),Edopt)[Pvv,0]  
                        
            Bdis = np.sum(b_dis, axis = 1, keepdims = True)
            m = (Adis+Bdis).reshape(Nvs,1)
        M[:,k] = m[:,0]
    
    M_opt_max = M[Pvv,:].reshape(samples,1)
    return Z_opt, M_opt_max, M




