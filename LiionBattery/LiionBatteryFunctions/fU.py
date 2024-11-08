import numpy as np

import LiionBattery.LiionBatteryFunctions.fICharge as fic

#Функция условий протекания процессов
def fU(t,#Моменты времени
       systemParameters#Параметры системы  
       ):
    #Выделяем параметры динамики и отдельно свойства веществ и процессов
    (Icur,otherSystemParameters) = fic.fICharge(np.array([t], dtype=np.double),#Моменты времени
                                                systemParameters#Параметры системы  
                                                )
    
    #Выводим результат
    return np.hstack((Icur,otherSystemParameters))
    