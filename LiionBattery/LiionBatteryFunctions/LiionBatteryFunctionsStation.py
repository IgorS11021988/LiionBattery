import numpy as np

#Вспомогательные функции
def funRI(alphaRI,dissU):#Мультипликативная корректировка по току через двойной слой
    #Проверяем падение напряжения
    _dissU = alphaRI*dissU
    if (np.fabs(_dissU) > 0.001):
        return 2*_dissU / (np.exp(_dissU) - np.exp(-_dissU))
    else:
        return 1
def funRQ(nRQ,Cnom,dissU,alphaRQ,nuLi):#Мультипликативная корректировка по заряду через внешнюю цепь
    #Приведенное число молей интеркалированных ионов лития
    nuLi = nuLi / (alphaRQ*Cnom)
    
    #Получаем коорректировочный коэффициент
    if dissU > 0:
        cf = 1 - nuLi
    elif dissU < 0:
        cf = nuLi
    else:
        cf = 1
        
    if (cf < 0):
        cf = 0
    cf = np.power(cf,nRQ)
    if (cf < 1e-51):
        return 1e+51
    else: 
        return 1 / cf
def funRT(alphaRT,bRT,TInAkk,rCRT):#Мультипликативная корректировка по температуре
    return np.exp(-alphaRT*(TInAkk - bRT)) + rCRT
def funCQbin(qbin,alphaCQ):#Емкость двойного слоя в зависимости от заряда
    return np.exp(-alphaCQ*np.abs(qbin))

#Функции для свойств веществ и процессов
def funNuEbin(nuLip,nuLin,rLiEpE,rLiEnE,Cnom,
              betaEQ2p,betaEQ2n,betaEQ3p,betaEQ3n):
    #Определяем коэффициент ЭДС, обусловленный переносом электрического заряда через внешнюю цепь
    EnuLip = nuLip / (rLiEpE*Cnom)#Положительный двойной слой
    EnuLin = nuLin / (rLiEnE*Cnom)#Отрицательный двойной слой
    
    #Добавляем довесочные члены к корректировкам сопротивления двойных слоев через перенесенный через внешнюю цепь электрический заряд
    EnuLip += betaEQ2p*np.power(EnuLip, 2) + betaEQ3p*np.power(EnuLip, 3)
    EnuLin += betaEQ2n*np.power(EnuLin, 2) + betaEQ3n*np.power(EnuLin, 3)
    
    #Выводим результат
    return (EnuLip,EnuLin)
def funEbin(EbinpC,EbinnC,EbinpD,EbinnD,nuLip,
            nuLin,rLiEpE,rLiEnE,Cnom,
            betaEQ2p,betaEQ2n,betaEQ3p,betaEQ3n):
    #Определяем коэффициент ЭДС, обусловленный переносом электрического заряда через внешнюю цепь
    (EnuLip,EnuLin) = funNuEbin(nuLip,nuLin,rLiEpE,rLiEnE,Cnom,
                                betaEQ2p,betaEQ2n,betaEQ3p,betaEQ3n)
    
    #Определяем ЭДС двойных слоев
    Ebinp = EbinpC - (EbinpC - EbinpD)*EnuLip#Положительный двойной слой
    Ebinn = EbinnC - (EbinnC - EbinnD)*EnuLin#Отрицательный двойной слой
    
    #Выводим результат
    return (Ebinp,Ebinn)
def funRbin(alphaRIp,alphaRIn,dissUbinp,dissUbinn,
            nRQp,nRQn,alphaRQp,alphaRQn,nuLip,nuLin,
            alphaRTp,alphaRTn,bRTp,bRTn,rCRTp,rCRTn,
            alphaRTm,bRTm,rCRTm,TInAkk,Cnom,Rbin0p,
            Rbin0n,Rm0,
            betaRI2p,betaRI2n,betaRI3p,betaRI3n,
            betaRQ2p,betaRQ2n,betaRQ3p,betaRQ3n,
            betaRT2p,betaRT2m,betaRT2n,betaRT3p,
            betaRT3m,betaRT3n):#Функция сопротивления
    #Определяем корректировку сопротивления двойных слоев через токи двойных слоев
    rIbinp = funRI(alphaRIp,dissUbinp)#Положительный двойной слой
    rIbinn = funRI(alphaRIn,dissUbinn)#Отрицательный двойной слой
    
    #Добавляем довесочные члены к корректировкам сопротивления двойных слоев через токи двойных слоев
    rIbinp += betaRI2p*np.power(rIbinp, 2) + betaRI3p*np.power(rIbinp, 3)
    rIbinn += betaRI2n*np.power(rIbinn, 2) + betaRI3n*np.power(rIbinn, 3)
    
    #Определяем корректировку сопротивления двойных слоев через перенесенный через внешнюю цепь электрический заряд
    rQbinp = funRQ(nRQp,Cnom,dissUbinp,alphaRQp,nuLip)#Положительный двойной слой
    rQbinn = funRQ(nRQn,Cnom,dissUbinn,alphaRQn,nuLin)#Отрицательный двойной слой
        
    #Добавляем довесочные члены к корректировкам сопротивления двойных слоев через перенесенный через внешнюю цепь электрический заряд
    rQbinp += betaRQ2p*np.power(rQbinp, 2) + betaRQ3p*np.power(rQbinp, 3)
    rQbinn += betaRQ2n*np.power(rQbinn, 2) + betaRQ3n*np.power(rQbinn, 3)
    
    #Определяем корректировку сопротивления двойных слоев через температуру
    rTbinp = funRT(alphaRTp,bRTp,TInAkk,rCRTp)
    rTbinn = funRT(alphaRTn,bRTn,TInAkk,rCRTn)
    
    #Добавляем довесочные члены к корректировкам сопротивления двойных слоев через температуру
    rTbinp += betaRT2p*np.power(rTbinp, 2) + betaRT3p*np.power(rTbinp, 3)
    rTbinn += betaRT2n*np.power(rTbinn, 2) + betaRT3n*np.power(rTbinn, 3)
    
    #Определяем сопротивления двойных слоев
    rbinp = rIbinp * rQbinp * rTbinp#Положительный двойной слой
    rbinn = rIbinn * rQbinn * rTbinn#Отрицательный двойной слой

    #Определяем сопротивление мембраны
    rm = funRT(-alphaRTm,bRTm,TInAkk,rCRTm)
    
    #Добавляем довесочные члены к корректировкам сопротивления мембраны через температуру
    rm += betaRT2m*np.power(rm, 2) + betaRT3m*np.power(rm, 3)
    
    #Выводим результат
    return (Rbin0p*rbinp,Rbin0n*rbinn,Rm0*rm)
def funCbin(qbinp,qbinn,alphaCQp,alphaCQn,Cbin0p,Cbin0n,
            betaCQ2p,betaCQ2n,betaCQ3p,betaCQ3n):#Функция емкостей двойных слоев
    #Определяем корректировочный коэффициент емкости двойного слоя
    rCbinQp = funCQbin(qbinp,alphaCQp)#Положительный двойной слой
    rCbinQn = funCQbin(qbinn,alphaCQn)#Отрицательный двойной слой
    
    #Учитываем довесочные слагаемые коэффициента емкости двойного слоя
    rCbinQp1 = rCbinQp - 1; rCbinQp += np.power(rCbinQp1,2) + np.power(rCbinQp1,3)#Положительный двойной слой
    rCbinQn1 = rCbinQn - 1; rCbinQn += np.power(rCbinQn1,2) + np.power(rCbinQn1,3)#Отрицательный двойной слой
    
    #Выводим результат
    return (Cbin0p*rCbinQp,Cbin0n*rCbinQn)
