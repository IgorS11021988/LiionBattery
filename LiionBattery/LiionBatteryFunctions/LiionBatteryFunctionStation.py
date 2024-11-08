import numpy as np

import LiionBattery.LiionBatteryFunctions.LiionBatteryFunctionsStation as lbfs
from MathProtEnergyProc import NonEqSystemQ

#Функция состояния для литий-ионного аккумулятора
def LiionBatteryStateFunction(stateCoordinates,
                              reducedTemp,
                              systemParameters):
    #получаем электрические заряды
    [qbinp,#Электрический заряд положительного двойного слоя
     qm,#Электрический заряд мембраны
     qbinn,#Электрический заряд отрицательного двойного слоя
     q#Перенесенный через внешнюю цепь электрический заряд
     ] = stateCoordinates
    
    #Получаем температуру
    [TInAkk,#Температура содержимого аккумулятора
     TBAkk#Температура корпуса аккумулятора
     ] = reducedTemp
    
    #Получаем параметры
    [I,#Ток во внешней цепи
     Tokr,#Температура окружающей среды
     EbinpC,#ЭДС положительного двойного слоя в заряженном состоянии
     EbinnC,#ЭДС отрицательного двойного слоя в заряженном состоянии
     EbinpD,#ЭДС положительного двойного слоя в разряженном состоянии
     EbinnD,#ЭДС отрицательного двойного слоя в разряженном состоянии
     Cbin0p,#Емкость положительного двойного слоя
     Cm,#Емкость мембраны
     Cbin0n,#Емкость отрицательного двойного слоя
     Rbin0p,#Сопротивление положительного двойного слоя
     Rm0,#Сопротивление мембраны
     Rbin0n,#Сопротивление отрицательного двойного слоя
     KInAkk,#Коэффициент теплопередачи содержимого литий-ионного аккумулятора
     CInAkk,#Теплоемкость содержимого литий-ионного аккумулятора
     KBAkk,#Коэффициент теплоотдачи корпуса литий-ионного аккумулятора
     CBAkk,#Теплоемкость корпуса литий-ионного аккумулятора
     Cnom,#Номинальная емкость литий-ионного аккумулятора
     rLiEpE,#Приведенное зарядовое число положительного электрода (по ЭДС)
     rLiEnE,#Приведенное зарядовое число отрицательного электрода (по ЭДС)
     alphaRIp,#Коэффициент сопротивления по току положительного двойного слоя
     alphaRIn,#Коэффициент сопротивления по току отрицательного двойного слоя
     alphaRQp,#Коэффициент сопротивления по перенесенному через внешнюю цепь заряду положительного двойного слоя
     alphaRQn,#Коэффициент сопротивления по перенесенному через внешнюю цепь заряду отрицательного двойного слоя
     nRQp,#Степенной коэффициент сопротивления по перенесенному через внешнюю цепь заряду положительного двойного слоя
     nRQn,#Степенной коэффициент сопротивления по перенесенному через внешнюю цепь заряду отрицательного двойного слоя
     alphaRTp,#Экспоненциальный коэффициент сопротивления по температуре положительного электрода
     alphaRTm,#Экспоненциальный коэффициент сопротивления по температуре мембраны
     alphaRTn,#Экспоненциальный коэффициент сопротивления по температуре отрицательного электрода
     bRTp,#Граничная температура по сопротивлению положительного электрода
     bRTm,#Граничная температура по сопротивлению мембраны
     bRTn,#Граничная температура по сопротивлению отрицательного электрода
     rCRTp,#Постоянный коэффициент температурной зависимости положительного электрода
     rCRTm,#Постоянный коэффициент температурной зависимости мембраны
     rCRTn,#Постоянный коэффициент температурной зависимости отрицательного электрода
     alphaCQp,#Зарядовый коэффициент емкости положительного электрода
     alphaCQn,#Зарядовый коэффициент емкости отрицательного электрода
     
     #Получаем довесочные коэффициенты
     betaRI2p,
     betaRI2n,
     betaRI3p,
     betaRI3n,
     betaRQ2p,
     betaRQ2n,
     betaRQ3p,
     betaRQ3n,
     betaRT2p,
     betaRT2m,
     betaRT2n,
     betaRT3p,
     betaRT3m,
     betaRT3n,
     betaCQ2p,
     betaCQ2n,
     betaCQ3p,
     betaCQ3n,
     betaEQ2p,
     betaEQ2n,
     betaEQ3p,
     betaEQ3n,
     
     Rkl#Сопротивление клемм
     ] = systemParameters
    
    #Матрица баланса
    balanceMatrix = np.array([])
    
    #Внешние потоки зарядов
    stateCoordinatesStreams = np.array([-I, -I, -I, I], dtype=np.double)
    
    #Внешние потоки теплоты
    heatStreambEnPow = Rkl*np.power(I, 2)#Внешний поток теплоты на корпус аккумулятора
    heatEnergyPowersStreams = np.array([heatStreambEnPow], dtype=np.double)
    
    #Выводим температуры
    energyPowerTemperatures = np.array([TInAkk,TBAkk,Tokr], dtype=np.double)
    
    #Рассчитываем приведенные (в единицах зарядовой емкости) числа молей интеркалированных в электроды ионов лития
    nuLip = qbinp + q#Положительный электрод
    nuLin = qbinn + q#Отрицательный электрод
    
    #Определяем ЭДС двойных слоев
    (Ebinp,Ebinn) = lbfs.funEbin(EbinpC,EbinnC,EbinpD,EbinnD,nuLip,
                                 nuLin,rLiEpE,rLiEnE,Cnom,
                                 betaEQ2p,betaEQ2n,betaEQ3p,betaEQ3n)
    
    #Определяем емкости двойных слоев
    (Cbinp,Cbinn) = lbfs.funCbin(qbinp,qbinn,alphaCQp,alphaCQn,Cbin0p,Cbin0n,
                                 betaCQ2p,betaCQ2n,betaCQ3p,betaCQ3n)
    
    #Определяем падения напряжения на двойных слоях
    dissUbinp = Ebinp - qbinp / Cbinp#Положительный двойной слой
    dissUbinn = Ebinn - qbinn / Cbinn#Отрицательный двойной слой
    
    #Потенциалы взаимодействия энергетических степеней свободы
    potentialInter = np.array([dissUbinp,-qm / Cm,dissUbinn,Ebinp + Ebinn], dtype=np.double)
    
    #Потенциалы взаимодействия между энергетическими степенями свободы
    potentialInterBet = np.array([])
    
    #Доли распределения некомпенсированной теплоты
    beta = np.array([])
    
    #Определяем сопротивления двойных слоев и мембраны
    (Rbinp,Rbinn,Rm) = lbfs.funRbin(alphaRIp,alphaRIn,dissUbinp,dissUbinn,
                                    nRQp,nRQn,alphaRQp,alphaRQn,nuLip,nuLin,
                                    alphaRTp,alphaRTn,bRTp,bRTn,rCRTp,rCRTn,
                                    alphaRTm,bRTm,rCRTm,TInAkk,Cnom,Rbin0p,
                                    Rbin0n,Rm0,
                                    betaRI2p,betaRI2n,betaRI3p,betaRI3n,
                                    betaRQ2p,betaRQ2n,betaRQ3p,betaRQ3n,
                                    betaRT2p,betaRT2m,betaRT2n,betaRT3p,
                                    betaRT3m,betaRT3n)
    
    #Главный блок кинетической матрицы по процессам
    kineticMatrixPCPC = np.array([1 / Rbinp, 1 / Rm, 1 / Rbinn], dtype=np.double) * TInAkk / NonEqSystemQ.GetTbase()
    
    #Перекрестные блоки кинетической матрицы по процессам
    kineticMatrixPCHeat = np.array([])
    kineticMatrixHeatPC = np.array([])
    
    #Главный блок кинетической матрицы по теплообмену
    kineticMatrixHeatHeat = np.array([KInAkk*TInAkk*TBAkk, KBAkk*TBAkk*Tokr], dtype=np.double) / NonEqSystemQ.GetTbase()
    
    #Обратная теплоемкость литий-ионного аккумулятора
    invHeatCapacityMatrixCf = np.array([1 / CInAkk, 1 / CBAkk], dtype=np.double)
    
    #Приведенные тепловые эффекты литий-ионного аккумулятора
    heatEffectMatrixCf = potentialInter / CInAkk
    
    #Выводим результат
    return (balanceMatrix,
            stateCoordinatesStreams,
            heatEnergyPowersStreams,
            energyPowerTemperatures,
            potentialInter,
            potentialInterBet,
            beta,kineticMatrixPCPC,
            kineticMatrixPCHeat,
            kineticMatrixHeatPC,
            kineticMatrixHeatHeat,
            invHeatCapacityMatrixCf,
            heatEffectMatrixCf)
