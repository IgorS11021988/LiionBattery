import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

from LiionBattery.LiionBatteryDynamics.ValuesGraphics import OneTimeValueGraphic, TimesValuesGraphics

#Функция расчета динамики
def LiionBatteryInputArrayCreate(Pars,#Параметры
                                 
                                 integrateAttributes#Аттрибуты интегрирования
                                 ):#Формирование массивов входных параметров
    #Корректируем аттрибуты токов
    Pars[["I",#Ток во внешней цепи
          "IStepDischarge",#Ток разряда после ступеньки, Cnom
          "AI"#Амплитуда колебаний тока во внешней цепи
          ]] *= Pars[["Cnom"]].to_numpy()
    integrateAttributes[["bI0DCh",#Граница разрядного нулевого тока
                         "bI0Ch"#Граница зарядного нулевого тока
                         ]] *= Pars[["Cnom"]].to_numpy()

    #Корректируем аттрибуты качающейся частоты
    Pars["AosI"] *= Pars["fI"] / 100#Амплитуда качающейся частоты

    #Относительная емкость
    rCnom = Pars["Cnom"] / 3.003

    #Корректируем частотные характеристики
    Pars[["AosI","fI","fosI"]] *= 2*np.pi

    #Емкости двойных слоев и мембраны
    Pars["Cbin0p"] *= rCnom#Емкость положительного двойного слоя, Ф
    Pars["Cm"] *= rCnom#Емкость мембраны, Ф
    Pars["Cbin0n"] *= rCnom#Емкость отрицательного двойного слоя, Ф

    #Емкости двойных слоев и мембраны
    Pars["Rbin0p"] /= rCnom#Сопротивление положительного двойного слоя, Ом
    Pars["Rm0"] /= rCnom#Сопротивление мембраны, Ом
    Pars["Rbin0n"] /= rCnom#Сопротивление отрицательного двойного слоя, Ом

    #Тепловые свойства литий-ионных аккумуляторов
    Pars["KInAkk"] *= rCnom#Коэффициент теплопередачи содержимого литий-ионного элемента, Вт/К
    Pars["CInAkk"] *= rCnom#Теплоемкость содержимого литий-ионного элемента, Дж/К
    Pars["KBAkk"] *= rCnom#Коэффициент теплопередачи корпуса литий-ионного элемента, Вт/К
    Pars["CBAkk"] *= rCnom#Теплоемкость корпуса литий-ионного элемента, Дж/К

    #Сопротивление клемм, Ом
    Pars["Rkl"] /= rCnom

    #Переводим номинальную емкость в кулоны
    Pars["Cnom"] *= 3600

    #Начальное состояние
    Pars["qbinp0"] *= Pars["EbinpC"]*Pars["Cbin0p"]#Заряд на положительном двойном слое, Кл
    Pars["qm0"] *= rCnom#Заряд мембраны, Кл
    Pars["qbinn0"] *= Pars["EbinnC"]*Pars["Cbin0n"]#Заряд на отрицательном двойном слое, Кл
    Pars["q0"] *= Pars["Cnom"]#Перенесенный через внешнюю цепь заряд, Кл
    Pars["TInAkk0"] += Pars["Tokr"]#Начальная температура сождержимого литий-ионного элемента, град С
    Pars["TBAkk0"] += Pars["Tokr"]#Начальная температура корпуса литий-ионного элемента, град С

    #Корректируем температуры
    Pars[["TInAkk0","TBAkk0","Tokr","bRTp","bRTm","bRTn"]] += 273.15

    #Время интегрирования
    Tints = np.array(integrateAttributes["TintI0"], dtype=np.double)#Времена интегрирования при нулевых токах
    bIZeros = np.logical_not(np.logical_and(Pars["I"] > -integrateAttributes["bI0Ch"],Pars["I"] < integrateAttributes["bI0DCh"]))
    if np.any(bIZeros):
        #Рассчитываем время интегрирования для ненулевого тока        
        sI = np.sign(Pars["I"][bIZeros])#Знак тока
        sc = (np.abs(sI) + sI) / 2#При заряде 0, при разряде 1
        sd = (np.abs(sI) - sI) / 2#При заряде 1, при разряде 0
        IZeros = 1.011*np.abs(Pars["I"][bIZeros])#Ненулевой ток (по модулю)
        Tints[bIZeros] = ((Pars["Cnom"][bIZeros] - Pars["q0"][bIZeros])*sc + Pars["q0"][bIZeros]*sd*integrateAttributes["cTimeCharge"]) / IZeros#Время интегрирования

    #Время начала конца заряда
    Pars["cTIStep"] = np.array(Pars["cTIStep"], dtype=np.double).reshape(-1,)#Приводим коэффициент времени ступенчатого перехода к массиву
    Pars["cTIStep"] *= Tints#Время начала конца разряда
    
    #Корректируем время разряда
    if np.any(bIZeros):
        bIZeros[bIZeros] = (sc > 0)#Режимы разряда
        if np.any(bIZeros):
            #Оставшееся время
            tost = Tints[bIZeros] - Pars["cTIStep"][bIZeros]
            bIZeros[bIZeros] = (tost > 0)#Индексы ненулевого оставшегося времени
            
            if np.any(bIZeros):
                #Оставшийся заряд
                cqost = Pars["Cnom"][bIZeros] - Pars["q0"][bIZeros] - Pars["I"][bIZeros]*Pars["cTIStep"][bIZeros]
                
                #Времена интегрирования
                IStepZeros = 1.011*np.abs(Pars["IStepDischarge"][bIZeros])#Ненулевой ток (по модулю)
                Tints[bIZeros] = Pars["cTIStep"][bIZeros] + cqost / IStepZeros#Обновленные времена интегрирования

    #Корректируем форму времени интегрирования
    if Tints.shape[0] > 1:
        Tints.shape = (-1,1)

    #Постоянная времени конца заряда
    Pars["cTauEndCharge"] *= Tints.reshape(-1,)

    #Массив параметров
    systemParameters = Pars[["I",#Ток во внешней цепи, А
                             "IStepDischarge",#Ток разряда после ступеньки, А
                             "fI",#Частота колебаний внешнего тока, Гц
                             "AI",#Амплитуда колебаний тока во внешней цепи, А
                             "fosI",#Частота качания частоты колебаний внешнего тока, Гц
                             "AosI",#Амплитуда качания частоты колебаний тока во внешней цепи, Гц
                             "cTIStep",#Время начала конца заряда
                             "cTauEndCharge",#Постоянная времени конца заряда
                             "alphaEndCharge",#Экспоненциальный коэффициент постоянной времени разряда
                             "powEndCharge",#Степень постоянной времени разряда
                             "Tokr",#Температура окружающей среды, град С
                             "EbinpC",#ЭДС положительного двойного слоя в заряженном состоянии, В
                             "EbinnC",#ЭДС отрицательного двойного слоя в заряженном состоянии, В
                             "EbinpD",#ЭДС положительного двойного слоя в разряженном состоянии, В
                             "EbinnD",#ЭДС отрицательного двойного слоя в разряженном состоянии, В
                             "Cbin0p",#Емкость положительного двойного слоя, Ф
                             "Cm",#Емкость мембраны, Ф
                             "Cbin0n",#Емкость отрицательного двойного слоя, Ф
                             "Rbin0p",#Сопротивление положительного двойного слоя, Ом
                             "Rm0",#Сопротивление мембраны, Ом
                             "Rbin0n",#Сопротивление отрицательного двойного слоя, Ом
                             "KInAkk",#Коэффициент теплопередачи литий-ионного элемента, Вт/К
                             "CInAkk",#Теплоемкость литий-ионного элемента, Дж/К
                             "KBAkk",#Коэффициент теплоотдачи корпуса литий-ионного элемента, Вт/К
                             "CBAkk",#Теплоемкость корпуса литий-ионного элемента, Дж/К
                             "Cnom",#Номинальная емкость элемента, А*ч
                             "rLiEpE",#Приведенное зарядовое число положительного электрода
                             "rLiEnE",#Приведенное зарядовое число отрицательного электрода
                             "alphaRIp",#Коэффициент сопротивления по току положительного двойного слоя, 1/В
                             "alphaRIn",#Коэффициент сопротивления по току отрицательного двойного слоя, 1/В
                             "alphaRQp",#Коэффициент сопротивления по перенесенному через внешнюю цепь заряду положительного двойного слоя
                             "alphaRQn",#Коэффициент сопротивления по перенесенному через внешнюю цепь заряду отрицательного двойного слоя
                             "nRQp",#Степенной коэффициент сопротивления по перенесенному через внешнюю цепь заряду положительного двойного слоя
                             "nRQn",#Степенной коэффициент сопротивления по перенесенному через внешнюю цепь заряду отрицательного двойного слоя
                             "alphaRTp",#Экспоненциальный коэффициент сопротивления по температуре положительного электрода, 1/К
                             "alphaRTm",#Экспоненциальный коэффициент сопротивления по температуре мембраны, 1/К
                             "alphaRTn",#Экспоненциальный коэффициент сопротивления по температуре отрицательного электрода, 1/К
                             "bRTp",#Граничная температура по сопротивлению положительного электрода, град С
                             "bRTm",#Граничная температура по сопротивлению мембраны, град С
                             "bRTn",#Граничная температура по сопротивлению отрицательного электрода, град С
                             "rCRTp",#Постоянный коэффициент температурной зависимости положительного электрода
                             "rCRTm",#Постоянный коэффициент температурной зависимости мембраны
                             "rCRTn",#Постоянный коэффициент температурной зависимости отрицательного электрода
                             "alphaCQp",#Зарядовый коэффициент емкости положительного электрода, 1/Кл
                             "alphaCQn",#Зарядовый коэффициент емкости отрицательного электрода, 1/Кл
                       
                             "betaRI2p",
                             "betaRI2n",
                             "betaRI3p",
                             "betaRI3n",
                             "betaRQ2p",
                             "betaRQ2n",
                             "betaRQ3p",
                             "betaRQ3n",
                             "betaRT2p",
                             "betaRT2m",
                             "betaRT2n",
                             "betaRT3p",
                             "betaRT3m",
                             "betaRT3n",
                             "betaCQ2p",
                             "betaCQ2n",
                             "betaCQ3p",
                             "betaCQ3n",
                             "betaEQ2p",
                             "betaEQ2n",
                             "betaEQ3p",
                             "betaEQ3n",
                             
                             "Rkl"
                             ]].to_numpy()

    #Массив начальных состояний
    stateCoordinates0 = Pars[["qbinp0","qm0","qbinn0","q0"]].to_numpy()
    reducedTemp0 = Pars[["TInAkk0","TBAkk0"]].to_numpy()

    #Моменты времени
    NPoints = np.array(integrateAttributes["NPoints"], dtype=np.int32)#Числа точек интегрирования
    tBegins = np.zeros_like(Tints).reshape(-1,)#Начальные моменты времени
    ts = []
    for ind in range(len(NPoints)):
        ts.append(np.linspace(tBegins[ind], Tints.reshape(-1,)[ind], NPoints[ind]))

    #Возвращаем исходные данные динамики системы
    return (Tints,
            stateCoordinates0,
            reducedTemp0,
            systemParameters,
            ts)
def LiionBatteryOutputValues(dyns,fileName,
                             sep,dec,
                             graphicFileDir,
                             graphicIndex,
                             plotGraphics = False#Необходимость построения графиков
                             ):
    #Получаем величины из кортежа
    (t, Ukl, Ubinp, Ubinn, Um, 
     TInAkk, TBAkk, q, Icur, Tokr) = dyns
    
    #Сохраняем динамику в файл
    DynamicDatas = pd.DataFrame({"Time": t.reshape(-1,),
                                 "Ukl": Ukl.reshape(-1,),
                                 "Ubinp": Ubinp.reshape(-1,),
                                 "Ubinn": Ubinn.reshape(-1,),
                                 "Um": Um.reshape(-1,),
                                 "TInAkk": TInAkk.reshape(-1,),
                                 "TBAkk": TBAkk.reshape(-1,),
                                 "q": q.reshape(-1,),
                                 "Icur": Icur.reshape(-1,),
                                 "Tokr": Tokr.reshape(-1,)
                                 })#Структура сохраняемых данных
    DynamicDatas.to_csv(fileName,sep=sep,decimal=dec,index=False)#Сохраняем в csv файл
    
    #Рисуем при необходимости график
    if plotGraphics:
        #Строим графики
        OneTimeValueGraphic(t,#Моменты времени
                            Ukl,#Величины в моменты времени
                            "Напряжение на клеммах",#Имя полотна
                            "Напряжение, В"#Имя оси
                            )#График напряжения на клеммах
        plt.savefig(graphicFileDir + "/Напряжение на клеммах " + str(graphicIndex) + ".jpg")#Сохраняем график
        TimesValuesGraphics(t,#Моменты времени
                            [TInAkk,TBAkk],#Список величин в моменты времени
                            ["Содержимое","Корпус"],#Список имен величин
                            "Температуры элемента",#Имя полотна
                            "Температура, град С",#Имя оси
                            )#Графики температуры содержимого и корпуса элемента
        plt.savefig(graphicFileDir + "/Температуры элемента " + str(graphicIndex) + ".jpg")#Сохраняем график
        TimesValuesGraphics(t,#Моменты времени
                            [Ubinn,Ubinp,Um],#Список величин в моменты времени
                            ["Отрицательный двойной слой","Положительный двойной слой","Мембрана"],#Список имен величин
                            "Напряжения в элементе",#Имя полотна
                            "Напряжение, В",#Имя оси
                            )#Графики напряжений двойных слоев и мембраны элемента
        plt.savefig(graphicFileDir + "/Напряжения в элементе " + str(graphicIndex) + ".jpg")#Сохраняем график
        OneTimeValueGraphic(t,#Моменты времени
                            Icur,#Величины в моменты времени
                            "Ток в цепи",#Имя полотна
                            "Ток, Cnom"#Имя оси
                            )#График тока во внешней цепи
        plt.savefig(graphicFileDir + "/Ток в цепи " + str(graphicIndex) + ".jpg")#Сохраняем график
def LiionBatteryParametersSave(dynamicIndex,#Индексы динамик
                               Pars,#Параметры
                               ParametersFileName,#Имя файла параметров
                               sep,#CSV разделитель
                               dec#Десятичный разделитель
                               ):#Сохранение параметров в файл
    dynamicIndex = pd.DataFrame({"dynamicIndex": dynamicIndex.reshape(-1,)})#Индекс динамики
    ParametersDatas = pd.concat([Pars,dynamicIndex],axis=1)#Параметры для 
    ParametersDatas.to_csv(ParametersFileName,sep=sep,decimal=dec,index=False)#Сохраняем в csv файл

#Формирование индекса из мультииндекса
def GetLiionAccumulatorMulIndex(MulIndex):#Получаем мультииндекс
    return MulIndex[["indexParameters","indexQCaps","indexTokrs","indexCurrents"]].to_numpy()
def LiionAccumulatorGenIndex(MulIndex,#Множественный индекс
                             nIndexes#Количества индексов
                             ):
    #Выделяем мультииндексы
    rMulIndex = GetLiionAccumulatorMulIndex(MulIndex)
    
    #Приводим количества индексов
    aNIndexes = np.array([1] + nIndexes, dtype=np.int32)#Приводим к массиву
    NIndexes = np.cumprod(aNIndexes[0:-1], dtype=np.int32).reshape(-1,1)#Получаем массив произведений
    
    #Удаляем индексы, большие соответсвующих чисел индексов
    bMoreRMulIndex = np.all((rMulIndex <= aNIndexes[1::]), axis=1)
    rMulIndex = rMulIndex[bMoreRMulIndex]
    
    #Вычисляем индекс
    return np.dot((rMulIndex - 1),NIndexes).reshape(-1,)

#Выделение ненулевых индексов
def SelectLiionAccumulatorNonZeroIndexes(MulIndex):
    #Получаем индексы нулевых и ненулевых индексов
    bNonZeroMultiIndex = (np.prod(GetLiionAccumulatorMulIndex(MulIndex), axis=1) > 0)
    
    #Факт наличия ненулевого мультииндекса
    isNonZeroMultiIndex = np.any(bNonZeroMultiIndex)
    
    #Выделяем ненулевые мультииндексы
    NonZeroMultiIndexes = MulIndex.loc[bNonZeroMultiIndex]
    
    #Выводим результат
    return (bNonZeroMultiIndex, isNonZeroMultiIndex, NonZeroMultiIndexes)

#Аттрибуты интегрирования
def IntegrateAttributes(integrateAttributes,#Аттрибуты интегрирования
                        arrNAllAccAttrs,#Массив чисел всех аттрибутов аккумулятора
                        nDyns#Число динамик
                        ):
    #Выделяем прочие аттрибуты
    (bIntegrateAttributesOther,
     isIntegrateAttributesOther,
     integrateAttributesOther) = SelectLiionAccumulatorNonZeroIndexes(integrateAttributes)
    
    #Выделяем базовые аттрибуты интегрирования
    integrateAttributesBase = integrateAttributes.loc[np.logical_not(bIntegrateAttributesOther)].loc[0]
    
    #Формируем аттрибуты интегрирования
    integrateAttributes = pd.DataFrame(np.full((nDyns,len(list(integrateAttributes))), integrateAttributesBase.to_numpy()),
                                       columns=list(integrateAttributes))#Приравниваем базовые аттрибуты интегрирования
    if isIntegrateAttributesOther:
        indIntegrateAttributesOther = LiionAccumulatorGenIndex(integrateAttributesOther, arrNAllAccAttrs)#Инденсы прочих аттрибутов интегрирования
        integrateAttributes.loc[indIntegrateAttributesOther] = integrateAttributesOther.to_numpy()#Приравниваем прочие аттрибуты интегрирования
    
    #Возвращаем аттрибуты интегрирования
    return integrateAttributes

#Индексы для построения графиков
def IndexesGraphics(indexesGraphics,#Аттрибуты интегрирования
                    arrNAllAccAttrs,#Массив чисел всех аттрибутов аккумулятора
                    ):
    #Выделяем индексы для построения графиков
    (bBuildGraphicsIndexes,
     isBuildGraphicsIndexes,
     buildGraphicsIndexes) = SelectLiionAccumulatorNonZeroIndexes(indexesGraphics)
    if isBuildGraphicsIndexes:
        #Формируем индексы графиков для построения
        buildGraphicsIndexes = LiionAccumulatorGenIndex(buildGraphicsIndexes, arrNAllAccAttrs) + 1#Индексы графиков для построения
        
    #Выводим результат
    return (buildGraphicsIndexes,isBuildGraphicsIndexes)
