import numpy as np
import pandas as pd
import json as js

from MathProtEnergyProc.DatasAugmentation import HStackMatrixRepeat
from MathProtEnergyProc.IndexedNames import IndexedNamesFromIndexes

from LiionBattery.LiionBatteryDynamics.LiionBatteryDynamics import LiionBatteryDynamicFunctionSaveVec
from LiionBattery.LiionBatteryDynamics.LiionBatteryValues import IntegrateAttributes, IndexesGraphics, LiionBatteryParametersSave

#Функция запуска проекта
def ExecuteProject(ProjectFileName,modelingDynamicIndicate):
    #Открываем файл проекта
    with open(ProjectFileName, 'r') as ProjFileName:
        ProjectsAttributes = js.load(ProjFileName)

    #Номинальная емкость аккумулятора, А*ч
    Cnom = ProjectsAttributes["Cnom"]

    #Температура окружающей среды, град С
    Tokr = ProjectsAttributes["Tokr"]

    #Степень разряда
    q0 = ProjectsAttributes["q0"]#Перенесенный через внешнюю цепь заряд (относительно номинальной емкости)

    #Метод интегрирования дифференциальных уравнений
    integrateMethod = ProjectsAttributes["integrateMethod"]
                               
    #Имя файла
    ParametersFileName = ProjectsAttributes["ParametersFileName"]#Файл csv параметров
    CurrentAttributesFileName = ProjectsAttributes["CurrentAttributesFileName"]#Файл csv аттрибутов тока
    AccumulatorAttributesFileName = ProjectsAttributes["AccumulatorAttributesFileName"]#Файл csv аттрибутов аккумулятора
    IntegrateAttributesFileName = ProjectsAttributes["IntegrateAttributesFileName"]#Файл csv аттрибутов интегрирования
    IndexesGraphicsFileName = ProjectsAttributes["IndexesGraphicsFileName"]#Файл csv индексов графиков, которые нужно построить
    DynamicFileName = ProjectsAttributes["DynamicFileName"]#Файл csv динамики
    sep = ProjectsAttributes["sep"]#Разделитель csv
    dec = ProjectsAttributes["dec"]#Десятичный разделитель
    GraphicFileDir = ProjectsAttributes["GraphicFileDir"]#Директория файлов графиков

    #Считываем файл аттрибутов тока
    currentAttributes = pd.read_csv(CurrentAttributesFileName,sep=sep,decimal=dec)

    #Считываем файл аттрибутов аккумулятора
    accumulatorAttributes = pd.read_csv(AccumulatorAttributesFileName,sep=sep,decimal=dec)

    #Размножаем параметры аккумулятора
    Pars = HStackMatrixRepeat([currentAttributes.to_numpy(),Tokr,  q0,accumulatorAttributes.to_numpy(),[Cnom]],
                              [                             True,True,                            True,  True])
    Pars = pd.DataFrame(Pars, columns=list(currentAttributes) + ["Tokr","q0"] + list(accumulatorAttributes) + ["Cnom"])

    #Считываем файл аттрибутов интегрирования
    integrateAttributes = pd.read_csv(IntegrateAttributesFileName,sep=sep,decimal=dec)

    #Получаем числа аттрибутов
    nDyns = len(Pars)#Число динамик
    nAccAttrs = len(accumulatorAttributes)#Число аттрибутов аккумулятора
    nCurAttrs = len(currentAttributes)#Число аттрибутов тока
    nTokrAttrs = len(Tokr)#Число аттрибутов температуры окружающей среди
    nQCapAttrs = len(q0)#Число аттрибутов отданной зарядовой емкости

    #Массив чисел аттрибутов
    arrAllAttrs = [nAccAttrs,nQCapAttrs,nTokrAttrs,nCurAttrs]

    #Вычисляем аттрибуты интегрирования
    integrateAttributes = IntegrateAttributes(integrateAttributes,
                                              arrAllAttrs,
                                              nDyns)#Приравниваем базовые аттрибуты интегрирования

    #Считываем файл индексов графиков
    indexesGraphics = pd.read_csv(IndexesGraphicsFileName,sep=sep,decimal=dec)

    #Формируем индексы динамик, графики которых мы будем строить
    (indexesGraphics,buildingGraphics) = IndexesGraphics(indexesGraphics,
                                                         arrAllAttrs)

    #Выполняем моделирование динамики системы
    dynamicIndex = LiionBatteryDynamicFunctionSaveVec(Pars,#Параметры
                                                      
                                                      #Параметры интегрирования
                                                      integrateAttributes,#Аттрибуты интегрирования
                                                      integrateMethod,#Метод интегрирования дифференциальных уравнений
                                                         
                                                      #Построение графиков
                                                      indexesGraphics,#Индексы графиков, которые нужно построить
                                                      buildingGraphics,#Необходимость построения графиков
                                                      GraphicFileDir,#Директория файлов графиков
                                                         
                                                      #Индикация динамики
                                                      modelingDynamicIndicate,
                                                         
                                                      #Имя файла
                                                      DynamicFileName,#Файл csv
                                                      sep,#Разделитель csv
                                                      dec#Десятичный разделитель
                                                      )

    #Сохраняем параметры
    LiionBatteryParametersSave(dynamicIndex,#Индексы динамик
                               Pars,#Параметры
                               ParametersFileName,#Имя файла параметров
                               sep,#CSV разделитель
                               dec#Десятичный разделитель
                               )
