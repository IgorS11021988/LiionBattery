import numpy as np

from LiionBattery.LiionBatteryDynamics.LiionBatteryStructure import LiionBatteryStructureDynamics
from LiionBattery.LiionBatteryDynamics.LiionBatteryValues import LiionBatteryInputArrayCreate, LiionBatteryOutputValues

from MathProtEnergyProc.RepeatLastRowsMatrix import RepeatLastRowsMatrix
from MathProtEnergyProc.IndexedNames import IndexedNamesFromIndexes

#Функция расчета динамик
def LiionBatteryDynamicFunctionSaveVec(Pars,#Параметры
                                       
                                       #Параметры интегрирования
                                       integrateAttributes,#Аттрибуты интегрирования
                                       integrateMethod,#Метод интегрирования дифференциальных уравнений
                                          
                                       #Построение графиков
                                       indexesGraphics,#Индексы графиков, которые нужно построить
                                       buildingGraphics,#Необходимость построения графиков
                                       graphicFileDir,#Директория файлов графиков
                                                
                                       #Индикация динамики
                                       modelingDynamicIndicate,
                                                               
                                       #Имя файла
                                       DynamicFileNameBase,#Базовое имя файла csv
                                       sep,#Разделитель csv
                                       dec#Десятичный разделитель
                                       ):
    #Исходные данные моделирования системы
    (Tints,
     stateCoordinates0s,
     reducedTemp0s,
     systemParameters,
     ts) = LiionBatteryInputArrayCreate(Pars,#Параметры
                                        
                                        integrateAttributes#Аттрибуты интегрирования
                                        )

    #Функция сохранения в файл
    nDyns = len(Tints)#Число динамик
    def SavedFinction(dyn, index):
        #Формируем имя файла
        fileName = IndexedNamesFromIndexes([index + 1],#Индексы
                                           DynamicFileNameBase,#Начало имени
                                           endName = ".csv",#Конец имени
                                           sepName = "_"#Раздлитель имени
                                           )[0]
        
        #Сохраняем данные в файл
        BuildGraphic = (buildingGraphics and np.any((index + 1) == indexesGraphics))#Необходимость построения графика
        modelingDynamicIndicate(index + 1,nDyns)#Выводим текущую динамику
        LiionBatteryOutputValues(dyn,fileName,sep,dec,
                                 graphicFileDir,index + 1,
                                 plotGraphics=BuildGraphic)
        
        #Возвращаем имя файла
        return index
    
    #Задаем класс динамик системы
    LiionAkkDyns = LiionBatteryStructureDynamics(integrateMethod,#Метод интегрирования дифференциальных уравнений
                                                 SavedFinction#Функция сохранения в файл
                                                 )
    
    #Получаем динамики
    return LiionAkkDyns.ComputingExperimentQ(Tints,
                                             stateCoordinates0s,
                                             reducedTemp0s,
                                             systemParameters,
                                             t_evals = ts)
 