from MathProtEnergyProc import NonEqSystemQ
from MathProtEnergyProc import NonEqSystemQDyn
from MathProtEnergyProc import CountDynamicsQ
from MathProtEnergyProc.IntegrateDyn import standartIntegrateDyn

from LiionBattery.LiionBatteryFunctions.LiionBatteryFunctionStation import LiionBatteryStateFunction
from LiionBattery.LiionBatteryFunctions.LiionBatteryFunctionCharacteristics import LiionBatteryCharacteristicsFunction
from LiionBattery.LiionBatteryFunctions.fU import fU

#Функция расчета динамики
def LiionBatteryStructureDynamic(integrateMethod#Метод интегрирования дифференциальных уравнений
                                 ):#Структура для расчета одной динамики
    #Описываем структуру литий-ионного элемента
    stateCoordinatesNames = ["qbinp","qm","qbinn","q"]#Имена координат состояния
    processCoordinatesNames = ["dqbinp","dqm","dqbinn"]#Имена координат процессов
    energyPowersNames = ["EnPowInAkk","EnPowBAkk","EnPowOkr"]#Имена энергетических степеней свободы
    reducedTemperaturesEnergyPowersNames = ["TInAkk","TBAkk"]#Имена приведенных температур энергетических степеней свободы
    energyPowersBetNames = []#Имена взаимодействий между энергетическими степенями свободы
    heatTransfersNames = ["QInBAkk","QBAkkExp"]#Имена потоков переноса теплоты
    heatTransfersOutputEnergyPowersNames = ["EnPowInAkk","EnPowBAkk"]#Имена энергетических степеней свободы, с которых уходит теплота
    heatTransfersInputEnergyPowersNames = ["EnPowBAkk","EnPowOkr"]#Имена энергетических степеней свободы, на которые приходит теплота
    stateCoordinatesStreamsNames = ["qbinp","qm","qbinn","q"]#Имена координат состояния, изменяемых в результате внешних потоков
    heatEnergyPowersStreamsNames = ["EnPowBAkk"]#Имена потоков теплоты на энергетические степени свободы
    stateFunction = LiionBatteryStateFunction#Функция состояния
    stateCoordinatesVarBalanceNames = []#Имена переменных коэффициентов матрицы баланса по координатам состояния
    processCoordinatesVarBalanceNames = []#Имена переменных коэффициентов матрицы баланса по координатам процессов
    energyPowersVarTemperatureNames = ["EnPowInAkk","EnPowBAkk","EnPowOkr"]#Имена переменных температур энергетических степеней свободы
    stateCoordinatesVarPotentialsInterNames = ["qbinp","qm","qbinn","q"]#Имена переменных потенциалов взаимодействия по координатам состояния
    energyPowersVarPotentialsInterNames = ["EnPowInAkk","EnPowInAkk","EnPowInAkk","EnPowInAkk"]#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
    stateCoordinatesVarPotentialsInterBetNames = []#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по координатам состояния
    energyPowersVarPotentialsInterBetNames = []#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по энергетическим степеням свободы
    energyPowersVarBetaNames = []#Имена переменных долей распределения некомпенсированной теплоты энергетических степеней свободы
    processCoordinatesVarBetaNames = []#Имена переменных долей распределения некомпенсированной теплоты координат процессов
    reducedTemperaturesEnergyPowersVarInvHeatCapacityNames = ["TInAkk","TBAkk"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
    energyPowersVarInvHeatCapacityNames = ["EnPowInAkk","EnPowBAkk"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
    reducedTemperaturesEnergyPowersVarHeatEffectNames = ["TInAkk","TInAkk","TInAkk","TInAkk"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
    stateCoordinatesVarHeatEffectNames = ["qbinp","qm","qbinn","q"]#Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
    varKineticPCPCNames = ["dqbinp","dqm","dqbinn"]#Имена сопряженностей между собой координат процессов
    varKineticPCPCAffNames = ["dqbinp","dqm","dqbinn"]#Имена сопряженностей между собой термодинамических сил
    varKineticPCHeatNames = []#Имена сопряженностей координат процессов с теплопереносами
    varKineticPCHeatAffNames = []#Имена сопряженностей термодинамических сил с теплопереносами
    varKineticHeatPCNames = []#Имена сопряженностей теплопереносов с координатами процессов
    varKineticHeatPCAffNames = []#Имена сопряженностей теплопереносов с термодинамическими силами
    varKineticHeatHeatNames = ["QInBAkk","QBAkkExp"]#Имена сопряженностей между собой перенесенных теплот
    varKineticHeatHeatAffNames = ["QInBAkk","QBAkkExp"]#Имена сопряженностей между собой термодинамических сил по переносу теплот
    stateCoordinatesVarStreamsNames = ["qbinp","qm","qbinn","q"]#Имена переменных внешних потоков
    heatEnergyPowersVarStreamsNames = ["EnPowBAkk"]#Имена переменных внешних потоков теплоты
                     
    #Литий-ионный аккумулятор
    LiionAkk = NonEqSystemQ(stateCoordinatesNames,#Имена координат состояния
                            processCoordinatesNames,#Имена координат процессов
                            energyPowersNames,#Имена энергетических степеней свободы
                            reducedTemperaturesEnergyPowersNames,#Имена приведенных температур энергетических степеней свободы
                            energyPowersBetNames,#Имена взаимодействий между энергетическими степенями свободы
                            heatTransfersNames,#Имена потоков переноса теплоты
                            heatTransfersOutputEnergyPowersNames,#Имена энергетических степеней свободы, с которых уходит теплота
                            heatTransfersInputEnergyPowersNames,#Имена энергетических степеней свободы, на которые приходит теплота
                             
                            stateCoordinatesStreamsNames,#Имена координат состояния, изменяемых в результате внешних потоков
                            heatEnergyPowersStreamsNames,#Имена потоков теплоты на энергетические степени свободы
                             
                            stateFunction,#Функция состояния
                             
                            stateCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам состояния
                            processCoordinatesVarBalanceNames,#Имена переменных коэффициентов матрицы баланса по координатам процессов
                            energyPowersVarTemperatureNames,#Имена переменных температур энергетических степеней свободы
                            stateCoordinatesVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по координатам состояния
                            energyPowersVarPotentialsInterNames,#Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
                            stateCoordinatesVarPotentialsInterBetNames,#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по координатам состояния
                            energyPowersVarPotentialsInterBetNames,#Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по энергетическим степеням свободы
                            energyPowersVarBetaNames,#Имена переменных долей распределения некомпенсированной теплоты энергетических степеней свободы
                            processCoordinatesVarBetaNames,#Имена переменных долей распределения некомпенсированной теплоты координат процессов
                            reducedTemperaturesEnergyPowersVarInvHeatCapacityNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
                            energyPowersVarInvHeatCapacityNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
                            reducedTemperaturesEnergyPowersVarHeatEffectNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
                            stateCoordinatesVarHeatEffectNames,#Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
                             
                            varKineticPCPCNames,#Имена сопряженностей между собой координат процессов
                            varKineticPCPCAffNames,#Имена сопряженностей между собой термодинамических сил
                            varKineticPCHeatNames,#Имена сопряженностей координат процессов с теплопереносами
                            varKineticPCHeatAffNames,#Имена сопряженностей термодинамических сил с теплопереносами
                            varKineticHeatPCNames,#Имена сопряженностей теплопереносов с координатами процессов
                            varKineticHeatPCAffNames,#Имена сопряженностей теплопереносов с термодинамическими силами
                            varKineticHeatHeatNames,#Имена сопряженностей между собой перенесенных теплот
                            varKineticHeatHeatAffNames,#Имена сопряженностей между собой термодинамических сил по переносу теплот
                             
                            stateCoordinatesVarStreamsNames,#Имена переменных внешних потоков
                            heatEnergyPowersVarStreamsNames#Имена переменных внешних потоков теплоты
                            )

    #Задаем связь между коордиинатами состояния и процессами
    LiionAkk.SetBalanceStateCoordinatesConstElement("qbinp","dqbinp",1)
    LiionAkk.SetBalanceStateCoordinatesConstElement("qm","dqm",1)
    LiionAkk.SetBalanceStateCoordinatesConstElement("qbinn","dqbinn",1)

    #Задаем доли распределения некомпенсированной теплоты
    LiionAkk.SetBetaConstElement("EnPowInAkk","dqbinp",1.0)
    LiionAkk.SetBetaConstElement("EnPowInAkk","dqm",1.0)
    LiionAkk.SetBetaConstElement("EnPowInAkk","dqbinn",1.0)

    #Задаем класс динамики системы
    integDynamic = standartIntegrateDyn(method = integrateMethod)
    return NonEqSystemQDyn(LiionAkk,#литий-ионный аккумулятор
                           fU,#Функция условий протекания процессов
                           LiionBatteryCharacteristicsFunction,#Функция внешних параметров
                           integDynamic#Метод интегрирования дифференциальных уравнений
                           )
def LiionBatteryStructureDynamics(integrateMethod,#Метод интегрирования дифференциальных уравнений
                                  SavedFinction#Функция сохранения в файл
                                  ):#Структура для расчета многих динамик
    #Создаем класс литий-ионного аккумулятора с расчетом хотя бы одной динамики
    LiionAkkDyn = LiionBatteryStructureDynamic(integrateMethod)
    
    #Возвращаем класс расчета динамик
    return CountDynamicsQ(LiionAkkDyn,SavedFinction)
