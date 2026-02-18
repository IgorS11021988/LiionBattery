import numpy as np

from .StationFunction import IndepStateFunction

from MathProtEnergyProc.HeatPowerValues import IntPotentialsOne, HeatValuesOne

from MathProtEnergyProc.CorrectionModel import PosLinearFilter


# Потенциалы взаимодействия в топливном элементе и камерах
potentialInterAkk = IntPotentialsOne(["qbinp", "qm", "qbinn", "q"],  # Имена координат состояния
                                     ["EnPowInAkk", "EnPowBAkk"],  # Имена энергетических степеней свободы

                                     [     "qbinp",         "qm",      "qbinn",          "q"],  # Имена переменных потенциалов взаимодействия по координатам состояния
                                     ["EnPowInAkk", "EnPowInAkk", "EnPowInAkk", "EnPowInAkk"]  # Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
                                     )

# Приведенные обратные теплоемкости и тепловые эффекты
heatValuesAkk = HeatValuesOne(["qbinp", "qm", "qbinn", "q"],  # Имена координат состояния
                              ["EnPowInAkk", "EnPowBAkk"],  # Имена энергетических степеней свободы
          
                              ["EnPowInAkk", "EnPowBAkk"],  # Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
                              ["EnPowInAkk", "EnPowInAkk", "EnPowInAkk", "EnPowInAkk"],  # Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
                              [     "qbinp",         "qm",      "qbinn",          "q"]  # Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
                              )


# Функция состояния для литий-ионного аккумулятора
def StateFunction(stateCoordinates,
                  reducedTemp,
                  systemParameters):
    # Получаем независимые составляющие свойств веществ и процессов
    (I, Tokr,
     heatStreambEnPow,
     JSq, JST, HSqT, HSTT,
     rAkk, KQAkk) = IndepStateFunction(stateCoordinates,
                                       reducedTemp,
                                       systemParameters)

    # Матрица баланса
    balanceMatrix = np.array([])

    # Внешние потоки зарядов
    stateCoordinatesStreams = np.array([-I, -I, -I, I], dtype=np.double)

    # Внешние потоки теплоты
    heatEnergyPowersStreams = np.array([heatStreambEnPow], dtype=np.double)

    # Выводим температуры
    energyPowerTemperatures = np.hstack([reducedTemp, [Tokr]])

    # Потенциалы взаимодействия энергетических степеней свободы
    potentialInter = potentialInterAkk(JSq, reducedTemp)

    # Потенциалы взаимодействия между энергетическими степенями свободы
    potentialInterBet = np.array([])

    # Доли распределения некомпенсированной теплоты
    beta = np.array([])

    # Главный блок кинетической матрицы по процессам
    kineticMatrixPCPC = 1 / PosLinearFilter(rAkk)

    # Перекрестные блоки кинетической матрицы по процессам
    kineticMatrixPCHeat = np.array([])
    kineticMatrixHeatPC = np.array([])

    # Главный блок кинетической матрицы по теплообмену
    kineticMatrixHeatHeat = PosLinearFilter(KQAkk)

    # Обратная теплоемкость и приведенные тепловые эффекты литий-ионного аккумулятора
    (invHeatCapacityMatrixCf,  # Обратная теплоемкость водородно-воздушного топливного элемента
     heatEffectMatrixCf  # Приведенные тепловые эффекты водородно-воздушного топливного элемента
     ) = heatValuesAkk(JST,  # Якобиан приведенной энтропии по температурам
                       HSTT,  # Матрица Гесса приведенной энтропии по температурам
                       HSqT,  # Матрица Гесса приведенной энтропии по температурам и координатам состояния
                       reducedTemp  # Температуры
                       )

    # Выводим результат
    return (balanceMatrix,
            stateCoordinatesStreams,
            heatEnergyPowersStreams,
            energyPowerTemperatures,
            potentialInter,
            potentialInterBet,
            beta, kineticMatrixPCPC,
            kineticMatrixPCHeat,
            kineticMatrixHeatPC,
            kineticMatrixHeatHeat,
            invHeatCapacityMatrixCf,
            heatEffectMatrixCf)


# Функция структуры аккумулятора
def StructureFunction():
    # Описываем структуру литий-ионного элемента
    stateCoordinatesNames = ["qbinp", "qm", "qbinn", "q"]  # Имена координат состояния
    processCoordinatesNames = ["dqbinp", "dqm", "dqbinn"]  # Имена координат процессов
    energyPowersNames = ["EnPowInAkk", "EnPowBAkk", "EnPowOkr"]  # Имена энергетических степеней свободы
    reducedTemperaturesEnergyPowersNames = ["TInAkk", "TBAkk"]  # Имена приведенных температур энергетических степеней свободы
    energyPowersBetNames = []  # Имена взаимодействий между энергетическими степенями свободы
    heatTransfersNames = ["QInBAkk", "QBAkkExp"]  # Имена потоков переноса теплоты
    heatTransfersOutputEnergyPowersNames = ["EnPowInAkk", "EnPowBAkk"]  # Имена энергетических степеней свободы, с которых уходит теплота
    heatTransfersInputEnergyPowersNames = ["EnPowBAkk", "EnPowOkr"]  # Имена энергетических степеней свободы, на которые приходит теплота
    stateCoordinatesStreamsNames = ["qbinp", "qm", "qbinn", "q"]  # Имена координат состояния, изменяемых в результате внешних потоков
    heatEnergyPowersStreamsNames = ["EnPowBAkk"]  # Имена потоков теплоты на энергетические степени свободы
    stateFunction = StateFunction  # Функция состояния
    stateCoordinatesVarBalanceNames = []  # Имена переменных коэффициентов матрицы баланса по координатам состояния
    processCoordinatesVarBalanceNames = []  # Имена переменных коэффициентов матрицы баланса по координатам процессов
    energyPowersVarTemperatureNames = ["EnPowInAkk", "EnPowBAkk", "EnPowOkr"]  # Имена переменных температур энергетических степеней свободы
    stateCoordinatesVarPotentialsInterNames = ["qbinp", "qm", "qbinn", "q"]  # Имена переменных потенциалов взаимодействия по координатам состояния
    energyPowersVarPotentialsInterNames = ["EnPowInAkk", "EnPowInAkk", "EnPowInAkk", "EnPowInAkk"]  # Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
    stateCoordinatesVarPotentialsInterBetNames = []  # Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по координатам состояния
    energyPowersVarPotentialsInterBetNames = []  # Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по энергетическим степеням свободы
    energyPowersVarBetaNames = []  # Имена переменных долей распределения некомпенсированной теплоты энергетических степеней свободы
    processCoordinatesVarBetaNames = []  # Имена переменных долей распределения некомпенсированной теплоты координат процессов
    reducedTemperaturesEnergyPowersVarInvHeatCapacityNames = ["TInAkk", "TBAkk"]  # Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
    energyPowersVarInvHeatCapacityNames = ["EnPowInAkk", "EnPowBAkk"]  # Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
    reducedTemperaturesEnergyPowersVarHeatEffectNames = ["TInAkk", "TInAkk", "TInAkk", "TInAkk"]  # Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
    stateCoordinatesVarHeatEffectNames = ["qbinp", "qm", "qbinn", "q"]  # Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
    varKineticPCPCNames = ["dqbinp", "dqm", "dqbinn"]  # Имена сопряженностей между собой координат процессов
    varKineticPCPCAffNames = ["dqbinp", "dqm", "dqbinn"]  # Имена сопряженностей между собой термодинамических сил
    varKineticPCHeatNames = []  # Имена сопряженностей координат процессов с теплопереносами
    varKineticPCHeatAffNames = []  # Имена сопряженностей термодинамических сил с теплопереносами
    varKineticHeatPCNames = []  # Имена сопряженностей теплопереносов с координатами процессов
    varKineticHeatPCAffNames = []  # Имена сопряженностей теплопереносов с термодинамическими силами
    varKineticHeatHeatNames = ["QInBAkk", "QBAkkExp"]  # Имена сопряженностей между собой перенесенных теплот
    varKineticHeatHeatAffNames = ["QInBAkk", "QBAkkExp"]  # Имена сопряженностей между собой термодинамических сил по переносу теплот
    stateCoordinatesVarStreamsNames = ["qbinp", "qm", "qbinn", "q"]  # Имена переменных внешних потоков
    heatEnergyPowersVarStreamsNames = ["EnPowBAkk"]  # Имена переменных внешних потоков теплоты

    # Выводим структуру литий-ионного аккумулятора
    return (stateCoordinatesNames,  # Имена координат состояния
            processCoordinatesNames,  # Имена координат процессов
            energyPowersNames,  # Имена энергетических степеней свободы
            reducedTemperaturesEnergyPowersNames,  # Имена приведенных температур энергетических степеней свободы
            energyPowersBetNames,  # Имена взаимодействий между энергетическими степенями свободы
            heatTransfersNames,  # Имена потоков переноса теплоты
            heatTransfersOutputEnergyPowersNames,  # Имена энергетических степеней свободы, с которых уходит теплота
            heatTransfersInputEnergyPowersNames,  # Имена энергетических степеней свободы, на которые приходит теплота
            stateCoordinatesStreamsNames,  # Имена координат состояния, изменяемых в результате внешних потоков
            heatEnergyPowersStreamsNames,  # Имена потоков теплоты на энергетические степени свободы
            stateFunction,  # Функция состояния
            stateCoordinatesVarBalanceNames,  # Имена переменных коэффициентов матрицы баланса по координатам состояния
            processCoordinatesVarBalanceNames,  # Имена переменных коэффициентов матрицы баланса по координатам процессов
            energyPowersVarTemperatureNames,  # Имена переменных температур энергетических степеней свободы
            stateCoordinatesVarPotentialsInterNames,  # Имена переменных потенциалов взаимодействия по координатам состояния
            energyPowersVarPotentialsInterNames,  # Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
            stateCoordinatesVarPotentialsInterBetNames,  # Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по координатам состояния
            energyPowersVarPotentialsInterBetNames,  # Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по энергетическим степеням свободы
            energyPowersVarBetaNames,  # Имена переменных долей распределения некомпенсированной теплоты энергетических степеней свободы
            processCoordinatesVarBetaNames,  # Имена переменных долей распределения некомпенсированной теплоты координат процессов
            reducedTemperaturesEnergyPowersVarInvHeatCapacityNames,  # Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
            energyPowersVarInvHeatCapacityNames,  # Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
            reducedTemperaturesEnergyPowersVarHeatEffectNames,  # Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
            stateCoordinatesVarHeatEffectNames,  # Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
            varKineticPCPCNames,  # Имена сопряженностей между собой координат процессов
            varKineticPCPCAffNames,  # Имена сопряженностей между собой термодинамических сил
            varKineticPCHeatNames,  # Имена сопряженностей координат процессов с теплопереносами
            varKineticPCHeatAffNames,  # Имена сопряженностей термодинамических сил с теплопереносами
            varKineticHeatPCNames,  # Имена сопряженностей теплопереносов с координатами процессов
            varKineticHeatPCAffNames,  # Имена сопряженностей теплопереносов с термодинамическими силами
            varKineticHeatHeatNames,  # Имена сопряженностей между собой перенесенных теплот
            varKineticHeatHeatAffNames,  # Имена сопряженностей между собой термодинамических сил по переносу теплот
            stateCoordinatesVarStreamsNames,  # Имена переменных внешних потоков
            heatEnergyPowersVarStreamsNames  # Имена переменных внешних потоков теплоты
            )


# Функция постоянных параметров литий-ионного аккумулятора
def ConstParametersFunction(sysStructure  # Структура системы
                            ):
    # Задаем связь между коордиинатами состояния и процессами
    sysStructure.SetBalanceStateCoordinatesConstElement("qbinp", "dqbinp", 1)
    sysStructure.SetBalanceStateCoordinatesConstElement("qm", "dqm", 1)
    sysStructure.SetBalanceStateCoordinatesConstElement("qbinn", "dqbinn", 1)

    # Задаем доли распределения некомпенсированной теплоты
    sysStructure.SetBetaConstElement("EnPowInAkk", "dqbinp", 1.0)
    sysStructure.SetBetaConstElement("EnPowInAkk", "dqm", 1.0)
    sysStructure.SetBetaConstElement("EnPowInAkk", "dqbinn", 1.0)
