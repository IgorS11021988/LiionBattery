import numpy as np

from .StationFunction import IndepStateFunction

from MathProtEnergyProc.HeatPowerValues import IntPotentialsOne, HeatValuesOne

from MathProtEnergyProc.CorrectionModel import PosLinearFilter, ReluFilter, KineticMatrixQ, KineticMatrixFromFacStreamEkvAff


# Потенциалы взаимодействия в топливном элементе и камерах
potentialInterAkk = IntPotentialsOne(["qbinp", "qm", "qbinn", "q", "qMatElp", "qMatEln", "qMatDegElp", "qMatDegEln", "qDegPosEl"],  # Имена координат состояния
                                     ["EnPowInAkk", "EnPowBAkk"],  # Имена энергетических степеней свободы

                                     [     "qbinp",         "qm",      "qbinn",          "q",    "qMatElp",    "qMatEln", "qMatDegElp", "qMatDegEln",  "qDegPosEl"],  # Имена переменных потенциалов взаимодействия по координатам состояния
                                     ["EnPowInAkk", "EnPowInAkk", "EnPowInAkk", "EnPowInAkk", "EnPowInAkk", "EnPowInAkk", "EnPowInAkk", "EnPowInAkk", "EnPowInAkk"]  # Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
                                     )

# Приведенные обратные теплоемкости и тепловые эффекты
heatValuesAkk = HeatValuesOne(["qbinp", "qm", "qbinn", "q", "qMatElp", "qMatEln", "qMatDegElp", "qMatDegEln", "qDegPosEl"],  # Имена координат состояния
                              ["EnPowInAkk", "EnPowBAkk"],  # Имена энергетических степеней свободы
          
                              ["EnPowInAkk", "EnPowBAkk"],  # Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
                              ["EnPowInAkk", "EnPowInAkk", "EnPowInAkk", "EnPowInAkk", "EnPowInAkk", "EnPowInAkk", "EnPowInAkk", "EnPowInAkk", "EnPowInAkk"],  # Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
                              [     "qbinp",         "qm",      "qbinn",          "q",    "qMatElp",    "qMatEln", "qMatDegElp", "qMatDegEln",  "qDegPosEl"]  # Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
                              )


# Кинетические матрицы по электродам
kinMatrixElp = KineticMatrixQ(["dqbinp", "deactp", "deactp"],  # Имена сопряженностей между собой координат процессов
                              ["dqbinp", "dqbinp", "deactp"],  # Имена сопряженностей между собой термодинамических сил
                              [],  # Имена сопряженностей координат процессов с теплопереносами
                              [],  # Имена сопряженностей термодинамических сил с теплопереносами
                              [],  # Имена сопряженностей теплопереносов с координатами процессов
                              [],  # Имена сопряженностей теплопереносов с термодинамическими силами
                              [],  # Имена сопряженностей между собой перенесенных теплот
                              [],  # Имена сопряженностей между собой термодинамических сил по переносу теплот

                              [["deactp", "dqbinp"]]  # Массив массивов имен координат процессов (в том числе и перенесенных теплот) по кинетической матрице
                              )  # Положительный электрод
kinMatrixEln = KineticMatrixQ(["dqbinn", "deactn", "deactn"],  # Имена сопряженностей между собой координат процессов
                              ["dqbinn", "dqbinn", "deactn"],  # Имена сопряженностей между собой термодинамических сил
                              [],  # Имена сопряженностей координат процессов с теплопереносами
                              [],  # Имена сопряженностей термодинамических сил с теплопереносами
                              [],  # Имена сопряженностей теплопереносов с координатами процессов
                              [],  # Имена сопряженностей теплопереносов с термодинамическими силами
                              [],  # Имена сопряженностей между собой перенесенных теплот
                              [],  # Имена сопряженностей между собой термодинамических сил по переносу теплот

                              [["deactn", "dqbinn"]]  # Массив массивов имен координат процессов (в том числе и перенесенных теплот) по кинетической матрице
                              )  # Отрицательный электрод


# Функция состояния для литий-ионного аккумулятора
def StateFunction(stateCoordinates,
                  reducedTemp,
                  systemParameters):
    # Получаем независимые составляющие свойств веществ и процессов
    (I, Tokr,
     heatStreambEnPow,
     JSq, JST, HSqT, HSTT,
     rbinp, rbinn, rm,
     aActp, aActn,
     kActp, kActn,
     KDegEl, KQAkk) = IndepStateFunction(stateCoordinates,
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

    # Кинетическая матрица положительного электрода
    sbinp = np.array([[1 / PosLinearFilter(rbinp)]], dtype=np.double)
    kMatrElp = KineticMatrixFromFacStreamEkvAff(ReluFilter(kActp),  # Необратимая составляющая кинетической матрицы
                                                aActp.reshape(1,1),  # Коэффициенты увлечения потоков
                                                np.zeros_like(aActp),  # Коэффициенты эквивалетность термодинаических сил
                                                sbinp  # Блок кинетической матрицы
                                                )
    (kineticMatrixPCPCElp,
     kineticMatrixPCHeatElp,
     kineticMatrixHeatPCElp,
     kineticMatrixHeatHeatElp) = kinMatrixElp([kMatrElp])

    # Кинетическая матрица отрицательного электрода
    sbinn = np.array([[1 / PosLinearFilter(rbinn)]], dtype=np.double)
    kMatrEln = KineticMatrixFromFacStreamEkvAff(ReluFilter(kActn),  # Необратимая составляющая кинетической матрицы
                                                aActn.reshape(1,1),  # Коэффициенты увлечения потоков
                                                np.zeros_like(aActn),  # Коэффициенты эквивалетность термодинаических сил
                                                sbinn  # Блок кинетической матрицы
                                                )
    (kineticMatrixPCPCEln,
     kineticMatrixPCHeatEln,
     kineticMatrixHeatPCEln,
     kineticMatrixHeatHeatEln) = kinMatrixEln([kMatrEln])

    # Главный блок кинетической матрицы по процессам
    kineticMatrixPCPC = np.hstack([kineticMatrixPCPCElp,
                                   kineticMatrixPCPCEln,
                                   [1 / PosLinearFilter(rm)],
                                   ReluFilter(KDegEl)])

    # Перекрестные блоки кинетической матрицы по процессам
    kineticMatrixPCHeat = np.hstack([kineticMatrixPCHeatElp,
                                     kineticMatrixPCHeatEln])
    kineticMatrixHeatPC = np.hstack([kineticMatrixHeatPCElp,
                                     kineticMatrixHeatPCEln])

    # Главный блок кинетической матрицы по теплообмену
    kineticMatrixHeatHeat = np.hstack([kineticMatrixHeatHeatElp,
                                       kineticMatrixHeatHeatEln,
                                       ReluFilter(KQAkk)])

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
    stateCoordinatesNames = ["qbinp", "qm", "qbinn", "q", "qMatElp", "qMatEln", "qMatDegElp", "qMatDegEln", "qDegPosEl"]  # Имена координат состояния
    processCoordinatesNames = ["dqbinp", "dqm", "dqbinn", "deactp", "deactn", "degp", "degn", "dqDegPosEl"]  # Имена координат процессов
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
    stateCoordinatesVarPotentialsInterNames = ["qbinp", "qm", "qbinn", "q", "qMatElp", "qMatEln", "qMatDegElp", "qMatDegEln", "qDegPosEl"]  # Имена переменных потенциалов взаимодействия по координатам состояния
    energyPowersVarPotentialsInterNames = ["EnPowInAkk", "EnPowInAkk", "EnPowInAkk", "EnPowInAkk", "EnPowInAkk", "EnPowInAkk", "EnPowInAkk", "EnPowInAkk", "EnPowInAkk"]  # Имена переменных потенциалов взаимодействия по энергетическим степеням свободы
    stateCoordinatesVarPotentialsInterBetNames = []  # Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по координатам состояния
    energyPowersVarPotentialsInterBetNames = []  # Имена переменных потенциалов взаимодействия для взаимодействий между энергетическими степенями свободы по энергетическим степеням свободы
    energyPowersVarBetaNames = []  # Имена переменных долей распределения некомпенсированной теплоты энергетических степеней свободы
    processCoordinatesVarBetaNames = []  # Имена переменных долей распределения некомпенсированной теплоты координат процессов
    reducedTemperaturesEnergyPowersVarInvHeatCapacityNames = ["TInAkk", "TBAkk"]  # Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
    energyPowersVarInvHeatCapacityNames = ["EnPowInAkk", "EnPowBAkk"]  # Имена переменных коэффициентов обратных теплоемкостей по отношению к энергетическим степеням свободы
    reducedTemperaturesEnergyPowersVarHeatEffectNames = ["TInAkk", "TInAkk", "TInAkk", "TInAkk", "TInAkk", "TInAkk", "TInAkk", "TInAkk", "TInAkk"]  # Имена переменных коэффициентов обратных теплоемкостей по отношению к приведенным температурам
    stateCoordinatesVarHeatEffectNames = ["qbinp", "qm", "qbinn", "q", "qMatElp", "qMatEln", "qMatDegElp", "qMatDegEln", "qDegPosEl"]  # Имена переменных коэффициентов обратных теплоемкостей по отношению к координатам состояния
    varKineticPCPCNames = ["dqbinp", "deactp", "deactp", "dqbinn", "deactn", "deactn", "dqm", "degp", "degn", "dqDegPosEl"]  # Имена сопряженностей между собой координат процессов
    varKineticPCPCAffNames = ["dqbinp", "dqbinp", "deactp", "dqbinn", "dqbinn", "deactn", "dqm", "degp", "degn", "dqDegPosEl"]  # Имена сопряженностей между собой термодинамических сил
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
    sysStructure.SetBalanceStateCoordinatesConstElement("qMatElp", "deactp", 1)
    sysStructure.SetBalanceStateCoordinatesConstElement("qMatEln", "deactn", 1)
    sysStructure.SetBalanceStateCoordinatesConstElement("qMatDegElp", "degp", 1)
    sysStructure.SetBalanceStateCoordinatesConstElement("qMatDegEln", "degn", 1)
    sysStructure.SetBalanceStateCoordinatesConstElement("qbinp", "dqDegPosEl", 1)
    sysStructure.SetBalanceStateCoordinatesConstElement("qDegPosEl", "dqDegPosEl", 1)

    # Задаем доли распределения некомпенсированной теплоты
    sysStructure.SetBetaConstElement("EnPowInAkk", "dqbinp", 1.0)
    sysStructure.SetBetaConstElement("EnPowInAkk", "dqm", 1.0)
    sysStructure.SetBetaConstElement("EnPowInAkk", "dqbinn", 1.0)
    sysStructure.SetBetaConstElement("EnPowInAkk", "deactp", 1.0)
    sysStructure.SetBetaConstElement("EnPowInAkk", "deactn", 1.0)
    sysStructure.SetBetaConstElement("EnPowInAkk", "degp", 1.0)
    sysStructure.SetBetaConstElement("EnPowInAkk", "degn", 1.0)
    sysStructure.SetBetaConstElement("EnPowInAkk", "dqDegPosEl", 1.0)
