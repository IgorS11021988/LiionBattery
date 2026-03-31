import numpy as np

from MathProtEnergyProcBase.IndexFunctions import GetIndex, GetIndexes

from .AttributesNames import stateCoordinatesNames, reducedTemperaturesEnergyPowersNames, USystemParametersNames, otherSystemParametersNames, processCoordinatesNames
from .StationFunctions import funCbin


# Индексы координат состояния
qInd = GetIndex(stateCoordinatesNames, "q")  # Индекс перенесенного через внешнюю цепь электрического заряда
qMatElpInd = GetIndex(stateCoordinatesNames, "qMatElp")  # Индекс зарядового числа молей недеградированного положительного электрода
qMatElnInd = GetIndex(stateCoordinatesNames, "qMatEln")  # Индекс зарядового числа молей недеградированного отрицательного электрода
qMatDegElpInd = GetIndex(stateCoordinatesNames, "qMatDegElp")  # Индекс зарядового числа молей деградированного положительного электрода
qMatDegElnInd = GetIndex(stateCoordinatesNames, "qMatDegEln")  # Индекс зарядового числа молей деградированного отрицательного электрода
qDegPosElInd = GetIndex(stateCoordinatesNames, "qDegPosEl")  # Индекс зарядового числа молей разрушенного положительного электрода

# Индексы приведенных температур
TInAkkInd = GetIndex(reducedTemperaturesEnergyPowersNames, "TInAkk")  # Индекс температуры содержимого аккумулятора
TBAkkInd = GetIndex(reducedTemperaturesEnergyPowersNames, "TBAkk")  # Индекс температуры корпуса аккумулятора

# Индексы токов в элементе
IElInd = GetIndexes(processCoordinatesNames, ["dqbinp", "dqm", "dqbinn"])

# Индексы параметров системы
otherSystemParametersInd = GetIndexes(otherSystemParametersNames, ["Cnom",
                                                                   "Rkl",
                                                                   "Tokr",
                                                                   "qMatAllp",
                                                                   "qMatAlln"])


# Функция состояния для литий-ионного аккумулятора
def CharacteristicsFunction(t,  # Моменты времени
                            stateCoordinates,  # Координаты состояния
                            reducedTemp,  # Приведенные температуры
                            USystemParameters,  # U-параметры системы
                            otherSystemParameters,  # Прочие параметры системы
                            nEqSysQ  # Указатель на функцию системы
                            ):
    # Рассчитываем аттрибуты системы
    def GetElAttr(ind):
        # Рассчитываем аттрибуты
        nEqSysQ.CountSystem(stateCoordinates[ind],  # Координаты состояния
                            reducedTemp[ind],  # Приведенные температуры энергетических степеней свободы
                            np.hstack([USystemParameters[ind], otherSystemParameters])  # Параметры системы
                            )

        # Получаем напряжения двойных слоев
        Ubinp = nEqSysQ.GetStateFunction().GetIndepStateFunction().GetUbinp()  # Положительный двойной слой
        Ubinn = nEqSysQ.GetStateFunction().GetIndepStateFunction().GetUbinn()  # Отрицательный двойной слой

        # Получаем напряжение мембраны
        Um = nEqSysQ.GetStateFunction().GetIndepStateFunction().GetUm()

        # Получаем напряжение внутри аккумулятора
        Uin = nEqSysQ.GetStateFunction().GetIndepStateFunction().GetUin()

        # Получаем ток во внешней цепи
        Icur = nEqSysQ.GetStateFunction().GetIndepStateFunction().GetIcur()

        # Получаем токи через двойные слои
        (Ibinp, Im, Ibinn) = nEqSysQ.GetVProcessCoordinates()[IElInd]

        # Выводим результат
        return (Ibinp, Im, Ibinn, Icur, Ubinp, Um, Ubinn, Uin)
    inds = np.arange(t.shape[0])  # Массив индексов
    GetAttrs = np.vectorize(GetElAttr)
    (Ibinp, Im, Ibinn, Icur, Ubinp, Um, Ubinn, Uin) = GetAttrs(inds)

    # Получаем координаты состояния
    q = stateCoordinates[:, qInd]  # Перенесенный через внешнюю цепь электричекий заряд
    qMatElp = stateCoordinates[:, qMatElpInd]  # Зарядовое число молей недеградированного положительного электрода
    qMatEln = stateCoordinates[:, qMatElnInd]  # Зарядовое число молей недеградированного отрицательного электрода
    qMatDegElp = stateCoordinates[:, qMatDegElpInd]  # Зарядовое число молей деградированного положительного электрода
    qMatDegEln = stateCoordinates[:, qMatDegElnInd]  # Зарядовое число молей деградированного отрицательного электрода
    qDegPosEl = stateCoordinates[:, qDegPosElInd]  # Зарядовое число молей разрушенного положительного электрода

    # Температура аккумулятора
    TInAkk = reducedTemp[:, TInAkkInd] - 273.15
    TBAkk = reducedTemp[:, TBAkkInd] - 273.15

    # Получаем параметры
    [Cnom,
     Rkl,
     Tokr,
     qMatAllp,
     qMatAlln
     ] = otherSystemParameters[otherSystemParametersInd]
    Tokr = np.full_like(Icur, Tokr - 273.15, dtype=np.double)  # Массив температур окружающей среды

    # Напряжение на клеммах (окончательно)
    Ukl = Uin - Icur * Rkl

    # Вывод результата
    return (t.reshape(-1,), Ukl, Ubinp, Ubinn, Um, TInAkk, TBAkk, q / Cnom,
            Ibinp * 3600 / Cnom, Im * 3600 / Cnom, Ibinn * 3600 / Cnom, Icur * 3600 / Cnom,
            Tokr, qMatElp / qMatAllp, qMatEln / qMatAlln,
            qMatDegElp / qMatAllp, qMatDegEln / qMatAlln, qDegPosEl / qMatAllp)
