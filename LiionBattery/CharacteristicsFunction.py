import numpy as np

from MathProtEnergyProcBase.IndexFunctions import GetIndex, GetIndexes

from .AttributesNames import stateCoordinatesNames, reducedTemperaturesEnergyPowersNames, USystemParametersNames, otherSystemParametersNames
from .StationFunctions import funCbin


# Индексы координат состояния
qbinpInd = GetIndex(stateCoordinatesNames, "qbinp")  # Индекс заряда положительного электрода
qmInd = GetIndex(stateCoordinatesNames, "qm")  # Индекс заряда мембраны
qbinnInd = GetIndex(stateCoordinatesNames, "qbinn")  # Индекс заряда отрицательного электрода
qInd = GetIndex(stateCoordinatesNames, "q")  # Индекс перенесенного через внешнюю цепь электрического заряда
qMatElpInd = GetIndex(stateCoordinatesNames, "qMatElp")  # Индекс зарядового числа молей недеградированного положительного электрода
qMatElnInd = GetIndex(stateCoordinatesNames, "qMatEln")  # Индекс зарядового числа молей недеградированного отрицательного электрода
qMatDegElpInd = GetIndex(stateCoordinatesNames, "qMatDegElp")  # Индекс зарядового числа молей деградированного положительного электрода
qMatDegElnInd = GetIndex(stateCoordinatesNames, "qMatDegEln")  # Индекс зарядового числа молей деградированного отрицательного электрода
qDegPosElInd = GetIndex(stateCoordinatesNames, "qDegPosEl")  # Индекс зарядового числа молей разрушенного положительного электрода

# Индексы приведенных температур
TInAkkInd = GetIndex(reducedTemperaturesEnergyPowersNames, "TInAkk")  # Индекс температуры содержимого аккумулятора
TBAkkInd = GetIndex(reducedTemperaturesEnergyPowersNames, "TBAkk")  # Индекс температуры корпуса аккумулятора

# Индексы переменных параметров системы
IInd = GetIndex(USystemParametersNames, "I")  # Индекс тока

# Индексы параметров системы
systemParametersIndexes = GetIndexes(otherSystemParametersNames, ["Tokr",  # Температура окружающей среды
                                                                  "Cbin0p",  # Емкость положительного двойного слоя
                                                                  "Cm",  # Емкость мембраны
                                                                  "Cbin0n",  # Емкость отрицательного двойного слоя
                                                                  "Cnom",  # Номинальная емкость
                                                                  "alphaCQp",  # Зарядовый коэффициент емкости положительного электрода
                                                                  "alphaCQn",  # Зарядовый коэффициент емкости отрицательного электрода
                                                                  "qMatAllp",  # Общее зарядовое число молей материала положительного электрода
                                                                  "qMatAlln",  # Общее зарядовое число молей материала отрицательного электрода

                                                                  "betaCQ2p",
                                                                  "betaCQ2n",
                                                                  "betaCQ3p",
                                                                  "betaCQ3n",

                                                                  "Rkl"  # Сопротивление клемм
                                                                  ])


# Функция состояния для литий-ионного аккумулятора
def CharacteristicsFunction(t,  # Моменты времени
                            stateCoordinates,  # Координаты состояния
                            reducedTemp,  # Приведенные температуры
                            USystemParameters,  # U-параметры системы
                            otherSystemParameters  # Прочие параметры системы
                            ):
    # Получаем динамику тока
    Icur = USystemParameters[:, IInd]  # Ток в текущие моменты времени

    # Получаем координаты состояния
    qbinp = stateCoordinates[:, qbinpInd]  # Заряд на положительном двойном слое
    qm = stateCoordinates[:, qmInd]  # Заряд на мембране
    qbinn = stateCoordinates[:, qbinnInd]  # Заряд на отрицательном двойном слое
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
    [Tokr,  # Температура окружающей среды
     Cbin0p,  # Емкость положительного двойного слоя
     Cm,  # Емкость мембраны
     Cbin0n,  # Емкость отрицательного двойного слоя
     Cnom,  # Номинальная емкость литий-ионного аккумулятора
     alphaCQp,  # Зарядовый коэффициент емкости положительного электрода
     alphaCQn,  # Зарядовый коэффициент емкости отрицательного электрода
     qMatAllp,  # Общее зарядовое число молей материала положительного электрода
     qMatAlln,  # Общее зарядовое число молей материала отрицательного электрода

     # Получаем довесочные коэффициенты
     betaCQ2p,
     betaCQ2n,
     betaCQ3p,
     betaCQ3n,

     # Получаем сопротивление клемм
     Rkl] = otherSystemParameters[systemParametersIndexes]
    Tokr = np.full_like(Icur, Tokr - 273.15, dtype=np.double)  # Массив температур окружающей среды

    # Определяем емкости двойных слоев
    (Cbinp, Cbinn) = funCbin(qbinp, qbinn, alphaCQp, alphaCQn, Cbin0p, Cbin0n,
                             betaCQ2p, betaCQ2n, betaCQ3p, betaCQ3n)

    # Рассчитываем напряжения двойных слоев
    Ubinp = qbinp / Cbinp  # Положительный двойной слой
    Um = qm / Cm  # Мембрана
    Ubinn = qbinn / Cbinn  # Отрицательный двойной слой

    # Напряжение на клеммах
    Ukl = Ubinp + Um + Ubinn - Icur * Rkl
    return (t.reshape(-1,), Ukl, Ubinp, Ubinn, Um, TInAkk, TBAkk, q / Cnom, Icur * 3600 / Cnom,
            Tokr, qMatElp / qMatAllp, qMatEln / qMatAlln,
            qMatDegElp / qMatAllp, qMatDegEln / qMatAlln, qDegPosEl / qMatAllp)
