import numpy as np

from .StationFunctions import funCbin
from .fICharge import fICharge


# Функция состояния для литий-ионного аккумулятора
def CharacteristicsFunction(t,  # Моменты времени
                            stateCoordinates,  # Координаты состояния
                            reducedTemp,  # Приведенные температуры
                            systemParameters  # Параметры системы
                            ):
    # Получаем динамику тока
    (Icur, otherSystemParameters) = fICharge(np.array(t, dtype=np.double),  # Моменты времени
                                             systemParameters  # Параметры системы
                                             )
    Icur = np.array(Icur, dtype=np.double).reshape(-1)  # Приводим токи к одномерному массиву

    # Получаем координаты состояния
    qbinp = stateCoordinates[:, 0]  # Заряд на положительном двойном слое
    qm = stateCoordinates[:, 1]  # Заряд на мембране
    qbinn = stateCoordinates[:, 2]  # Заряд на отрицательном двойном слое
    q = stateCoordinates[:, 3]  # Перенесенный через внешнюю цепь электричекий заряд
    qMatElp = stateCoordinates[:, 4]  # Зарядовое сило молей недеградированного положительного электрода
    qMatEln = stateCoordinates[:, 5]  # Зарядовое сило молей недеградированного отрицательного электрода
    qMatDegElp = stateCoordinates[:, 6]  # Зарядовое сило молей деградированного положительного электрода
    qMatDegEln = stateCoordinates[:, 7]  # Зарядовое сило молей деградированного отрицательного электрода
    qDegPosEl = stateCoordinates[:, 8]  # Зарядовое сило молей разрушенного положительного электрода

    # Температура содержимого аккумулятора
    TInAkk = reducedTemp[:, 0] - 273.15
    TBAkk = reducedTemp[:, 1] - 273.15

    # Температура окружающей среды
    Tokr = otherSystemParameters[0] - 273.15
    Tokr = np.full_like(Icur, Tokr, dtype=np.double)  # Массив температур окружающей среды

    # Получаем параметры
    Cbin0p = otherSystemParameters[5]  # Емкость положительного двойного слоя
    Cm = otherSystemParameters[6]  # Емкость мембраны
    Cbin0n = otherSystemParameters[7]  # Емкость отрицательного двойного слоя
    Cnom = otherSystemParameters[15]  # Номинальная емкость литий-ионного аккумулятора
    alphaCQp = otherSystemParameters[33]  # Зарядовый коэффициент емкости положительного электрода
    alphaCQn = otherSystemParameters[34]  # Зарядовый коэффициент емкости отрицательного электрода
    qMatAllp = otherSystemParameters[35]  # Общее зарядовое число молей материала положительного электрода
    qMatAlln = otherSystemParameters[36]  # Общее зарядовое число молей материала отрицательного электрода

    # Получаем довесочные коэффициенты
    betaCQ2p = otherSystemParameters[75]
    betaCQ2n = otherSystemParameters[76]
    betaCQ3p = otherSystemParameters[77]
    betaCQ3n = otherSystemParameters[78]

    # Получаем сопротивление клемм
    Rkl = otherSystemParameters[-1]

    # Определяем емкости двойных слоев
    (Cbinp, Cbinn) = funCbin(qbinp, qbinn, alphaCQp, alphaCQn, Cbin0p, Cbin0n,
                             betaCQ2p, betaCQ2n, betaCQ3p, betaCQ3n)

    # Рассчитываем напряжения двойных слоев
    Ubinp = qbinp / Cbinp  # Положительный двойной слой
    Um = qm / Cm  # Мембрана
    Ubinn = qbinn / Cbinn  # Отрицательный двойной слой

    # Напряжение на клеммах
    Ukl = Ubinp + Um + Ubinn - Icur * Rkl
    return (t, Ukl, Ubinp, Ubinn, Um, TInAkk, TBAkk, q / Cnom, Icur * 3600 / Cnom, Tokr,
            qMatElp / qMatAllp, qMatEln / qMatAlln,
            qMatDegElp / qMatAllp, qMatDegEln / qMatAlln, qDegPosEl / qMatAllp)
