import numpy as np

import .LiionBatteryFunctionsStation as lbfs
import .fICharge as fic


# Функция состояния для литий-ионного аккумулятора
def LiionBatteryCharacteristicsFunction(t,  # Моменты времени
                                        stateCoordinates,  # Координаты состояния
                                        reducedTemp,  # Приведенные температуры
                                        systemParameters  # Параметры системы
                                        ):
    # Получаем динамику тока
    (Icur, otherSystemParameters) = fic.fICharge(np.array(t, dtype=np.double),  # Моменты времени
                                                 systemParameters  # Параметры системы
                                                 )
    Icur = np.array(Icur, dtype=np.double).reshape(-1)  # Приводим токи к одномерному массиву

    # Получаем координаты состояния
    qbinp = stateCoordinates[:, 0]  # Заряд на положительном двойном слое
    qm = stateCoordinates[:, 1]  # Заряд на мембране
    qbinn = stateCoordinates[:, 2]  # Заряд на отрицательном двойном слое
    q = stateCoordinates[:, 3]  # Перенесенный через внешнюю цепь электричекий заряд

    # Температура содержимого аккумулятора
    TInAkk = reducedTemp[:, 0] - 273.15
    TBAkk = reducedTemp[:, 1] - 273.15

    # Температура окружающей среды
    Tokr = otherSystemParameters[0] - 273.15
    Tokr = np.full_like(Icur, Tokr, dtype=np.double)  # Массив температур окружающей среды

    # Получаем параметры
    EbinpC = otherSystemParameters[1]  # ЭДС положительного двойного слоя в заряженном состоянии
    EbinnC = otherSystemParameters[2]  # ЭДС отрицательного двойного слоя в заряженном состоянии
    EbinpD = otherSystemParameters[3]  # ЭДС положительного двойного слоя в разряженном состоянии
    EbinnD = otherSystemParameters[4]  # ЭДС отрицательного двойного слоя в разряженном состоянии
    Cbin0p = otherSystemParameters[5]  # Емкость положительного двойного слоя
    Cm = otherSystemParameters[6]  # Емкость мембраны
    Cbin0n = otherSystemParameters[7]  # Емкость отрицательного двойного слоя
    Cnom = otherSystemParameters[15]  # Номинальная емкость литий-ионного аккумулятора
    rLiEpE = otherSystemParameters[16]  # Приведенное зарядовое число положительного электрода (по ЭДС)
    rLiEnE = otherSystemParameters[17]  # Приведенное зарядовое число отрицательного электрода (по ЭДС)
    alphaCQp = otherSystemParameters[33]  # Зарядовый коэффициент емкости положительного электрода, 1/Кл
    alphaCQn = otherSystemParameters[34]  # Зарядовый коэффициент емкости отрицательного электрода, 1/Кл

    # Получаем довесочные коэффициенты
    betaCQ2p = otherSystemParameters[49]
    betaCQ2n = otherSystemParameters[50]
    betaCQ3p = otherSystemParameters[51]
    betaCQ3n = otherSystemParameters[52]
    betaEQ2p = otherSystemParameters[53]
    betaEQ2n = otherSystemParameters[54]
    betaEQ3p = otherSystemParameters[55]
    betaEQ3n = otherSystemParameters[56]

    # Получаем сопротивление клемм
    Rkl = otherSystemParameters[-1]

    # Рассчитываем приведенные (в единицах зарядовой емкости) числа молей интеркалированных в электроды ионов лития
    nuLip = qbinp + q  # Положительный электрод
    nuLin = qbinn + q  # Отрицательный электрод

    # Определяем ЭДС двойных слоев
    (Ebinp, Ebinn) = lbfs.funEbin(EbinpC, EbinnC, EbinpD, EbinnD, nuLip,
                                  nuLin, rLiEpE, rLiEnE, Cnom,
                                  betaEQ2p, betaEQ2n, betaEQ3p, betaEQ3n)

    # Определяем емкости двойных слоев
    (Cbinp, Cbinn) = lbfs.funCbin(qbinp, qbinn, alphaCQp, alphaCQn, Cbin0p, Cbin0n,
                                  betaCQ2p, betaCQ2n, betaCQ3p, betaCQ3n)

    # Рассчитываем напряжения двойных слоев
    Ubinp = qbinp / Cbinp  # Положительный двойной слой
    Um = qm / Cm  # Мембрана
    Ubinn = qbinn / Cbinn  # Отрицательный двойной слой

    # Напряжение на клеммах
    Ukl = Ubinp + Um + Ubinn - Icur * Rkl
    return (t, Ukl, Ubinp, Ubinn, Um, TInAkk, TBAkk, q / Cnom, Icur * 3600 / Cnom, Tokr)
