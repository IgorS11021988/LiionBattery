import numpy as np


# Функция условий протекания процессов
def fICharge(t,  # Моменты времени
             systemParameters  # Параметры системы
             ):
    # Выделяем параметры динамики
    [I,  # Постоянная составляющая тока
     IStepDischarge,  # Ток разряда после ступеньки
     OmegaI,  # Частота колебаний внешнего тока
     AI,  # Амплитуда колебаний тока во внешней цепи
     OmegaOsI,  # Частота колебаний внешнего тока
     AOsI,  # Амплитуда колебаний тока во внешней цепи
     tIStep,  # Время начала конца заряда
     tauEndCharge,  # Постоянная времени конца разряда
     alphaEndCharge,  # Экспоненциальный коэффициент постоянной времени разряда
     powEndCharge  # Степень постоянной времени разряда
     ] = systemParameters[0:10]

    # Прочие параметры системы
    otherSystemParameters = systemParameters[10::]

    # Рассчитываем ток во внешней цепи
    Icur = np.full_like(t, I, dtype=np.double)  # Массив токов во внешней цепи
    bIcur = (t > tIStep)  # Моменты времени после ступеньки
    if np.any(bIcur):  # Есть моменты времени после ступеньки
        if I < 0:
            # Относительные моменты времени конца разряда
            rt = (t[bIcur] - tIStep) / tauEndCharge

            # Экспоненциальный коэффициент
            rLambda = 1 + alphaEndCharge * np.power(rt, powEndCharge)

            # Ток во внешней цепи
            Icur[bIcur] *= np.exp(-rLambda * rt)
        else:
            # Токи после ступеньки
            Icur[bIcur] = IStepDischarge
    OmegaI += AOsI * np.sin(OmegaOsI * t)  # Колебания частоты
    Icur += AI * np.sin(OmegaI * t)  # Учитываем колебания

    # Выводим результат
    return (Icur, otherSystemParameters)
