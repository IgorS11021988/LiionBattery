import numpy as np

from .fICharge import fICharge


# Функция условий протекания процессов
def fU(t,  # Моменты времени
       systemParameters  # Параметры системы
       ):
    # Выделяем параметры динамики и отдельно свойства веществ и процессов
    (Icur, otherSystemParameters) = fICharge(np.array([t], dtype=np.double),  # Моменты времени
                                             systemParameters  # Параметры системы
                                             )

    # Выводим результат
    return np.hstack((Icur, otherSystemParameters))
