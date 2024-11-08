from distutils.core import setup

#Функция установки
setup(name = "LiionBattery",
      version = "1.0",
      author = "Igor Starostin",
      author_email = "starostinigo@yandex.ru",
      description="Liion accumulators modeling by mathematical prototiping method",
      packages = ["LiionBattery",
                  "LiionBattery.LiionBatteryDynamics",
                  "LiionBattery.LiionBatteryFunctions"],
      scripts = [],
      dependency_links = ["numpy","pandas","matplotlib"]
      )