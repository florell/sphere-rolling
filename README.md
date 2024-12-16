# Задача

Смоделировать движение шара в перевёрнутой полусферической чаше в условиях трения, визуализировать это в 3D и обеспечить управление параметрами и скоростью симуляции. При этом требуется:

1. Изначально шар располагается близко к краю перевёрнутой полусферы.
2. Под действием гравитации шар скатывается вниз.
3. Модель учитывает вращение шара без проскальзывания, но для стабилизации и конечной остановки вводится трение (демпфирование).
4. Должна быть предусмотрена возможность изменения параметров системы (массы шара, радиусов, ускорения свободного падения и др.) и скорости проигрывания симуляции.
5. По достижении шаром дна и остановки симуляция автоматически прекращается, при этом фиксируются максимальная достигнутая скорость и время её достижения.

# Формализация

- **Геометрия:**
  - Полусфера радиуса $\( R \)$ перевёрнута так, что её край лежит в плоскости $\( z = 0 \)$, а центр (нижняя точка чаши) в точке $\( z = -R \)$.
  - Положение шара задаётся углами $\(\theta\)$ (от вертикали вниз) и $\(\varphi\)$ (азимутальный угол в горизонтальной плоскости).
  - Координаты центра шара:
    
    $x = (R - r_{\text{ball}})\sin\theta\cos\varphi,\quad$
    $y = (R - r_{\text{ball}})\sin\theta\sin\varphi,\quad$ 
    $z = -(R - r_{\text{ball}})\cos\theta$.
  
- **Динамика без трения:**
  - Потенциальная энергия: $\( V = mg(-R_{\text{eff}}\cos\theta) \)$, где $\( R_{\text{eff}} = R - r_{\text{ball}} \)$.
  - Кинетическая энергия для качения без проскальзывания: $\( T = \tfrac{7}{10} m R_{\text{eff}}^2 (\dot{\theta^2} + \sin^2\theta \dot{\varphi^2}) \)$.
  - Уравнение для $\(\theta\)$:
    $\ddot{\theta} = \sin\theta\cos\theta \dot{\varphi^2} - \frac{5}{7}\frac{g}{R_{\text{eff}}}\sin\theta$.

- **Добавление трения (демпфирования):**
  - Демпфирующий член для угловой скорости $\(\dot{\theta}\): \(-\gamma_{\theta} \dot{\theta}\)$.
  - Демпфирование по $\(\varphi\)$: $\(\dot{\varphi} \leftarrow \dot{\varphi}(1 - \gamma_{\varphi} \Delta t)\)$ для экспоненциального затухания.
  - Итоговое уравнение для $\(\theta\)$:
    $\ddot{\theta} = \sin\theta\cos\theta \dot{\varphi^2} - \frac{5}{7}\frac{g}{R_{\text{eff}}}\sin\theta - \gamma_{\theta}\dot{\theta}$.

- **Численное решение:**
  - Используется численный метод (например, метод Рунге-Кутты второго порядка) для интегрирования уравнений по $\(\theta(t)\)$ и простое обновление $\(\varphi(t)\)$.
  - На каждом шаге вычисляются текущая скорость шара:
    $v = R_{\text{eff}}\sqrt{\dot{\theta^2} + \sin^2\theta \dot{\varphi^2}}$.

- **Условия остановки:**
  - Если $\(\theta \approx 0\)$ (шар у дна чаши) и $\(\dot{\theta}, \dot{\varphi}\)$ близки к нулю, считаем шар остановившимся.
  - Фиксируем достигнутую максимальную скорость и время её достижения, затем останавливаем симуляцию.

# Этапы выполнения

1. **Постановка задачи и выбор координат:**
   - Определить систему координат и параметризацию положения шара на внутренней поверхности полусферы.
   - Задать граничные условия и начальное положение.

2. **Физическая модель без трения:**
   - Записать выражения для энергии.
   - Определить дифференциальные уравнения движения $\(\theta(t)\)$ и $\(\varphi(t)\)$.

3. **Введение трения:**
   - Добавить в уравнения движения демпфирующие члены, обеспечивающие постепенное затухание движения.

4. **Численное решение:**
   - Выбрать шаг интегрирования $\(\Delta t\)$.
   - Реализовать численный метод интегрирования (например, РК2).

5. **Реализация интерфейса (UI):**
   - Создать окно с настройками параметров (масса, радиусы, $g$, начальные условия).
   - Добавить контроль скорости проигрывания симуляции.
   - Добавить кнопку старта, паузы и сброса.

6. **3D-визуализация:**
   - Отобразить полусферу и шар в 3D.
   - На каждом шаге анимации обновлять положение шара.

7. **Автоматическая остановка и фиксация результатов:**
   - При достижении шаром состояния покоя внизу чаши остановить симуляцию.
   - Запомнить максимальную скорость и время её достижения, отобразить эту информацию пользователю.

8. **Тестирование и отладка:**
   - Проверить работу для различных параметров.
   - Убедиться, что при включенном трении шар в итоге останавливается.