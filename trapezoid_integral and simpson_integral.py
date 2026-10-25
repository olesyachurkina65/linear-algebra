import math

# Підінтегральна функція f(x)
def f(x):
    return math.log(1 + x) / (1 + x * x)

# Метод трапецій із подвоєнням кроку та оцінкою похибки за Рунге
def trapezoid_integral(a, b, eps):
    n = 1
    h = b - a
    I_prev = h * (f(a) + f(b)) / 2   # перше наближення

    while True:
        n *= 2
        h /= 2

        # сума значень у нових середніх точках
        s = sum(f(a + i * h) for i in range(1, n, 2))

        # нове наближення
        I_curr = I_prev / 2 + h * s

        # похибка Рунге для методу трапецій
        r = (I_curr - I_prev) / 3

        if abs(r) <= eps:
            return I_curr + r

        I_prev = I_curr

# Метод парабол (Сімпсона) з подвоєнням кроку та оцінкою похибки за Рунге
def simpson_integral(a, b, eps):
    n = 1
    h = (b - a) / 2

    s0 = f(a) + f(b)     # значення на кінцях
    s1 = f(a + h)        # середина
    s2 = 0               # парні вузли
    I_prev = (h / 3) * (s0 + 4 * s1 + 2 * s2)

    while True:
        n *= 2
        h /= 2

        # сума для парних вузлів
        s2 += s1

        # сума для непарних вузлів
        s1 = sum(f(a + i * h) for i in range(1, 2 * n, 2))

        # нове наближення
        I_curr = (h / 3) * (s0 + 4 * s1 + 2 * s2)

        # похибка Рунге для Сімпсона
        r = (I_curr - I_prev) / 15

        if abs(r) <= eps:
            return I_curr + r

        I_prev = I_curr

# Вхідні дані
eps = 1e-6
a = 0
b = 1

# Обчислення
I_trap = trapezoid_integral(a, b, eps)
I_simp = simpson_integral(a, b, eps)

# Вивід результатів
print("Інтеграл (метод трапецій):", I_trap)
print("Інтеграл (метод парабол / Сімпсона):", I_simp)
