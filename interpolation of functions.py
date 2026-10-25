from math import log, sin, cos, pi

# --- Функція з варіанту ---
def f(x):
    return log(5*x) + (x**2)*sin(x)

# --- Рівномірні вузли ---
def uniform_nodes(a, b, n):
    step = (b - a) / n
    return [a + i*step for i in range(n+1)]

# --- Вузли Чебишова ---
def chebyshev_nodes(a, b, n):
    nodes = []
    for k in range(n+1):
        xk = 0.5*(a+b) - 0.5*(b-a)*cos((2*k+1)/(2*(n+1))*pi)
        nodes.append(xk)
    return nodes

# --- Лагранж ---
def lagrange(x, nodes, values):
    total = 0
    n = len(nodes)
    for i in range(n):
        term = values[i]
        for j in range(n):
            if i != j:
                term *= (x - nodes[j])/(nodes[i]-nodes[j])
        total += term
    return total

# --- Розділені різниці ---
def divided_diffs(nodes, values):
    n = len(nodes)
    table = [[0]*n for _ in range(n)]
    for i in range(n):
        table[i][0] = values[i]
    for j in range(1, n):
        for i in range(n-j):
            table[i][j] = (table[i+1][j-1] - table[i][j-1]) / (nodes[i+j] - nodes[i])
    return [table[0][j] for j in range(n)]

# --- Ньютон вперед ---
def newton_forward(x, nodes, coeffs):
    n = len(coeffs)
    result = coeffs[0]
    product = 1
    for i in range(1,n):
        product *= (x - nodes[i-1])
        result += coeffs[i]*product
    return result

# --- Ньютон назад ---
def newton_backward(x, nodes, coeffs):
    n = len(coeffs)
    result = coeffs[-1]
    product = 1
    for i in range(1,n):
        product *= (x - nodes[-i])
        result += coeffs[-(i+1)]*product
    return result

# --- Таблиця результатів ---
def print_table(title, nodes, approx_func):
    print("\n" + title)
    print(f"{'x':>10} {'f(x)':>15} {'P(x)':>15} {'|f-P|':>15}")
    for x in nodes:
        fx = f(x)
        px = approx_func(x)
        err = abs(fx-px)
        print(f"{x:10.6f} {fx:15.8f} {px:15.8f} {err:15.8e}")

# --- Основна програма ---
a, b = 1, 2
n = 4   # кількість відрізків (n+1 вузлів)

# Рівномірні вузли
nodes_u = uniform_nodes(a, b, n)
vals_u = [f(x) for x in nodes_u]

# Вузли Чебишова
nodes_c = chebyshev_nodes(a, b, n)
vals_c = [f(x) for x in nodes_c]

# Лагранж
print_table("Лагранж (рівномірні вузли)", nodes_u, lambda x: lagrange(x, nodes_u, vals_u))
print_table("Лагранж (Чебишов)", nodes_c, lambda x: lagrange(x, nodes_c, vals_c))

# Ньютон вперед
coeffs_u = divided_diffs(nodes_u, vals_u)
coeffs_c = divided_diffs(nodes_c, vals_c)

print_table("Ньютон вперед (рівномірні вузли)", nodes_u, lambda x: newton_forward(x, nodes_u, coeffs_u))
print_table("Ньютон вперед (Чебишов)", nodes_c, lambda x: newton_forward(x, nodes_c, coeffs_c))

# Ньютон назад
print_table("Ньютон назад (рівномірні вузли)", nodes_u, lambda x: newton_backward(x, nodes_u, coeffs_u))
print_table("Ньютон назад (Чебишов)", nodes_c, lambda x: newton_backward(x, nodes_c, coeffs_c))