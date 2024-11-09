"""
    Problema 2.
"""
def nearest_neighbor(graph, start, n):
    """
    Encuentra una ruta aproximada de menor costo usando el algoritmo de vecino más cercano.

    Parámetros:
        graph (list): Matriz de adyacencia que representa el costo entre nodos.
        start (int): Nodo inicial desde el cual comenzar la ruta.
        n (int): Número total de nodos en el grafo.

    Returns:
        list: Lista que representa el camino recorrido, terminando en el nodo inicial.
    """
    visited = [False] * n
    path = [start]
    visited[start] = True
    current = start

    for _ in range(n - 1):
        next_city = None
        min_dist = float('inf')
        for city in range(n):
            if not visited[city] and graph[current][city] < min_dist:
                min_dist = graph[current][city]
                next_city = city
        path.append(next_city)
        visited[next_city] = True
        current = next_city

    path.append(start)
    return path

def repeated_nearest_neighbor(graph, n):
    """
    Encuentra la mejor ruta de menor costo utilizando el algoritmo de vecino más cercano repetido 
    desde cada nodo como punto de partida.

    Parámetros:
        graph (list): Matriz de adyacencia que representa el costo entre nodos.
        n (int): Número total de nodos en el grafo.

    Returns:
        tuple: Contiene la mejor ruta encontrada (list) y el costo de la mejor ruta (int).
    """
    best_path = None
    best_cost = float('inf')

    for start in range(n):
        path = nearest_neighbor(graph, start, n)
        cost = sum(graph[path[i]][path[i + 1]] for i in range(n))
        if cost < best_cost:
            best_cost = cost
            best_path = path

    return best_path, best_cost


"""
    Problema 3.
"""
def bfs(graph, s, t, parent, n):
    """
    Función: bfs
    Parámetros: graph, s, t, parent, n
            graph: Matriz de flujos
            s: Nodo de inicio
            t: Nodo de destino
            parent: Lista de nodos padres
            n: Número de nodos
    Regresa: True si nodo de destino es alcanzable desde el nodo de inicio,
            False en caso contrario
    Descripción: Realiza un recorrido en anchura del grafo para encontrar un camino
                desde el nodo de inicio hasta el nodo de destino.
    """
    visited = [False] * n
    queue = []
    queue.append(s)
    visited[s] = True

    while queue:
        u = queue.pop(0)

        for i in range(len(graph[u])):
            val =  graph[u][i]
            if visited[i] == False and val > 0:
                queue.append(i)
                visited[i] = True
                parent[i] = u

    if visited[t]:
        return True
    else:
        return False


def ford_fulkerson(graph, inicio, destino, n):
    """
    Función: ford_fulkerson
    Parámetros: graph, inicio, destino, n
            graph: Matriz de flujos
            inicio: Nodo de inicio
            destino: Nodo de destino
            n: Número de nodos
    """
    # Lista de nodos padres
    parent = [-1] * n
    # Flujo máximo
    max_flow = 0

    # Mientras haya caminos existentes de nodo inicio a destino
    while bfs(graph, inicio, destino, parent, n):
        # Flujo del camino
        path_flow = float("Inf")
        actual = destino

        # Recorre desde nodo destino hasta nodo inicio
        while actual != inicio:
            # Mínimo entre flujo del camino y flujo de trayecto actual
            path_flow = min(path_flow, graph[parent[actual]][actual])
            actual = parent[actual]

        v = destino
        # Actualizar flujo de los nodos del camino
        while v != inicio:
            u = parent[v]
            graph[u][v] -= path_flow
            graph[v][u] += path_flow
            v = parent[v]

        # Acumular flujo máximo de todos los caminos
        max_flow += path_flow

    return max_flow

    

def leer_archivo(filename):
    """
    Lee un archivo de texto y convierte su contenido en una matriz de adyacencia para un grafo.

    Parámetros:
        filename (str): Nombre del archivo que contiene el número de nodos en la primera línea 
                        y la matriz de adyacencia en las siguientes líneas.

    Returns:
        tuple: Contiene el número de nodos (int) y la matriz de adyacencia del grafo (list).
    """
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    n = int(lines[0].strip())
    graph = []
    flujos = []

    for i in range(1, n + 1):
        graph.append(list(map(int, lines[i].strip().split())))
    
    for i in range(n + 1, 2 * n + 1):
        flujos.append(list(map(int, lines[i].strip().split())))
    
    return n, graph, flujos

filename = input("Nombre del archivo: ")
n, graph, flujos = leer_archivo(filename)

path, cost = repeated_nearest_neighbor(graph, n)
max_flow = ford_fulkerson(flujos, 0, n - 1, n)

print("Mejor ruta encontrada:", path)
print("Costo de la mejor ruta:", cost)
print("Flujo máximo:", max_flow)
