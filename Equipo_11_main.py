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
    for line in lines[1:]:
        graph.append(list(map(int, line.strip().split())))
    
    return n, graph

filename = input("Nombre del archivo: ")
n, graph = leer_archivo(filename)
path, cost = repeated_nearest_neighbor(graph, n)
print("Mejor ruta encontrada:", path)
print("Costo de la mejor ruta:", cost)
