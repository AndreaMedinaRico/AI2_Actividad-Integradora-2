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
        path (list): Lista que representa el camino recorrido, terminando en el nodo inicial.
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
        tuple: best_path(lis) Mejor ruta encontrada 
                best_cost(int) Costo de la mejor ruta
    """
    best_path = None
    best_cost = float('inf')

    for start in range(n):
        path = nearest_neighbor(graph, start, n)
        cost = sum(graph[path[i]][path[i + 1]] for i in range(n))
        if cost < best_cost:
            best_cost = cost
            best_path = path

    index_to_letter = {i: chr(65 + i) for i in range(n)}
    best_path = [index_to_letter[i] for i in best_path]

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
    Returns: max_flow(int) Flujo máximo
    Descripción: Algoritmo de fold-fulkerson para encontrar el flujo máximo
                de un grafo.
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



"""
    Problema 4.
"""
    
import math

class NodoKDTree:
    def __init__(self, punto, izquierda = None, derecha = None):
        self.punto = punto          # Coordenada del nodo
        self.izquierda = izquierda  # Subárbol izquierdo
        self.derecha = derecha      # Subárbol derecho

def construir_kdtree(puntos, profundidad = 0):
    """
    Función: construir_kdtree
    Parámetros: 
            puntos (list): Lista de puntos en el plano
            profundidad(int): Nivel de profundidad del árbol
    Returns: NodoKDTree
    Descripción: Construye un árbol KD a partir de una lista de puntos en el plano.        
    """
    if not puntos:
        return None

    # Elegir el eje (x o y) según la profundidad
    k = 2     # Coordenadas -> 2
    eje = profundidad % k   # Coordenada 0 -> x

    # Ordena los puntos y elige punto medio
    puntos.sort(key = lambda punto: punto[eje])
    mediana = len(puntos) // 2

    # Crea nodo y construye subárboles recursivamente
    return NodoKDTree(
        punto = puntos[mediana],
        izquierda = construir_kdtree(puntos[:mediana], profundidad + 1),
        derecha = construir_kdtree(puntos[mediana + 1:], profundidad + 1)
    )

def calcular_distancia(punto1, punto2):
    """
    Función: distancia_euclidiana
    Parámetros: 
            punto1 (tuple): Coordenadas del primer punto
            punto2 (tuple): Coordenadas del segundo punto
    Returns: float de la distancia
    Descripción: Calcula la distancia eucladiana entre dos puntos.
    """
    suma = 0
    for i in range(len(punto1)):
        suma += (punto1[i] - punto2[i]) ** 2

    return math.sqrt(suma)

def busqueda_kdtree(nodo, punto, profundidad = 0, mejor = None):
    """
    Función: busqueda_kdtree
    Parámetros:
            nodo (NodoKDTree): Nodo actual
            punto (tuple): Coordenadas del punto a buscar
            profundidad (int): Nivel de profundidad del árbol
            mejor (tuple): Mejor punto encontrado
    Returns: tuple del punto más cercano
    Descripción: Busca el punto más cercano a un punto dado en un árbol KD
    """
    if nodo is None:
        return mejor

    # Actualiza la mejor central
    if (mejor is None) or (calcular_distancia(punto, nodo.punto) < calcular_distancia(punto, mejor)):
        mejor = nodo.punto

    k = 2
    eje = profundidad % k

    # Decide subárbol para búsqueda
    if punto[eje] < nodo.punto[eje]:
        siguiente_rama = nodo.izquierda 
        otra_rama = nodo.derecha 
    else:
        siguiente_rama = nodo.derecha
        otra_rama = nodo.izquierda

    # Búsqueda recursiva en subárbol más probable
    mejor = busqueda_kdtree(siguiente_rama, punto, profundidad + 1, mejor)

    # Si es necesario, buscar en el otro subárbol
    if abs(punto[eje] - nodo.punto[eje]) < calcular_distancia(punto, mejor):
        mejor = busqueda_kdtree(otra_rama, punto, profundidad + 1, mejor)

    return mejor


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
    coordenadas = []

    for i in range(1, n + 1):
        graph.append(list(map(int, lines[i].strip().split())))
    
    for i in range(n + 1, 2 * n + 1):
        flujos.append(list(map(int, lines[i].strip().split())))

    for i in range(2 * n + 1, len(lines) - 1):
        coordenadas.append(tuple(map(int, lines[i].strip().strip('()').split(','))))

    nueva_central = tuple(map(int, lines[-1].strip().strip('()').split(',')))
    
    return n, graph, flujos, coordenadas, nueva_central

filename = input("Nombre del archivo: ")
n, graph, flujos, coordenadas, nueva_central = leer_archivo(filename)

path, cost = repeated_nearest_neighbor(graph, n)
max_flow = ford_fulkerson(flujos, 0, n - 1, n)

arbol = construir_kdtree(coordenadas)
central_cercana = busqueda_kdtree(arbol, nueva_central)
distancia = calcular_distancia(nueva_central, central_cercana)

print("\nParte 2:\nMejor ruta encontrada:", path)
print("Costo de la mejor ruta:", cost)
print("\nParte 3:\nFlujo máximo:", max_flow)
print("\nParte 4:\nCentral más cercana:", central_cercana, "\nDistancia: ", distancia, "km")
