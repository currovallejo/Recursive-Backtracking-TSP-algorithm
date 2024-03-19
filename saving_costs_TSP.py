# Imports
from math import atan2, cos, sin, sqrt, radians
from matplotlib import pyplot as plt
import matplotlib.animation

# Problem Data
cities_andalucia = {'Almeria': (36.83892362, -2.46413188),
                    'Cadiz': (36.52171152, -6.28414575),
                    'Cordoba': (37.87954225, -4.78032455), 
                    'Granada':(37.17641932, -3.60001883),
                    'Huelva': (37.26004113, -6.95040588),
                    'Jaen': (37.7651913, -3.790359),
                    'Málaga': (36.72034267, -4.41997511),
                    'Sevilla': (37.38620512, -5.99251368)}

coord_cities = {'Huelva': (37.26638, -6.94004),
                'Barcelona': (41.3887900, 2.1589900),
                'Bilbao': (43.25721957, -2.92390606), 
                'La Coruña':(43.37012643, -8.39114853),
                'Lisboa': (38.71667, -9.13333),
                'Madrid': (40.40841191, -3.68760088),
                'Murcia': (37.98436361, -1.1285408),
                'Salamanca': (40.96736822, -5.66538084),
                'Sevilla': (37.38620512, -5.99251368),
                'Valencia': (39.47534441, -0.37565717),
                'Zaragoza': (41.65645655, -0.87928652),
                'Oporto': (41.14961, -8.61099),
                'Cáceres': (39.47316762, -6.37121092),
                'Toulouse': (43.60426, 1.44367),
                'León': (42.59912097, -5.56707631), 
                'Ciudad Real': (38.98651781, -3.93131981)
            }

# Functions
def distance_earth(a, b):
    '''
    The great circle
    '''
    EARTH_RADIUS = 6370

    lat1, lng1 = radians(a[0]), radians(a[1])
    lat2, lng2 = radians(b[0]), radians(b[1])

    sin_lat1, cos_lat1 = sin(lat1), cos(lat1)
    sin_lat2, cos_lat2 = sin(lat2), cos(lat2)

    delta_lng = lng2 - lng1
    cos_delta_lng, sin_delta_lng = cos(delta_lng), sin(delta_lng)

    d = atan2(sqrt((cos_lat2 * sin_delta_lng) ** 2 + (cos_lat1 * sin_lat2 - sin_lat1 * cos_lat2 * cos_delta_lng) ** 2),
              sin_lat1 * sin_lat2 + cos_lat1 * cos_lat2 * cos_delta_lng)

    return EARTH_RADIUS * d

def plot_problem(dic_cities: dict):
    cities = list(dic_cities.keys()) # List of cities to visit
    X = [dic_cities[city][1] for city in cities] # List of X coordenates
    Y = [dic_cities[city][0] for city in cities] # List of Y coordenates
    X.append(X[0]) # Initial position X 
    Y.append(Y[0]) # Initial position Y
    plt.plot(X,Y,'o')
    for city in cities:
        plt.annotate(city, (dic_cities[city][1]+0.01, dic_cities[city][0]))
        
def plot_solution(solution: list, dic_cities: dict):
    X = [dic_cities[city][1] for city in solution]
    Y = [dic_cities[city][0] for city in solution]
    X.append(X[0])
    Y.append(Y[0])
    plt.plot(X,Y)

def animate_solution_arcs(cities: list, arcs, cities_dict: dict):
    plt.rcParams["animation.html"] = "jshtml"
    plt.rcParams['figure.dpi'] = 150  
    # Initialize the figure and axis
    plt.ioff()
    fig, ax = plt.subplots()
    X = [cities_dict[city][1] for city in cities] # List of X coordenates
    Y = [cities_dict[city][0] for city in cities] # List of Y coordenates
    X.append(X[0]) # Initial position X 
    Y.append(Y[0]) # Initial position Y

    def animate(frame):
        plt.cla()
        plt.plot(X,Y,'o')
        for city in cities:
            plt.annotate(city, (cities_dict[city][1]+0.01, cities_dict[city][0]))
        for arc in arcs[:frame]:
            node1 = cities_dict[arc[1]]
            node2 = cities_dict[arc[0]]
            plt.plot([node1[1], node2[1]], [node1[0], node2[0]], 'r')  # Plotting arcs
        
        plt.xlim(min(X)-1,max(X)+1)
        plt.ylim(min(Y)-1,max(Y)+1)

    return matplotlib.animation.FuncAnimation(fig, animate, frames=len(arcs)+1)

# Saving cost algorithm 
def saving_cost_tsp(coord_cities:dict,func_distance,start=None):

    print("number of nodes ", len(coord_cities))
    
    def distances_between_cities():
        distancesCities={
            (city1,city2):func_distance(coord_cities[city1],coord_cities[city2])
            for i,city1 in enumerate(coord_cities.keys())
            for j,city2 in enumerate(coord_cities.keys())
            if (i!=j and j>i)
        }
        print("number of arcs: ", len(distancesCities))

        return distancesCities
    
    def order_arcs(distancesCities):
        sorted_distances = {cities: distance for cities, distance in sorted(distancesCities.items(), key=lambda x: x[1])}

        # for cities, distance in sorted_distances.items():
        #     print(f"Distance between {cities[0]} and {cities[1]}: {distance}")

        return sorted_distances
    
    def has_overconnected_city(arcs):
        visited = []
        for arc in arcs:
            for node in arc:
                visited.append(node)
                if visited.count(node)>2:
                    print("La ciudad ", visited[-1], "está sobreconectada")
                    arcs.pop()
                    return True

                #print("visitadas arcos", visited)
                #print(node, "se ha visitado ", visited.count(node))
        return False
    
    from collections import defaultdict

    def check_closed_path(arcs):
        # Step 1: Represent arcs and their connections
        graph = defaultdict(list)
        for arc in arcs:
            start, end = arc
            graph[start].append(end)
            graph[end].append(start)
        
        print("El grafo es: ",graph)
        
        # Step 2: Create a graph
        
        def dfs(node, visited, is_closed_path, path):
            if is_closed_path:
                #print(2)
                return True
            else:
                #print(3)
                visited.add(node)
                path.append(node)
                print("Añado al ciclo ",path[-1], path)
                for neighbor in graph[node]:
                    if neighbor not in visited and neighbor != path[-1]:
                        #print("3a")
                        is_closed_path = dfs(neighbor, visited, is_closed_path, path)
                    if neighbor == path[0] and len(path)>2:
                        #print("3b")
                        print("en el nodo ",node, "cierro el ciclo con ", neighbor)
                        path.append(neighbor)
                        return True
                    if is_closed_path:
                        return True
                    
                print("Elimino del ciclo ", path[-1], path[0:-1])
                path.pop()
                    
            return False
         
        visited = set()
        is_closed_path = False
        for node in graph:
            #print(1)
            if node not in visited:
                #print("1a")
                path = []
                is_closed_path = dfs(node, visited, is_closed_path, path)
                #print("salgo del dfs")
                if is_closed_path:  # Closed path condition
                    print("1aa")
                    return True, path
                
        return False, []
                
    distancesCities = distances_between_cities() # Calcula distancia entre ciudades
    orderedArcs=order_arcs(distancesCities) # Ordena los arcos de menor a mayor
    selectedArcs=[]
    path=[]
    for i,arc in enumerate(orderedArcs):
        print("\n     Iteración ", i)
        print("Añado el arco ", arc)
        selectedArcs.append(arc)
        print("Los arcos seleccionados son: ", selectedArcs)
        if has_overconnected_city(selectedArcs)==False:
            print("hay ciudades sobreconectadas?" ,has_overconnected_city(selectedArcs))
            is_closed_path,path = check_closed_path(selectedArcs)
            print("los arcos forman algún ciclo? ", is_closed_path)
            if is_closed_path and len(path)<len(coord_cities)+1:
                print("El ciclo ", path, "tiene dimensión menor a N")
                print("Eliminando último arco añadido...")
                selectedArcs.pop()
            else:
                print("Los arcos seleccionados no forman un ciclo")

        if len(path)==len(coord_cities)+1:
            print("El ciclo ", path, "tiene dimension N")
            break
    
    return selectedArcs, path

# main
if __name__=='__main__':
    solution, path = saving_cost_tsp(coord_cities, distance_earth)
    print("number of arcs in solution ", len(solution))
    print("number of nodes in path ", len(path))
    print("selected arcs: ",solution)
    print("path selected: ", path)

    plot_problem(dic_cities=coord_cities)
    plot_solution(solution=path,dic_cities=coord_cities)
    #plt.show()

    anim = animate_solution_arcs(cities=path, arcs=solution, cities_dict=coord_cities)
    matplotlib.rcParams['animation.ffmpeg_path'] = r'C:\PATH_Programs\ffmpeg.exe'
    FFwriter = matplotlib.animation.FFMpegWriter(fps=4)
    writergif = matplotlib.animation.PillowWriter(fps=4) 

    anim.save(r'C:\Users\Usuario\Documents\MII_MOIGE\MOD Y OPT PROBLEMAS GESTION\E1_Introduction_notebooks\animation.gif', writer=writergif)