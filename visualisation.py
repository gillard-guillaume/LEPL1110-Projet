import matplotlib.pyplot as plt

def read_joukowsky_points(filename):
    airfoil_points = []
    circles = [[] for _ in range(3)]
    
    with open(filename, 'r') as file:
        lines = file.readlines()
        circle_index = -1
        for line in lines:
            line = line.strip()
            if line.startswith("# Circles"):
                circle_index += 1
                continue
            if line and not line.startswith("#"):
                x, y = map(float, line.split())
                if circle_index == -1:
                    airfoil_points.append((x, y))
                else:
                    circles[circle_index].append((x, y))
    
    return airfoil_points, circles

def plot_joukowsky_shapes(filename):
    airfoil_points, circles = read_joukowsky_points(filename)
    
    # Extract X and Y coordinates
    x_airfoil, y_airfoil = zip(*airfoil_points)
    circle_coords = [list(zip(*circle)) for circle in circles]
    
    # Plot airfoil
    plt.figure(figsize=(8, 6))
    plt.plot(x_airfoil, y_airfoil, label='Joukowsky Airfoil', color='b')
    plt.scatter(x_airfoil, y_airfoil, color='b', s=10)  # Add points for airfoil
    
 
    
    plt.axhline(0, color='gray', linewidth=0.5)
    plt.axvline(0, color='gray', linewidth=0.5)
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.legend()
    plt.title("Joukowsky Transformation and Inner Circles")
    plt.grid(True)
    plt.axis("equal")
    plt.show()

# Run the visualization
plot_joukowsky_shapes("data/joukowsky.dat")
