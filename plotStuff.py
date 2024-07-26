import matplotlib.pyplot as plt

def plot_curves(x, y, z, labely, labelz):
    """
    Plot two curves using a list of values (x, y, z).

    Parameters:
    - x: list of values for the x-axis.
    - y: list of values for the y-axis of the first curve.
    - z: list of values for the y-axis of the second curve.
    """
    plt.figure(figsize=(10, 6))

    # Plot the first curve (x, y)
    plt.plot(x, y, label=labely, marker='o', linestyle='-', color='blue')

    # Plot the second curve (x, z)
    plt.plot(x, z, label=labelz, marker='s', linestyle='--', color='red')

    # Add title and labels
    plt.title('Plot of Two Curves')
    plt.xlabel('run')
    plt.ylabel('value')

    # Add legend
    plt.legend()

    # Show grid
    plt.grid(True)

    # Show plot
    plt.show()

