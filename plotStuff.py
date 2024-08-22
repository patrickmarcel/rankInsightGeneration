import matplotlib.pyplot as plt


def plot_curves(data, curve_names, xlabel, ylabel, title):
    """
    Plot an arbitrary number of curves from data and assign names to the curves.

    Parameters:
    - data: list of tuples, where the first value of each tuple is the x-axis values,
      the second value is the y-axis values for the first curve, the third value is the y-axis values
      for the second curve, and so on.
    - curve_names: list of strings, names for the curves to plot.
    """
    # Ensure the number of curve names matches the number of y-axis datasets
    num_curves = len(data[0]) - 1
    if len(curve_names) != num_curves:
        raise ValueError("The number of curve names must match the number of y-axis datasets in data.")

    # Extract x and y values for each curve
    x_values = [t[0] for t in data]
    y_values = [[] for _ in range(num_curves)]

    for t in data:
        for i in range(num_curves):
            y_values[i].append(t[i + 1])

    # Plot the curves
    for i in range(num_curves):
        plt.plot(x_values, y_values[i], label=curve_names[i])

    # Add labels and legend
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.title(title)
    plt.grid(True)
    plt.show()




def plot_2curves(x, y, z, labely, labelz, title):
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
    plt.title(title)
    plt.xlabel('run')
    plt.ylabel('value')

    # Add legend
    plt.legend()

    # Show grid
    plt.grid(True)

    # Show plot
    plt.show()

