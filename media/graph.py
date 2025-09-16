import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

def create_sample_csv(filename="data.csv", num_rows=820, num_cols=16):
    """
    Generates a sample CSV file with random data if it doesn't already exist.
    This allows the script to be run without needing a pre-existing data file.
    """
    if os.path.exists(filename):
        print(f"'{filename}' already exists. Skipping creation.")
        return

    print(f"Generating a sample CSV file: '{filename}'...")
    # Create column headers like "Column 1", "Column 2", etc.
    headers = [f"Column_{i+1}" for i in range(num_cols)]
    
    # Generate random data
    # We use a cumulative sum to make the lines look more like time-series data
    data = np.random.randn(num_rows, num_cols).cumsum(axis=0)
    
    # Create a pandas DataFrame
    df = pd.DataFrame(data, columns=headers)
    
    # Save the DataFrame to a CSV file
    df.to_csv(filename, index=False)
    print("Sample CSV file created successfully.")

def plot_data_from_csv(filename="data.csv"):
    """
    Reads data from a CSV file and plots each column as a separate line graph.
    """
    try:
        # Read the data from the CSV file into a pandas DataFrame
        df = pd.read_csv(filename)
    except FileNotFoundError:
        print(f"Error: The file '{filename}' was not found.")
        print("Please make sure the CSV file is in the same directory as the script,")
        print("or run the script once to generate a sample 'data.csv' file.")
        return

    # Check if there are columns to plot
    if df.empty:
        print("The CSV file is empty. Nothing to plot.")
        return

    # Create a figure and axes for the plot
    # plt.figure() creates a container for the plot
    # figsize=(width, height) in inches
    plt.figure(figsize=(14, 8))

    # Loop through each column in the DataFrame and plot it
    for column in df.columns:
        plt.plot(df.index, df[column], label=column)

    # --- Customize the plot's appearance ---
    
    # Add a title to the graph
    plt.title("Graph of Harmonic Data", fontsize=16)
    
    # Add labels to the X and Y axes
    plt.xlabel("Row Number (Index)", fontsize=12)
    plt.ylabel("Value", fontsize=12)
    
    # Add a legend to identify which line corresponds to which column
    # We adjust the position to be outside the plot area to avoid covering lines
    plt.legend(title="Columns", bbox_to_anchor=(1.02, 1), loc='upper left')
    
    # Add a grid for better readability
    plt.grid(True, linestyle='--', alpha=0.6)
    
    # Adjust layout to prevent labels from being cut off
    plt.tight_layout(rect=[0, 0, 0.88, 1]) # Adjust right margin for legend
    
    # Display the plot
    print("Displaying the plot...")
    plt.show()

# --- Main execution block ---
if __name__ == "__main__":
    # Define the name of the CSV file
    csv_file = "AMPEGuit6_82_F.csv"
    
    # 1. Create a sample CSV file to ensure the script is runnable
    # create_sample_csv(csv_file)
    
    # 2. Plot the data from the CSV file
    plot_data_from_csv(csv_file)
