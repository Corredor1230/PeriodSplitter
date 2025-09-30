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
    First column is treated as row headers (x-axis labels).
    """
    try:
        # Read CSV and use the first column as index
        df = pd.read_csv(filename, index_col=0)
    except FileNotFoundError:
        print(f"Error: The file '{filename}' was not found.")
        return

    if df.empty:
        print("The CSV file is empty. Nothing to plot.")
        return

    plt.figure(figsize=(14, 8))

    # Plot each column
    for column in df.columns:
        plt.plot(df.index, df[column], label=column)

    # Customize plot
    plt.title("Graph of Harmonic Data", fontsize=16)
    plt.xlabel("Row Header", fontsize=12)   # now x-axis shows your row headers
    plt.ylabel("Value", fontsize=12)
    plt.legend(title="Columns", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.tight_layout(rect=[0, 0, 0.88, 1])

    print("Displaying the plot...")
    plt.show()

# --- Main execution block ---
if __name__ == "__main__":
    # Define the name of the CSV file
    csv_file = "AMPEGuit6_110_F.csv"
    
    # 1. Create a sample CSV file to ensure the script is runnable
    # create_sample_csv(csv_file)
    
    # 2. Plot the data from the CSV file
    plot_data_from_csv(csv_file)
