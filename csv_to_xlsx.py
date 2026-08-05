import os
import pandas as pd


def csv_to_xlsx(directory_path):
    # Loop through each file in the directory
    for filename in os.listdir(directory_path):
        # Check if the file is a .txt file
        if filename.endswith('.csv'):
            # Define the file paths
            csv_file_path = os.path.join(directory_path, filename)
            xlsx_file_path = os.path.join(directory_path, f"{os.path.splitext(filename)[0]}.xlsx")

            # Read the text file into a DataFrame
            try:
                df = pd.read_csv(csv_file_path)
                # Save DataFrame to an Excel file
                df.to_excel(xlsx_file_path, index=False)
                print(f"Converted '{filename}' to '{xlsx_file_path}'")
            except Exception as e:
                print(f"Failed to convert '{filename}': {e}")


# Example usage with a Windows-style path
directory_path = "C:\\Users\\boxx_\\Desktop\\gnomAD\\New folder"  # Replace with your actual directory path
csv_to_xlsx(directory_path)
