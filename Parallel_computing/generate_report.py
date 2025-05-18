#!/usr/bin/env python3
import re
import pandas as pd
import matplotlib.pyplot as plt

# Function to extract runtime from files
def extract_time(filename):
    with open(filename, 'r') as f:
        content = f.read()
        # Look for real time in the format "real    0m12.345s"
        real_time = re.search(r'real\s+(\d+)m(\d+\.\d+)s', content)
        if real_time:
            minutes = int(real_time.group(1))
            seconds = float(real_time.group(2))
            return minutes * 60 + seconds
    return None

# Extract runtimes
files = ['serial.time', 'p10.time', 'p25.time', 'p50.time', 'p100.time']
methods = ['Serial', '10 Processes', '25 Processes', '50 Processes', '100 Processes']
times = [extract_time(f) for f in files]

# Calculate speedup (serial_time / parallel_time)
serial_time = times[0]
speedups = [serial_time / t if t else 0 for t in times]

# Create dataframe
data = {
    'Method': methods,
    'Runtime (seconds)': times,
    'Speedup': speedups
}
df = pd.DataFrame(data)

# Print report
print("\n--- PERFORMANCE REPORT ---")
print(df.to_string(index=False))
print("\n")

# Generate plots
plt.figure(figsize=(10, 5))

plt.subplot(1, 2, 1)
plt.bar(methods, times)
plt.title('Runtime Comparison')
plt.ylabel('Time (seconds)')
plt.xticks(rotation=45)

plt.subplot(1, 2, 2)
plt.bar(methods, speedups)
plt.title('Speedup Comparison')
plt.ylabel('Speedup (x times faster)')
plt.xticks(rotation=45)

plt.tight_layout()
plt.savefig('performance_report.png')
print("Generated performance_report.png")