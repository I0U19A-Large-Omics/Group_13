try:
    import pandas as pd
    import plotly.express as px
except ImportError:
    print("Required packages not found. Please install them with:")
    print("    pip3 install --user pandas plotly")
    exit(1)

# Load CSV
df = pd.read_csv('test_time/benchmark_results.csv')

# Plot Runtime vs Number of Processes
fig_runtime = px.line(df, x='Processes', y='Runtime_seconds', markers=True,
                      title='Runtime vs Number of Processes',
                      labels={'Processes': 'Number of Parallel Processes', 'Runtime_seconds': 'Runtime (seconds)'})
fig_runtime.update_layout(yaxis=dict(title='Runtime (seconds)'), xaxis=dict(dtick=1))

# Plot Speedup vs Number of Processes
fig_speedup = px.line(df, x='Processes', y='Speedup', markers=True, 
                      title='Speedup vs Number of Processes',
                      labels={'Processes': 'Number of Parallel Processes', 'Speedup': 'Speedup (relative to serial)'})
fig_speedup.update_layout(yaxis=dict(title='Speedup (relative to serial)'), xaxis=dict(dtick=1))

# Save plots as interactive HTML files
print("Saving runtime_plot.html...")
fig_runtime.write_html("plots/runtime_plot.html")
print("Saved runtime_plot.html!")

print("Saving speedup_plot.html...")
fig_speedup.write_html("plots/speedup_plot.html")
print("Saved speedup_plot.html!")

# Optional: show plots in interactive window (uncomment if you want this)
#fig_runtime.show()
#fig_speedup.show()
