import matplotlib.pyplot as plt

# HP Summary Results (Cutoff, Best Model, AUC, StdErr)
hp_summary = [
    (2.5, "power-split-1-4", 0.9600, 0.0030),
    (3.0, "power-split-0-2", 0.9600, 0.0020),
    (3.5, "power-split-0-2", 0.9500, 0.0020),
    (4.0, "power-split-2-5", 0.9400, 0.0030)
]

# Extract data for plotting
cutoffs = [x[0] for x in hp_summary]
auc_values = [x[2] for x in hp_summary]
stderr_values = [x[3] for x in hp_summary]

# Create the plot
plt.figure(figsize=(6, 4))

# Scatter plot with error bars
plt.errorbar(
    cutoffs,
    auc_values,
    yerr=stderr_values,
    fmt='o-',  # Circle markers with lines connecting
    color='blue',
    ecolor='black',
    elinewidth=1.2,
    capsize=5,
    label='HP model performance vs. contact cutoff threshold'
)

# Customize labels and title
plt.title("HP model performance vs. cutoff distance", fontsize=14)
plt.xlabel("Contact map cutoff distance", fontsize=12)
plt.ylabel("Average test AUC across 185 randomized test-train splits", fontsize=12)

# Grid and ticks formatting
plt.grid(True, linestyle='--', alpha=0.6)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)

# Add legend
plt.legend(fontsize=12)

# Ensure a compact layout
plt.tight_layout()

# Save the plot
plt.savefig("hp_model_performance_vs_cutoff_distance_plot.png", dpi=300)

# Display the plot
plt.show()
