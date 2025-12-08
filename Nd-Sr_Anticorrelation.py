import numpy as np
import matplotlib.pyplot as plt

# Time array from 4.5 Ga to present (0 Ga)
t = np.linspace(4.5, 0, 500)  
extract_time = 2.9

# Convert to years before present (4.5 Ga = 4.5e9 years ago)
years_bp = t * 1e9

# Decay constants
lambda_hf = 1.9e-11  # for 176Lu -> 176Hf
lambda_nd = 6.54e-12  # for 147Sm -> 143Nd

# Initial ratios (approximate values)
R0_hf = 0.2798  # Initial 176Hf/177Hf
R0_nd = 0.5068  # Initial 143Nd/144Nd

# Parent/daughter ratios (schematic values)
PD_bse_hf = 0.033  # Lu/Hf for BSE
PD_bse_nd = 0.196  # Sm/Nd for BSE

PD_depleted_hf = 0.050  # Higher Lu/Hf for depleted mantle
PD_depleted_nd = 0.250  # Higher Sm/Nd for depleted mantle  
PD_crust_hf = 0.010    # Lower Lu/Hf for continental crust
PD_crust_nd = 0.120    # Lower Sm/Nd for continental crust

# Initialize arrays
R_bse_hf = np.zeros_like(t)
R_bse_nd = np.zeros_like(t)
R_depleted_hf = np.zeros_like(t)
R_depleted_nd = np.zeros_like(t)
R_crust_hf = np.zeros_like(t)
R_crust_nd = np.zeros_like(t)

# Find extraction index
idx_extract = np.argmin(np.abs(t - extract_time))

# Calculate evolution
for i, time_bp in enumerate(years_bp):
    # Time elapsed from 4.5 Ga to current time
    t_elapsed = (4.5e9 - time_bp)
    
    # BSE evolution (always follows average path)
    R_bse_hf[i] = R0_hf + PD_bse_hf * (np.exp(lambda_hf * t_elapsed) - 1)
    R_bse_nd[i] = R0_nd + PD_bse_nd * (np.exp(lambda_nd * t_elapsed) - 1)
    
    if t[i] >= extract_time:
        # Before extraction - all reservoirs are the same as BSE
        R_depleted_hf[i] = R_bse_hf[i]
        R_depleted_nd[i] = R_bse_nd[i]
        R_crust_hf[i] = R_bse_hf[i]
        R_crust_nd[i] = R_bse_nd[i]
    else:
        # After extraction - separate evolution
        # Time since extraction
        t_since_extract = (extract_time - t[i]) * 1e9
        
        # Start from the value at extraction time
        R_extract_hf = R_bse_hf[idx_extract]
        R_extract_nd = R_bse_nd[idx_extract]
        
        # Evolve separately with different parent/daughter ratios
        R_depleted_hf[i] = R_extract_hf + PD_depleted_hf * (np.exp(lambda_hf * t_since_extract) - 1)
        R_depleted_nd[i] = R_extract_nd + PD_depleted_nd * (np.exp(lambda_nd * t_since_extract) - 1)
        R_crust_hf[i] = R_extract_hf + PD_crust_hf * (np.exp(lambda_hf * t_since_extract) - 1)
        R_crust_nd[i] = R_extract_nd + PD_crust_nd * (np.exp(lambda_nd * t_since_extract) - 1)

# Plot Hf evolution
fig1, ax1 = plt.subplots(figsize=(8, 5))
ax1.plot(t, R_bse_hf, 'k-', label='Bulk Silicate Earth (BSE)', linewidth=2)
ax1.plot(t, R_depleted_hf, 'r--', label='Depleted Mantle', linewidth=2)
ax1.plot(t, R_crust_hf, 'b:', label='Continental Crust', linewidth=2)
ax1.axvline(extract_time, color='gray', linestyle='-', alpha=0.7, linewidth=1)
ax1.text(extract_time + 0.05, R0_hf + 0.001, 'Crust Extraction\n2.9 Ga', va='bottom')
ax1.set_xlabel('Time (Ga before present)')
ax1.set_ylabel('$^{176}$Hf/$^{177}$Hf')
ax1.set_title('Evolution of $^{176}$Hf/$^{177}$Hf Isotopic Ratio')
ax1.set_xlim(4.5, 0)
ax1.invert_xaxis()
ax1.legend()
ax1.grid(True, alpha=0.3)
fig1.tight_layout()

# Plot Nd evolution
fig2, ax2 = plt.subplots(figsize=(8, 5))
ax2.plot(t, R_bse_nd, 'k-', label='Bulk Silicate Earth (BSE)', linewidth=2)
ax2.plot(t, R_depleted_nd, 'r--', label='Depleted Mantle', linewidth=2)
ax2.plot(t, R_crust_nd, 'b:', label='Continental Crust', linewidth=2)
ax2.axvline(extract_time, color='gray', linestyle='-', alpha=0.7, linewidth=1)
ax2.text(extract_time + 0.05, R0_nd + 0.003, 'Crust Extraction\n2.9 Ga', va='bottom')
ax2.set_xlabel('Time (Ga before present)')
ax2.set_ylabel('$^{143}$Nd/$^{144}$Nd')
ax2.set_title('Evolution of $^{143}$Nd/$^{144}$Nd Isotopic Ratio')
ax2.set_xlim(4.5, 0)
ax2.invert_xaxis()
ax2.legend()
ax2.grid(True, alpha=0.3)
fig2.tight_layout()

plt.show()

# Save figures
fig1.savefig('evolution_Hf_corrected.png', dpi=150, bbox_inches='tight')
fig2.savefig('evolution_Nd_corrected.png', dpi=150, bbox_inches='tight')