import os
import time
os.environ["TK_SILENCE_DEPRECATION"] = "1"

import sys
import numpy as np
import tkinter as tk
from tkinter import ttk, messagebox

import matplotlib
if sys.platform == "darwin":
    matplotlib.use("MacOSX")
else:
    matplotlib.use("TkAgg")

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


def compute_evolution(extract_time: float = 2.9):
    t = np.linspace(4.5, 0, 500)
    years_bp = t * 1e9

    lambda_hf = 1.9e-11
    lambda_nd = 6.54e-12

    R0_hf = 0.2798
    R0_nd = 0.5068

    PD_bse_hf = 0.033
    PD_bse_nd = 0.196
    PD_depleted_hf = 0.050
    PD_depleted_nd = 0.250
    PD_crust_hf = 0.010
    PD_crust_nd = 0.120

    R_bse_hf = np.zeros_like(t)
    R_bse_nd = np.zeros_like(t)
    R_depleted_hf = np.zeros_like(t)
    R_depleted_nd = np.zeros_like(t)
    R_crust_hf = np.zeros_like(t)
    R_crust_nd = np.zeros_like(t)

    idx_extract = np.argmin(np.abs(t - extract_time))

    for i, time_bp in enumerate(years_bp):
        t_elapsed = 4.5e9 - time_bp

        R_bse_hf[i] = R0_hf + PD_bse_hf * (np.exp(lambda_hf * t_elapsed) - 1)
        R_bse_nd[i] = R0_nd + PD_bse_nd * (np.exp(lambda_nd * t_elapsed) - 1)

        if t[i] >= extract_time:
            R_depleted_hf[i] = R_bse_hf[i]
            R_depleted_nd[i] = R_bse_nd[i]
            R_crust_hf[i] = R_bse_hf[i]
            R_crust_nd[i] = R_bse_nd[i]
        else:
            t_since_extract = (extract_time - t[i]) * 1e9
            R_extract_hf = R_bse_hf[idx_extract]
            R_extract_nd = R_bse_nd[idx_extract]

            R_depleted_hf[i] = R_extract_hf + PD_depleted_hf * (np.exp(lambda_hf * t_since_extract) - 1)
            R_depleted_nd[i] = R_extract_nd + PD_depleted_nd * (np.exp(lambda_nd * t_since_extract) - 1)

            R_crust_hf[i] = R_extract_hf + PD_crust_hf * (np.exp(lambda_hf * t_since_extract) - 1)
            R_crust_nd[i] = R_extract_nd + PD_crust_nd * (np.exp(lambda_nd * t_since_extract) - 1)

    return {
        "t": t,
        "R_bse_hf": R_bse_hf,
        "R_depleted_hf": R_depleted_hf,
        "R_crust_hf": R_crust_hf,
        "R_bse_nd": R_bse_nd,
        "R_depleted_nd": R_depleted_nd,
        "R_crust_nd": R_crust_nd,
        "extract_time": extract_time,
        "R0_hf": R0_hf,
        "R0_nd": R0_nd,
    }


def build_figures(data):
    t = data["t"]
    extract_time = data["extract_time"]

    fig1, ax1 = plt.subplots(figsize=(8, 5))
    ax1.plot(t, data["R_bse_hf"], "k-", label="Bulk Silicate Earth (BSE)", linewidth=2)
    ax1.plot(t, data["R_depleted_hf"], "r--", label="Depleted Mantle", linewidth=2)
    ax1.plot(t, data["R_crust_hf"], "b:", label="Continental Crust", linewidth=2)
    ax1.axvline(extract_time, color="gray", linestyle="-", alpha=0.7, linewidth=1)
    ax1.text(
        extract_time + 0.05,
        np.interp(extract_time, t[::-1], data["R_bse_hf"][::-1]),
        f"Crust Extraction\n{extract_time:.2f} Ga",
        va="bottom",
    )
    ax1.set_xlabel("Time (Ga before present)")
    ax1.set_ylabel("176Hf/177Hf")
    ax1.set_title("Evolution of 176Hf/177Hf Isotopic Ratio")
    ax1.legend(loc="upper right")
    ax1.set_xlim(4.5, 0)
    ax1.invert_xaxis()
    ax1.grid(True, alpha=0.3)
    fig1.tight_layout()

    fig2, ax2 = plt.subplots(figsize=(8, 5))
    ax2.plot(t, data["R_bse_nd"], "k-", label="Bulk Silicate Earth (BSE)", linewidth=2)
    ax2.plot(t, data["R_depleted_nd"], "r--", label="Depleted Mantle", linewidth=2)
    ax2.plot(t, data["R_crust_nd"], "b:", label="Continental Crust", linewidth=2)
    ax2.axvline(extract_time, color="gray", linestyle="-", alpha=0.7, linewidth=1)
    ax2.text(
        extract_time + 0.05,
        np.interp(extract_time, t[::-1], data["R_bse_nd"][::-1]),
        f"Crust Extraction\n{extract_time:.2f} Ga",
        va="bottom",
    )
    ax2.set_xlabel("Time (Ga before present)")
    ax2.set_ylabel("143Nd/144Nd")
    ax2.set_title("Evolution of 143Nd/144Nd Isotopic Ratio")
    ax2.legend(loc="upper right")
    ax2.set_xlim(4.5, 0)
    ax2.invert_xaxis()
    ax2.grid(True, alpha=0.3)
    fig2.tight_layout()

    return fig1, fig2


class IsoEvolutionApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Hf/Nd Isotopic Evolution")

        control_frame = ttk.Frame(root, padding=10)
        control_frame.pack(side=tk.TOP, fill=tk.X)

        ttk.Label(control_frame, text="Crust extraction time (Ga):").pack(side=tk.LEFT)
        self.extract_var = tk.DoubleVar(value=2.9)
        ttk.Entry(control_frame, textvariable=self.extract_var, width=6).pack(side=tk.LEFT, padx=6)

        ttk.Button(control_frame, text="Generate plots", command=self.generate).pack(side=tk.LEFT, padx=5)
        ttk.Button(control_frame, text="Save PNGs", command=self.save_pngs).pack(side=tk.LEFT, padx=5)

        ttk.Label(control_frame, text="Show:").pack(side=tk.LEFT, padx=(10, 2))
        self.plot_selector = tk.StringVar(value="hf")
        ttk.Radiobutton(control_frame, text="Hf", variable=self.plot_selector, value="hf", command=self.show_selected).pack(side=tk.LEFT)
        ttk.Radiobutton(control_frame, text="Nd", variable=self.plot_selector, value="nd", command=self.show_selected).pack(side=tk.LEFT)

        self.plot_frame = ttk.Frame(root, padding=10)
        self.plot_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.canvas1 = None 
        self.canvas2 = None 
        self.fig1 = None
        self.fig2 = None

        self.generate()

    def clear_canvases(self):
        if self.canvas1:
            self.canvas1.get_tk_widget().destroy()
        if self.canvas2:
            self.canvas2.get_tk_widget().destroy()
        self.canvas1 = None
        self.canvas2 = None

    def generate(self):
        try:
            extract_time = float(self.extract_var.get())
        except ValueError:
            messagebox.showerror("Invalid input", "Extraction time must be a number.")
            return

        data = compute_evolution(extract_time)
        fig1, fig2 = build_figures(data)

        self.clear_canvases()
        self.fig1, self.fig2 = fig1, fig2

        self.canvas1 = FigureCanvasTkAgg(self.fig1, master=self.plot_frame)
        self.canvas1.draw()
        self.canvas2 = FigureCanvasTkAgg(self.fig2, master=self.plot_frame)
        self.canvas2.draw()

        self.show_selected()

    def save_pngs(self):
        reports_dir = os.path.join(os.getcwd(), "photo-reports")
        os.makedirs(reports_dir, exist_ok=True)
        timestamp = time.strftime("%Y%m%d-%H%M%S")
        hf_path = os.path.join(reports_dir, f"evolution_Hf({timestamp}).png")
        nd_path = os.path.join(reports_dir, f"evolution_Nd({timestamp}).png")
        self.fig1.savefig(hf_path, dpi=150, bbox_inches="tight")
        self.fig2.savefig(nd_path, dpi=150, bbox_inches="tight")
        messagebox.showinfo("Saved", f"Saved:\n{hf_path}\n{nd_path}")

    def show_selected(self):
        """Show only the selected plot (Hf or Nd) to avoid a huge window."""
        if self.canvas1:
            self.canvas1.get_tk_widget().pack_forget()
        if self.canvas2:
            self.canvas2.get_tk_widget().pack_forget()

        if self.plot_selector.get() == "nd":
            self.canvas2.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        else:
            self.canvas1.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)


def main():
    root = tk.Tk()
    app = IsoEvolutionApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
