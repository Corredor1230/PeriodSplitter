import struct
import tkinter as tk
from tkinter import filedialog, ttk
import os

def read_val(f, fmt):
    """Helper to read a single unpacked value from the binary file."""
    size = struct.calcsize(fmt)
    data = f.read(size)
    if not data or len(data) < size:
        return None
    return struct.unpack(fmt, data)[0]

def skip_bytes(f, count):
    f.seek(count, os.SEEK_CUR)

def load_sihat_to_tree(filepath, tree):
    """Parses the .sihat file and populates the tkinter Treeview."""
    # Clear existing tree
    for item in tree.get_children():
        tree.delete(item)

    root_node = tree.insert("", "end", text=f"📄 {os.path.basename(filepath)}", open=True)

    try:
        with open(filepath, 'rb') as f:
            # --- HEADER ---
            head_node = tree.insert(root_node, "end", text="Header", open=True)
            sr = read_val(f, 'I')
            fund_freq = read_val(f, 'f')
            tree.insert(head_node, "end", text=f"Sample Rate (uint32): {sr} Hz")
            tree.insert(head_node, "end", text=f"Fund. Freq (float): {fund_freq:.2f} Hz")

            # --- BLOCK 1: INDICES ---
            b1_node = tree.insert(root_node, "end", text="Block 1: Indices")
            n_indices = read_val(f, 'I')
            tree.insert(b1_node, "end", text=f"Count (uint32): {n_indices}")
            if n_indices > 0:
                skip_bytes(f, n_indices * 4)
                tree.insert(b1_node, "end", text=f"Data: [Array of {n_indices} uint32_t]")

            # --- BLOCK 2: AMPS ---
            b2_node = tree.insert(root_node, "end", text="Block 2: Amplitudes")
            n_harmonics = read_val(f, 'I')
            tree.insert(b2_node, "end", text=f"Num Harmonics (uint32): {n_harmonics}")
            for i in range(n_harmonics):
                h_len = read_val(f, 'I')
                if h_len > 0: skip_bytes(f, h_len * 4)
                tree.insert(b2_node, "end", text=f"Harmonic {i}: [Array of {h_len} floats]")

            # --- BLOCK 2b: PHASES ---
            b2b_node = tree.insert(root_node, "end", text="Block 2b: Phases")
            n_phases = read_val(f, 'I')
            tree.insert(b2b_node, "end", text=f"Num Phases (uint32): {n_phases}")
            for i in range(n_phases):
                p_len = read_val(f, 'I')
                if p_len > 0: skip_bytes(f, p_len * 4)
                tree.insert(b2b_node, "end", text=f"Phase List {i}: [Array of {p_len} floats]")

            # --- BLOCK 3: FREQUENCIES ---
            b3_node = tree.insert(root_node, "end", text="Block 3: Frequencies (ChangePoints)")
            n_freq_harms = read_val(f, 'I')
            tree.insert(b3_node, "end", text=f"Num Harmonics (uint32): {n_freq_harms}")
            for i in range(n_freq_harms):
                cp_len = read_val(f, 'I')
                # ChangePoint: index(uint32), value(float), ratio(float) = 12 bytes
                if cp_len > 0: skip_bytes(f, cp_len * 12) 
                tree.insert(b3_node, "end", text=f"Harmonic {i}: [Array of {cp_len} ChangePoints]")

            # --- BLOCK 4: TRANSIENTS ---
            b4_node = tree.insert(root_node, "end", text="Block 4: Transients")
            t_start, t_end = read_val(f, 'I'), read_val(f, 'I')
            tree.insert(b4_node, "end", text=f"Range (uint32): {t_start} -> {t_end}")
            
            env_hop = read_val(f, 'I')
            spec_hop = read_val(f, 'I')
            spec_win = read_val(f, 'I')
            spec_bins = read_val(f, 'I')
            params_node = tree.insert(b4_node, "end", text="FFT Parameters")
            tree.insert(params_node, "end", text=f"envHopSize: {env_hop}, specHopSize: {spec_hop}")
            tree.insert(params_node, "end", text=f"specWindowSize: {spec_win}, specNumBins: {spec_bins}")

            rise_time = read_val(f, 'i')
            peak_amp = read_val(f, 'f')
            tree.insert(b4_node, "end", text=f"Rise Time (int32): {rise_time}")
            tree.insert(b4_node, "end", text=f"Peak Amp (float): {peak_amp:.4f}")

            # Optimized Envelope
            env_size = read_val(f, 'I')
            env_node = tree.insert(b4_node, "end", text=f"Optimized Amp Envelope (Size: {env_size})")
            if env_size > 0:
                first_idx = read_val(f, 'i')
                skip_bytes(f, env_size * 4)
                tree.insert(env_node, "end", text=f"First Index (int32): {first_idx}")
                tree.insert(env_node, "end", text=f"Data: [Array of {env_size} floats]")

            # Float Vectors
            for name in ["Centroid", "Flatness", "MainOvertones"]:
                v_size = read_val(f, 'I')
                if v_size > 0: skip_bytes(f, v_size * 4)
                tree.insert(b4_node, "end", text=f"{name}: [Array of {v_size} floats]")

            # Partials
            num_partials = read_val(f, 'I')
            band_env_size = read_val(f, 'I')
            if band_env_size > 0: skip_bytes(f, band_env_size * 4)
            tree.insert(b4_node, "end", text=f"Num Partials (uint32): {num_partials}")
            tree.insert(b4_node, "end", text=f"Band Envelopes: [Array of {band_env_size} floats]")

            # Trajectories (UPDATED)
            num_traj = read_val(f, 'I')
            traj_node = tree.insert(b4_node, "end", text=f"Trajectories (Count: {num_traj})")
            for i in range(num_traj):
                # Read the target peak variables explicitly (double, double, double)
                t_freq = read_val(f, 'd')
                t_mag = read_val(f, 'd')
                t_amp = read_val(f, 'd')
                
                # Read the trajectory floor (float)
                t_floor = read_val(f, 'f')
                
                env_len = read_val(f, 'I')
                
                # TrackedPoint: sampleIndex(int32), freq(float), crestFactor(float) = 12 bytes
                if env_len > 0:
                    skip_bytes(f, env_len * 12) 
                
                t_sub = tree.insert(traj_node, "end", text=f"Trajectory {i}")
                tree.insert(t_sub, "end", text=f"Target Peak: Freq={t_freq:.1f}Hz, Mag={t_mag:.4f}, Amp={t_amp:.4f}")
                tree.insert(t_sub, "end", text=f"Floor: {t_floor:.4f}")
                tree.insert(t_sub, "end", text=f"Envelope: [Array of {env_len} TrackedPoints (int32, float, float)]")

            # --- BLOCK 5: RMS ---
            b5_node = tree.insert(root_node, "end", text="Block 5: RMS", open=True)
            h_rms = read_val(f, 'f')
            t_rms = read_val(f, 'f')
            tree.insert(b5_node, "end", text=f"Harmonic RMS (float): {h_rms:.4f}")
            tree.insert(b5_node, "end", text=f"Transient RMS (float): {t_rms:.4f}")

    except Exception as e:
        tree.insert(root_node, "end", text=f"❌ Error parsing file: {e}")

def open_file_dialog(tree):
    filepath = filedialog.askopenfilename(
        title="Select a .sihat File",
        filetypes=[("Sihat Files", "*.sihat"), ("All Files", "*.*")]
    )
    if filepath:
        load_sihat_to_tree(filepath, tree)

def main():
    # Setup standard Window
    root = tk.Tk()
    root.title("Sihat Binary Viewer")
    root.geometry("650x750")

    # Layout styling
    style = ttk.Style()
    style.configure("Treeview", font=('Consolas', 10), rowheight=25)
    style.configure("Treeview.Heading", font=('Consolas', 11, 'bold'))

    # Top Frame for Button
    top_frame = tk.Frame(root, pady=10)
    top_frame.pack(fill="x")
    
    # Treeview setup
    tree_frame = tk.Frame(root)
    tree_frame.pack(fill="both", expand=True, padx=10, pady=10)
    
    tree = ttk.Treeview(tree_frame, show="tree")
    tree.pack(side="left", fill="both", expand=True)
    
    scrollbar = ttk.Scrollbar(tree_frame, orient="vertical", command=tree.yview)
    scrollbar.pack(side="right", fill="y")
    tree.configure(yscrollcommand=scrollbar.set)

    # Open button
    btn = tk.Button(top_frame, text="📁 Open .sihat File", font=('Arial', 12), 
                    command=lambda: open_file_dialog(tree), bg="#e0e0e0", padx=10, pady=5)
    btn.pack()

    tree.insert("", "end", text="Waiting for file... Click 'Open .sihat File' above.")
    root.mainloop()

if __name__ == "__main__":
    main()