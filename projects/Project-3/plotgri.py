from pathlib import Path
import numpy as np
import matplotlib.tri as mtri
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import os


mpl.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],

    "font.size": 30, 
    "axes.titlesize": 30,
    "axes.labelsize": 30,
    "xtick.labelsize": 24,
    "ytick.labelsize": 24,
    "legend.fontsize": 24,

    "axes.linewidth": 2.0,
    "lines.linewidth": 2.5,

    "figure.titlesize": 32
})

#-----------------------------------------------------------
def readgri(fname):
    with open(fname, 'r') as f:
        Nn, Ne, dim = [int(s) for s in f.readline().split()]
        if dim != 2:
            raise ValueError(f"Expected dim=2, got {dim}")

        V = np.array([[float(s) for s in f.readline().split()] for _ in range(Nn)])

        NB = int(f.readline())
        B_raw = []
        Bname = []
        for _ in range(NB):
            s = f.readline().split()
            Nb = int(s[0])
            Bname.append(s[2])
            Bi = np.array([[int(t) for t in f.readline().split()] for _ in range(Nb)], dtype=int)
            B_raw.append(Bi)

        Ne0 = 0
        E_blocks = []
        while Ne0 < Ne:
            s = f.readline().split()
            ne_blk = int(s[0])
            Ei = np.array([[int(t) for t in f.readline().split()] for _ in range(ne_blk)], dtype=int)
            E_blocks.append(Ei)
            Ne0 += ne_blk

    E = np.concatenate(E_blocks, axis=0) if E_blocks else np.empty((0, 3), dtype=int)

    # Detect 0-based vs 1-based indexing
    all_arrays = [E] + B_raw
    amin = min(arr.min() for arr in all_arrays if arr.size > 0)
    amax = max(arr.max() for arr in all_arrays if arr.size > 0)

    if amin == 0 and amax <= Nn - 1:
        shift = 0
    elif amin >= 1 and amax <= Nn:
        shift = 1
    else:
        raise ValueError(
            f"Could not determine indexing convention in {fname}. "
            f"Index range is [{amin}, {amax}] but node count is {Nn}."
        )

    E = E - shift
    B = [Bi - shift for Bi in B_raw]

    # Final sanity checks
    if E.min() < 0 or E.max() >= Nn:
        raise ValueError(
            f"Connectivity out of range after shifting in {fname}: "
            f"[{E.min()}, {E.max()}], expected within [0, {Nn-1}]"
        )

    return {'V': V, 'E': E, 'B': B, 'Bname': Bname}
def read_xy_txt(path: str):
    arr = np.loadtxt(path, dtype=float)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    return arr[:, :2]

#-----------------------------------------------------------
def read_scalar_field_maybe_header(path):
    """
    Accepts either:
      A) N lines of floats
      B) first line is integer N, followed by N floats
    Returns a 1D float array.
    """
    vals = np.loadtxt(path, dtype=float, ndmin=1)
    if vals.size == 0:
        return vals

    # If first entry is an integer-like count and matches remaining length, drop it.
    n0 = int(round(vals[0]))
    if abs(vals[0] - n0) < 1e-12 and (vals.size - 1) == n0:
        return vals[1:]
    return vals

#-----------------------------------------------------------
def plot_wall_distance(Mesh, dist, fname, show_mesh=True, use_log=False, clim=None, plot_sizing=False):
    V = Mesh['V']; E = Mesh['E']

    if dist.shape[0] != V.shape[0]:
        raise ValueError(f"dist length {dist.shape[0]} != number of nodes {V.shape[0]}")

    if plot_sizing:
        if use_log:
            field = np.log10(np.maximum(field, 1e-16))
            title = "log10(size function h)"
        else:
            field = dist.copy()
            title = "Size function h"
    else:
        if use_log:
            field = np.log10(np.maximum(field, 1e-16))
            title = "log10(wall distance)"
        else:
            field = dist.copy()
            title = "Wall distance"
    

    fig = plt.figure(figsize=(12, 8))
    tpc = plt.tripcolor(V[:, 0], V[:, 1], E, field, shading='gouraud')

    if clim is not None:
        tpc.set_clim(clim[0], clim[1])

    plt.colorbar(tpc, shrink=0.85, label=title)

    if show_mesh:
        plt.triplot(V[:, 0], V[:, 1], E, 'k-', linewidth=0.2)

    plt.axis('equal')
    plt.tight_layout()
    plt.savefig(fname, dpi=300)
    plt.close(fig)

def plot_mesh_with_blades(mesh, blade_upper, blade_lower, out_png):
    V = mesh["V"]; E = mesh["E"]

    fig = plt.figure(figsize=(10, 10))
    plt.triplot(V[:,0], V[:,1], E, linewidth=0.3, color='black', alpha=1)
    up = read_xy_txt(blade_upper)
    lo = read_xy_txt(blade_lower)
    plt.legend(loc="best")
    plt.axis("equal")
    plt.tight_layout()
    plt.savefig(out_png, dpi=400)
    plt.close(fig)

# -----------------------------------------------------------
def read_curved_edges_txt(path: str):
    """
    Reads files of the form:

    Curve5 44 45
    x0 y0
    x1 y1
    x2 y2

    Curve5 45 46
    x0 y0
    x1 y1
    x2 y2
    ...

    Returns a list of dicts:
      [{"curve": str, "n1": int, "n2": int, "pts": np.ndarray}, ...]
    """
    edges = []
    with open(path, "r") as f:
        lines = [ln.strip() for ln in f if ln.strip()]

    i = 0
    while i < len(lines):
        header = lines[i].split()
        if len(header) != 3:
            raise ValueError(f"Bad header line in {path}: {lines[i]}")
        curve_name = header[0]
        n1 = int(header[1])
        n2 = int(header[2])
        i += 1

        pts = []
        while i < len(lines):
            parts = lines[i].split()
            if len(parts) == 3:
                # next header
                break
            if len(parts) != 2:
                raise ValueError(f"Bad point line in {path}: {lines[i]}")
            pts.append([float(parts[0]), float(parts[1])])
            i += 1

        edges.append({
            "curve": curve_name,
            "n1": n1,
            "n2": n2,
            "pts": np.array(pts, dtype=float)
        })

    return edges

# -----------------------------------------------------------
# -----------------------------------------------------------
def lagrange_basis_1d(xi_eval, xi_nodes):
    """
    Evaluate 1D Lagrange basis polynomials at xi_eval.

    Parameters
    ----------
    xi_eval : (m,) array
        Evaluation points in reference 1D space.
    xi_nodes : (n,) array
        Interpolation nodes in reference 1D space.

    Returns
    -------
    L : (m, n) array
        L[k, i] = l_i(xi_eval[k])
    """
    xi_eval = np.asarray(xi_eval, dtype=float)
    xi_nodes = np.asarray(xi_nodes, dtype=float)

    m = xi_eval.size
    n = xi_nodes.size
    L = np.ones((m, n), dtype=float)

    for i in range(n):
        for j in range(n):
            if j == i:
                continue
            L[:, i] *= (xi_eval - xi_nodes[j]) / (xi_nodes[i] - xi_nodes[j])

    return L


# -----------------------------------------------------------
def evaluate_curved_edge(edge_pts, nplot=200):
    """
    Reconstruct a curved edge from its Lagrange nodes using 1D shape functions.

    Parameters
    ----------
    edge_pts : (q+1, 2) array
        Physical coordinates of edge interpolation nodes.
        Assumed ordered along the edge.
    nplot : int
        Number of dense plotting points.

    Returns
    -------
    xy : (nplot, 2) array
        Interpolated curved edge coordinates.
    """
    edge_pts = np.asarray(edge_pts, dtype=float)
    q = edge_pts.shape[0] - 1

    if q < 1:
        raise ValueError("Need at least 2 points on an edge.")

    # Equally spaced reference nodes on [-1, 1]
    xi_nodes = np.linspace(-1.0, 1.0, q + 1)
    xi_plot = np.linspace(-1.0, 1.0, nplot)

    L = lagrange_basis_1d(xi_plot, xi_nodes)   # (nplot, q+1)
    xy = L @ edge_pts                          # (nplot, 2)

    return xy



# -----------------------------------------------------------
def plot_curved_boundary_multiQ(
    mesh,
    blade_upper,
    blade_lower,
    curved_files_dict,   # {"Q1": (upper, lower), "Q2": (...), "Q3": (...)}
    out_png,
    nplot_per_edge=300,
    show_mesh=True,
    show_edge_nodes=False
):
    """
    Plot Q1, Q2, Q3 curved edges side-by-side (3 columns).
    """
    V = mesh["V"]
    E = mesh["E"]

    up = read_xy_txt(blade_upper)
    lo = read_xy_txt(blade_lower)

    ncols = len(curved_files_dict)
    fig, axs = plt.subplots(1, ncols, figsize=(6*ncols, 6), constrained_layout=True)

    if ncols == 1:
        axs = [axs]

    for ax, (label, (upper_file, lower_file)) in zip(axs, curved_files_dict.items()):

        upper_edges = read_curved_edges_txt(upper_file)
        lower_edges = read_curved_edges_txt(lower_file)

        # background mesh
        if show_mesh:
            ax.triplot(V[:, 0], V[:, 1], E, linewidth=0.3, color="black", alpha=1, label="linear mesh")

        # reference blade (optional)
        ax.plot(up[:, 0], up[:, 1], "k--", linewidth=1.5, label="reference")
        ax.plot(lo[:, 0], lo[:, 1] + 18.0, "k--", linewidth=1.5)

        # upper curved edges
        # for k, edge in enumerate(upper_edges):
        #     pts = edge["pts"]
        #     xy = evaluate_curved_edge(pts, nplot=nplot_per_edge)

        #     ax.plot(
        #         xy[:, 0], xy[:, 1],
        #         color="tab:blue",
        #         linewidth=2.5,
        #         label="upper" if k == 0 else None
        #     )

        #     if show_edge_nodes:
        #         ax.scatter(
        #             pts[:, 0], pts[:, 1],
        #             s=25,
        #             facecolors="white",
        #             edgecolors="tab:blue",
        #             linewidths=1.5,
        #             zorder=5,
        #             label="upper nodes" if k == 0 else None
        #         )

        # lower curved edges
        for k, edge in enumerate(lower_edges):
            pts = edge["pts"]
            xy = evaluate_curved_edge(pts, nplot=nplot_per_edge)

            ax.plot(
                xy[:, 0], xy[:, 1],
                color="tab:red",
                linewidth=2.5,
                label="curved edge" if k == 0 else None
            )

            if show_edge_nodes:
                ax.scatter(
                    pts[:, 0], pts[:, 1],
                    s=25,
                    facecolors="white",
                    edgecolors="tab:red",
                    linewidths=1.5,
                    zorder=5,
                    label="curved edge nodes" if k == 0 else None
                )

            if show_edge_nodes:
                pts = edge["pts"]
                ax.plot(pts[:,0], pts[:,1], "s", color="tab:red", markersize=4)
        ax.set_xlim(-12, -5)
        ax.set_ylim(14, 20)
        ax.set_title(label, fontsize=28)
        ax.set_aspect("equal")

    axs[0].legend(loc="best")

    fig.suptitle("Curved boundary reconstruction", fontsize=32)
    plt.savefig(out_png, dpi=500, bbox_inches="tight")
    plt.close(fig)

# -----------------------------------------------------------
def read_dg_solution_txt(path: str, ncomp: int = 4):
    """
    Reads DG solution stored as plain text with columns:
        rho, rho*u, rho*v, rho*E

    Returns:
        U : (Nrows, 4) array
    """
    U = np.loadtxt(path, dtype=float)
    if U.ndim == 1:
        U = U.reshape(1, -1)

    if U.shape[1] < ncomp:
        raise ValueError(f"Expected at least {ncomp} columns, got {U.shape[1]} in {path}")

    return U[:, :ncomp]


# -----------------------------------------------------------
def n_lagrange_tri(p: int) -> int:
    """Number of Lagrange nodes on a triangle of order p."""
    return (p + 1) * (p + 2) // 2


# -----------------------------------------------------------
def tri_lagrange_reference_nodes(p: int):
    """
    Reference-triangle Lagrange nodes on:
        (0,0), (1,0), (0,1)

    Ordering:
      p=1 -> [(0,0), (1,0), (0,1)]
      for p>1 -> standard barycentric lattice ordering

    Returns:
        xi : (Np, 2)
    """
    if p == 1:
        return np.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
        ], dtype=float)

    pts = []
    for i in range(p + 1):
        for j in range(p + 1 - i):
            # barycentric: l1 = 1 - xi - eta, l2 = xi, l3 = eta
            xi = j / p
            eta = i / p
            pts.append([xi, eta])
    return np.array(pts, dtype=float)


# -----------------------------------------------------------
def map_reference_to_physical_triangle(xi_eta, xtri):
    """
    Maps reference points (xi,eta) to a physical straight triangle.

    xtri : (3,2) array with physical vertices
           ordered consistently with reference triangle:
           v0=(0,0), v1=(1,0), v2=(0,1)

    Returns:
        X : (N,2) array
    """
    xi = xi_eta[:, 0]
    eta = xi_eta[:, 1]

    N0 = 1.0 - xi - eta
    N1 = xi
    N2 = eta

    X = (
        N0[:, None] * xtri[0][None, :] +
        N1[:, None] * xtri[1][None, :] +
        N2[:, None] * xtri[2][None, :]
    )
    return X


# -----------------------------------------------------------
def build_broken_p1_dg_mesh(mesh, U, p=1, one_based=False):
    """
    Build a broken mesh for DG visualization.

    For p=1:
      - each element contributes its own 3 vertices
      - each element contributes one triangle connectivity
      - values remain discontinuous across elements

    mesh['V'] : (Nn,2)
    mesh['E'] : (Ne,3)
    U         : (Ne*Np, 4), with Np = 3 for p=1

    Returns:
        Xb : (Ne*Np, 2) broken coordinates
        Tb : (Ne, 3) connectivity on broken mesh
        Ub : (Ne*Np, 4) broken solution values
        elem_ids : (Ne,) mapping each broken triangle -> original element id
    """
    if p != 1:
        raise NotImplementedError("This helper currently supports p=1 exactly. See note below for p=2/3 extension.")

    V = mesh["V"]
    E = mesh["E"].copy()

    if one_based:
        E = E - 1

    Np = n_lagrange_tri(p)
    Ne = E.shape[0]

    if U.shape[0] != Ne * Np:
        raise ValueError(
            f"Solution rows ({U.shape[0]}) do not match Ne*Np = {Ne}*{Np} = {Ne*Np}."
        )

    Xb = np.zeros((Ne * Np, 2), dtype=float)
    Ub = np.zeros((Ne * Np, U.shape[1]), dtype=float)
    Tb = np.zeros((Ne, 3), dtype=int)
    elem_ids = np.arange(Ne, dtype=int)

    for e in range(Ne):
        gnodes = E[e]              # global vertex ids of this element
        xtri = V[gnodes, :]        # (3,2)

        i0 = e * Np
        i1 = i0 + Np

        # For p=1, DG nodes coincide with triangle vertices
        Xb[i0:i1, :] = xtri
        Ub[i0:i1, :] = U[i0:i1, :]
        Tb[e, :] = np.array([i0, i0 + 1, i0 + 2], dtype=int)

    return Xb, Tb, Ub, elem_ids


# -----------------------------------------------------------
def plot_dg_solution_fields(
    mesh,
    solution_file,
    out_png="dg_solution_fields.png",
    p=1,
    one_based=False,
    show_mesh=True,
    show_element_ids=False,
    annotate_every=1
):
    """
    Plot 4 conserved-variable fields of DG Euler solution:
        rho, rho*u, rho*v, rho*E

    For p=1, linear interpolation within each element is obtained using
    gouraud shading on the broken mesh.
    """
    U = read_dg_solution_txt(solution_file, ncomp=4)

    Xb, Tb, Ub, elem_ids = build_broken_p1_dg_mesh(mesh, U, p=p, one_based=one_based)

    triang = mtri.Triangulation(Xb[:, 0], Xb[:, 1], Tb)

    labels = [r"$\rho$", r"$\rho u$", r"$\rho v$", r"$\rho E$"]

    fig, axs = plt.subplots(2, 2, figsize=(18, 14), constrained_layout=True)
    axs = axs.ravel()

    for k in range(4):
        ax = axs[k]

        # DG broken-mesh visualization:
        # gouraud -> linear interpolation inside each triangle
        tpc = ax.tripcolor(
            triang,
            Ub[:, k],
            shading="gouraud",
            cmap="turbo"
        )
        cbar = fig.colorbar(tpc, ax=ax, shrink=0.85)
        cbar.set_label(labels[k], fontsize=30)
        cbar.ax.tick_params(labelsize=26)

        if show_mesh:
            # Overlay original mesh for geometric reference
            V = mesh["V"]
            E = mesh["E"].copy()
            if one_based:
                E = E - 1
            ax.triplot(V[:, 0], V[:, 1], E, color="k", linewidth=0.35, alpha=0.65)

        if show_element_ids:
            V = mesh["V"]
            E = mesh["E"].copy()
            if one_based:
                E = E - 1

            for e, conn in enumerate(E):
                if e % annotate_every != 0:
                    continue
                xc = V[conn, 0].mean()
                yc = V[conn, 1].mean()
                ax.text(
                    xc, yc, str(e),
                    fontsize=20,
                    ha="center", va="center",
                    bbox=dict(facecolor="white", edgecolor="none", alpha=0.6, pad=0.5)
                )

        ax.set_title(labels[k], fontsize=30)
        ax.set_aspect("equal")

    fig.suptitle(f"p={p}", fontsize=34)
    fig.savefig(out_png, dpi=500)
    plt.close(fig)


# -----------------------------------------------------------
def inspect_dg_solution_layout(mesh, solution_file, p=1):
    """
    Small helper to verify element-to-solution indexing.
    Prints how rows in the solution file map to element ids.
    """
    U = read_dg_solution_txt(solution_file, ncomp=4)
    Ne = mesh["E"].shape[0]
    Np = n_lagrange_tri(p)

    print(f"Number of elements      : {Ne}")
    print(f"Polynomial order p      : {p}")
    print(f"Lagrange nodes/element  : {Np}")
    print(f"Expected solution rows  : {Ne * Np}")
    print(f"Actual solution rows    : {U.shape[0]}")

    if U.shape[0] != Ne * Np:
        print("WARNING: solution size does not match Ne*Np.")
        return

    for e in range(Ne):
        i0 = e * Np
        i1 = i0 + Np
        print(f"Element {e:4d} -> solution rows [{i0}:{i1})")

# -----------------------------------------------------------
def conserved_to_primitive(U, gamma=1.4):
    """
    U: (N,4) with columns [rho, rho*u, rho*v, rho*E]

    Returns:
        rho, u, v, p, a, M
    """
    rho = U[:, 0]
    rhou = U[:, 1]
    rhov = U[:, 2]
    rhoE = U[:, 3]

    eps = 1.0e-14
    rho = np.maximum(rho, eps)

    u = rhou / rho
    v = rhov / rho
    vel2 = u*u + v*v

    p = (gamma - 1.0) * (rhoE - 0.5 * rho * vel2)
    p = np.maximum(p, eps)

    a = np.sqrt(gamma * p / rho)
    M = np.sqrt(vel2) / np.maximum(a, eps)

    return rho, u, v, p, a, M


# -----------------------------------------------------------
def pressure_coefficient_from_U(U, gamma=1.4, p_out=1.0/1.4, M_out=0.1):
    _, _, _, p, _, _ = conserved_to_primitive(U, gamma=gamma)
    q_out = 0.5 * gamma * p_out * M_out**2
    Cp = (p - p_out) / q_out
    return p, Cp

# -----------------------------------------------------------
# -----------------------------------------------------------
def plot_dg_mach_field(
    mesh,
    solution_file,
    out_png="dg_mach_field.png",
    p=1,
    one_based=False,
    gamma=1.4,
    show_mesh=True,
    show_element_ids=False,
    annotate_every=1
):
    """
    Plot Mach number only.
    """
    U = read_dg_solution_txt(solution_file, ncomp=4)

    Xb, Tb, Ub, elem_ids = build_broken_p1_dg_mesh(mesh, U, p=p, one_based=one_based)
    triang = mtri.Triangulation(Xb[:, 0], Xb[:, 1], Tb)

    _, _, _, _, _, Mach = conserved_to_primitive(Ub, gamma=gamma)

    fig, ax = plt.subplots(1, 1, figsize=(9, 7.5), constrained_layout=True)

    tpc = ax.tripcolor(
        triang,
        Mach,
        shading="gouraud",
        cmap="turbo"
    )

    cbar = fig.colorbar(tpc, ax=ax, shrink=0.88)
    cbar.set_label(r"Ma", fontsize=30)
    cbar.ax.tick_params(labelsize=24)

    if show_mesh:
        V = mesh["V"]
        E = mesh["E"].copy()
        if one_based:
            E = E - 1
        ax.triplot(V[:, 0], V[:, 1], E, color="k", linewidth=0.35, alpha=0.65)

    if show_element_ids:
        V = mesh["V"]
        E = mesh["E"].copy()
        if one_based:
            E = E - 1

        for e, conn in enumerate(E):
            if e % annotate_every != 0:
                continue
            xc = V[conn, 0].mean()
            yc = V[conn, 1].mean()
            ax.text(
                xc, yc, str(e),
                fontsize=16,
                ha="center", va="center",
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.6, pad=0.4)
            )

    ax.set_title(r"Ma", fontsize=30)
    ax.set_aspect("equal")
    fig.savefig(out_png, dpi=500, bbox_inches="tight")
    plt.close(fig)

# -----------------------------------------------------------
def plot_p0_mach(
    mesh,
    solution_file,
    out_png="p0_mach.png",
    gamma=1.4,
    show_mesh=True,
    vmin=None,
    vmax=None
):
    U = read_dg_solution_txt(solution_file, ncomp=4)

    V = mesh["V"]
    E = mesh["E"]

    Ne = E.shape[0]
    if U.shape[0] != Ne:
        raise ValueError(f"Expected {Ne} rows, got {U.shape[0]}")

    # Compute Mach
    _, _, _, _, _, Mach = conserved_to_primitive(U, gamma=gamma)

    triang = mtri.Triangulation(V[:, 0], V[:, 1], E)

    fig, ax = plt.subplots(figsize=(10, 8))

    tpc = ax.tripcolor(
        triang,
        facecolors=Mach,
        shading="flat",
        cmap="turbo",
        vmin=vmin,
        vmax=vmax
    )

    cbar = fig.colorbar(tpc, ax=ax)
    cbar.set_label("Ma", fontsize=30)
    cbar.ax.tick_params(labelsize=24)

    if show_mesh:
        ax.triplot(V[:, 0], V[:, 1], E, color="k", linewidth=0.2, alpha=0.5)

    ax.set_title("Ma", fontsize=30)
    ax.set_aspect("equal")

    fig.savefig(out_png, dpi=500, bbox_inches="tight")
    plt.close(fig)

def get_mach_range(mesh, solution_file, gamma=1.4):
    U = read_dg_solution_txt(solution_file, ncomp=4)
    _, _, _, _, _, Mach = conserved_to_primitive(U, gamma=gamma)
    return Mach.min(), Mach.max()

# -----------------------------------------------------------
def main1():
    base = os.getcwd()

    mesh_file = base + "/mesh_coarse.gri"
    mesh = readgri(mesh_file)

    blade_upper = base + "/bladeupper.txt"
    blade_lower = base + "/bladelower.txt"

    curved_files = {
        "q=1": (
            base + "/upper_curved_edges_Q1_coarse.txt",
            base + "/lower_curved_edges_Q1_coarse.txt",
        ),
        "q=2": (
            base + "/upper_curved_edges_Q2_coarse.txt",
            base + "/lower_curved_edges_Q2_coarse.txt",
        ),
        "q=3": (
            base + "/upper_curved_edges_Q3_coarse.txt",
            base + "/lower_curved_edges_Q3_coarse.txt",
        ),
    }

    out_png = base + "/mesh_curved.png"

    plot_curved_boundary_multiQ(
        mesh,
        blade_upper,
        blade_lower,
        curved_files,
        out_png,
        nplot_per_edge=300,
        show_mesh=True,
        show_edge_nodes=True
    )
def main2():
    base = Path(__file__).resolve().parent

    mesh_file = base / "mesh_test.gri"
    sol_file  = base / "results" / "test0322.txt"

    out_png_cons = base / "dg_solution_fields_p1_test.png"
    out_png_mach = base / "dg_mach_field_p1_test.png"

    mesh = readgri(str(mesh_file))

    inspect_dg_solution_layout(mesh, str(sol_file), p=1)

    plot_dg_solution_fields(
        mesh,
        str(sol_file),
        out_png=str(out_png_cons),
        p=1,
        one_based=False,
        show_mesh=True,
        show_element_ids=False
    )

    plot_dg_mach_field(
        mesh,
        str(sol_file),
        out_png=str(out_png_mach),
        p=1,
        one_based=False,
        gamma=1.4,
        show_mesh=True,
        show_element_ids=False
    )

def main3():
    base = Path(__file__).resolve().parent

    mesh_file = base / "mesh_refined_2394.gri"

    roe_file  = base / "results" / "coarse_mesh_steady_p0_q1_RK4_roe.txt"
    hlle_file = base / "results" / "coarse_mesh_steady_p0_q1_RK4_hlle.txt"

    mesh = readgri(str(mesh_file))

    # ---- compute global Mach range ----
    mmin1, mmax1 = get_mach_range(mesh, str(roe_file))
    mmin2, mmax2 = get_mach_range(mesh, str(hlle_file))

    vmin = min(mmin1, mmin2)
    vmax = max(mmax1, mmax2)

    print(f"Global Mach range: [{vmin:.4f}, {vmax:.4f}]")

    # ---- plot both with SAME scale ----
    plot_p0_mach(
        mesh,
        str(roe_file),
        out_png=str(base / "mach_p0_roe.png"),
        vmin=vmin,
        vmax=vmax
    )

    plot_p0_mach(
        mesh,
        str(hlle_file),
        out_png=str(base / "mach_p0_hlle.png"),
        vmin=vmin,
        vmax=vmax
    )

if __name__ == "__main__":
    main1()
    # main2()
    main3()