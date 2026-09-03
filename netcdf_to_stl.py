#!/usr/bin/env python3
"""
Convert glacier NetCDF data to printable STL files.

Creates:
1. Topography STL (solid base terrain)
2. Glacier STL (watertight ice volume)
"""

from pathlib import Path
from collections import defaultdict
import argparse

import numpy as np
from netCDF4 import Dataset
from scipy.spatial import Delaunay
from skimage import measure
from shapely.geometry import Polygon, Point
from stl import mesh


def grid_to_triangles(nx, ny, offset=0):
    for i in range(ny - 1):
        for j in range(nx - 1):
            idx = i * nx + j + offset
            yield [idx, idx + 1, idx + nx]
            yield [idx + 1, idx + nx + 1, idx + nx]


def save_stl(vertices, faces, output_path):
    stl_mesh = mesh.Mesh(np.zeros(len(faces), dtype=mesh.Mesh.dtype))

    for i, face in enumerate(faces):
        for j in range(3):
            stl_mesh.vectors[i][j] = vertices[face[j]]

    stl_mesh.save(str(output_path))
    print(f"Saved: {output_path}")


def create_topography_stl(ds, output_path, vertical_offset=100):
    thk = np.asarray(ds.variables["thk"][:])
    usurf = np.asarray(ds.variables["usurf"][:])

    bedrock = usurf - thk

    x = np.asarray(ds.variables["x"][:])
    y = np.asarray(ds.variables["y"][:])
    ny, nx = bedrock.shape

    X, Y = np.meshgrid(x, y)

    Z = bedrock.flatten()
    Z -= np.nanmin(Z)
    Z += vertical_offset

    top_vertices = np.column_stack((X.flatten(), Y.flatten(), Z)).astype(np.float32)
    bottom_vertices = np.column_stack((X.flatten(), Y.flatten(), np.zeros_like(Z))).astype(np.float32)

    top_faces = list(grid_to_triangles(nx, ny, offset=0))
    bottom_faces = list(grid_to_triangles(nx, ny, offset=len(top_vertices)))

    wall_faces = []

    # vertical walls
    for i in range(ny - 1):
        # left
        a, b = i * nx, (i + 1) * nx
        wall_faces += [[a, a + len(top_vertices), b],
                       [b, a + len(top_vertices), b + len(top_vertices)]]

        # right
        a, b = i * nx + (nx - 1), (i + 1) * nx + (nx - 1)
        wall_faces += [[a, b, a + len(top_vertices)],
                       [b, b + len(top_vertices), a + len(top_vertices)]]

    for j in range(nx - 1):
        # bottom
        a, b = j, j + 1
        wall_faces += [[a, a + len(top_vertices), b],
                       [b, a + len(top_vertices), b + len(top_vertices)]]

        # top
        a, b = (ny - 1) * nx + j, (ny - 1) * nx + j + 1
        wall_faces += [[a, b, a + len(top_vertices)],
                       [b, b + len(top_vertices), a + len(top_vertices)]]

    vertices = np.vstack((top_vertices, bottom_vertices))
    faces = np.asarray(top_faces + bottom_faces + wall_faces)

    save_stl(vertices, faces, output_path)


def get_boundary_edges(triangles):
    edge_count = defaultdict(int)

    for tri in triangles:
        for i in range(3):
            edge = tuple(sorted((tri[i], tri[(i + 1) % 3])))
            edge_count[edge] += 1

    return [e for e, c in edge_count.items() if c == 1]


def filter_by_edge_length(vertices, triangles, max_edge_length):
    valid = []
    for tri in triangles:
        a, b, c = vertices[tri]
        edges = [
            np.linalg.norm(a[:2] - b[:2]),
            np.linalg.norm(b[:2] - c[:2]),
            np.linalg.norm(c[:2] - a[:2]),
        ]
        if max(edges) < max_edge_length:
            valid.append(tri)

    return np.asarray(valid, dtype=np.int64).reshape(-1, 3)


def create_glacier_stl(ds, output_path, min_thickness=1.0, max_edge_length=2000.0):

    thk = np.asarray(ds.variables["thk"][:])
    usurf = np.asarray(ds.variables["usurf"][:])
    topg = usurf - thk

    x = np.asarray(ds.variables["x"][:])
    y = np.asarray(ds.variables["y"][:])

    X, Y = np.meshgrid(x, y)
    mask = thk > min_thickness

    contours = measure.find_contours(mask, 0.5)
    if not contours:
        raise ValueError("No glacier found")

    main = max(contours, key=len)

    cy = np.interp(main[:, 0], np.arange(len(y)), y)
    cx = np.interp(main[:, 1], np.arange(len(x)), x)

    poly = Polygon(zip(cx, cy))

    Xf, Yf = X.flatten(), Y.flatten()
    inside = np.array([poly.contains(Point(p)) for p in zip(Xf, Yf)])

    valid = mask.flatten() & inside

    x_v = Xf[valid]
    y_v = Yf[valid]
    z_top = usurf.flatten()[valid] + 50
    z_bot = topg.flatten()[valid]

    pts2d = np.column_stack((x_v, y_v))

    if len(pts2d) < 3:
        raise ValueError(
            f"Not enough glacier points to triangulate. "
            f"Found {len(pts2d)} points. Try lowering --min-thickness."
        )

    tri = Delaunay(pts2d).simplices

    v_top = np.column_stack((x_v, y_v, z_top))
    v_bot = np.column_stack((x_v, y_v, z_bot))

    tri_top = filter_by_edge_length(v_top, tri, max_edge_length)
    tri_bot = filter_by_edge_length(v_bot, tri, max_edge_length)

    if len(tri_top) == 0:
        raise ValueError(
            f"All glacier triangles were removed by --max-edge-length={max_edge_length}. "
            f"Try increasing --max-edge-length."
        )

    boundary = get_boundary_edges(tri_top)

    n = len(v_top)
    walls = []

    for a, b in boundary:
        walls += [[a, a + n, b],
                  [b, a + n, b + n]]

    wall_faces = np.asarray(walls, dtype=np.int64).reshape(-1, 3)

    vertices = np.vstack((v_top, v_bot))
    faces = np.vstack((tri_top, tri_bot + n, wall_faces))

    save_stl(vertices, faces, output_path)


def main():
    parser = argparse.ArgumentParser(description="Convert glacier NetCDF to STL")

    parser.add_argument("input", type=str, help="Input NetCDF file")
    parser.add_argument("--output-dir", type=str, default="STLS", help="Output directory")

    parser.add_argument("--vertical-offset", type=float, default=100.0)
    parser.add_argument("--min-thickness", type=float, default=1.0)
    parser.add_argument("--max-edge-length", type=float, default=100.0)

    args = parser.parse_args()

    input_file = Path(args.input)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)

    with Dataset(input_file) as ds:
        print("Variables:", list(ds.variables.keys()))

        create_topography_stl(
            ds,
            output_dir / "topography.stl",
            vertical_offset=args.vertical_offset,
        )

        create_glacier_stl(
            ds,
            output_dir / "glacier.stl",
            min_thickness=args.min_thickness,
            max_edge_length=args.max_edge_length,
        )


if __name__ == "__main__":
    main()