#!/usr/bin/env python3
"""
Convert glacier NetCDF data to printable STL files.

This script creates:
1. A solid topography model with a flat base and vertical walls.
2. A watertight glacier model using the glacier surface and bed topography.

Expected NetCDF variables:
- x
- y
- thk
- usurf
"""

from pathlib import Path
from collections import defaultdict

import numpy as np
from netCDF4 import Dataset
from scipy.spatial import Delaunay
from skimage import measure
from shapely.geometry import Polygon, Point
from stl import mesh


def grid_to_triangles(nx, ny, offset=0):
    """Create triangular faces for a regular grid."""
    for i in range(ny - 1):
        for j in range(nx - 1):
            idx = i * nx + j + offset
            yield [idx, idx + 1, idx + nx]
            yield [idx + 1, idx + nx + 1, idx + nx]


def save_stl(vertices, faces, output_path):
    """Save vertices and triangular faces as STL."""
    stl_mesh = mesh.Mesh(np.zeros(len(faces), dtype=mesh.Mesh.dtype))

    for i, face in enumerate(faces):
        for j in range(3):
            stl_mesh.vectors[i][j] = vertices[face[j]]

    stl_mesh.save(str(output_path))
    print(f"Saved: {output_path}")


def create_topography_stl(ds, output_path, vertical_offset=100):
    """Create a solid topography STL with a flat base."""
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

    top_vertices = np.column_stack(
        (X.flatten(), Y.flatten(), Z)
    ).astype(np.float32)

    bottom_vertices = np.column_stack(
        (X.flatten(), Y.flatten(), np.zeros_like(Z))
    ).astype(np.float32)

    top_faces = list(grid_to_triangles(nx, ny, offset=0))
    bottom_faces = list(grid_to_triangles(nx, ny, offset=len(top_vertices)))

    wall_faces = []

    # Left and right walls
    for i in range(ny - 1):
        # Left edge
        top_a = i * nx
        top_b = (i + 1) * nx
        bot_a = top_a + len(top_vertices)
        bot_b = top_b + len(top_vertices)

        wall_faces.append([top_a, bot_a, top_b])
        wall_faces.append([top_b, bot_a, bot_b])

        # Right edge
        top_a = i * nx + (nx - 1)
        top_b = (i + 1) * nx + (nx - 1)
        bot_a = top_a + len(top_vertices)
        bot_b = top_b + len(top_vertices)

        wall_faces.append([top_a, top_b, bot_a])
        wall_faces.append([top_b, bot_b, bot_a])

    # Front and back walls
    for j in range(nx - 1):
        # Bottom edge
        top_a = j
        top_b = j + 1
        bot_a = top_a + len(top_vertices)
        bot_b = top_b + len(top_vertices)

        wall_faces.append([top_a, bot_a, top_b])
        wall_faces.append([top_b, bot_a, bot_b])

        # Top edge
        top_a = (ny - 1) * nx + j
        top_b = (ny - 1) * nx + j + 1
        bot_a = top_a + len(top_vertices)
        bot_b = top_b + len(top_vertices)

        wall_faces.append([top_a, top_b, bot_a])
        wall_faces.append([top_b, bot_b, bot_a])

    vertices = np.vstack((top_vertices, bottom_vertices))
    faces = np.asarray(top_faces + bottom_faces + wall_faces)

    save_stl(vertices, faces, output_path)


def get_boundary_edges(triangles):
    """Return all triangle edges that occur only once."""
    edge_count = defaultdict(int)

    for tri in triangles:
        for i in range(3):
            edge = tuple(sorted((tri[i], tri[(i + 1) % 3])))
            edge_count[edge] += 1

    return [edge for edge, count in edge_count.items() if count == 1]


def filter_by_edge_length(vertices, triangles, max_edge_length):
    """Remove triangles with edges longer than max_edge_length."""
    valid_triangles = []

    for tri in triangles:
        a, b, c = vertices[tri]
        edges = [
            np.linalg.norm(a[:2] - b[:2]),
            np.linalg.norm(b[:2] - c[:2]),
            np.linalg.norm(c[:2] - a[:2]),
        ]

        if max(edges) < max_edge_length:
            valid_triangles.append(tri)

    return np.asarray(valid_triangles)


def create_glacier_stl(
    ds,
    output_path,
    min_thickness=1.0,
    max_edge_length=200.0,
):
    """Create a watertight glacier STL from surface and bed elevation."""
    thk = np.asarray(ds.variables["thk"][:])
    usurf = np.asarray(ds.variables["usurf"][:])
    topg = usurf - thk

    x = np.asarray(ds.variables["x"][:])
    y = np.asarray(ds.variables["y"][:])

    X, Y = np.meshgrid(x, y)
    icemask = thk > min_thickness

    contours = measure.find_contours(icemask, 0.5)

    if len(contours) == 0:
        raise ValueError("No glacier outline found. Check the thickness threshold.")

    main_contour = max(contours, key=len)

    contour_y = np.interp(main_contour[:, 0], np.arange(len(y)), y)
    contour_x = np.interp(main_contour[:, 1], np.arange(len(x)), x)

    glacier_polygon = Polygon(zip(contour_x, contour_y))

    X_flat = X.flatten()
    Y_flat = Y.flatten()
    usurf_flat = usurf.flatten()
    topg_flat = topg.flatten()
    mask_flat = icemask.flatten()

    points_xy = np.column_stack((X_flat, Y_flat))

    inside_polygon = np.asarray(
        [glacier_polygon.contains(Point(point)) for point in points_xy]
    )

    valid = mask_flat & inside_polygon

    x_valid = X_flat[valid]
    y_valid = Y_flat[valid]
    z_top = usurf_flat[valid]
    z_bottom = topg_flat[valid]

    points_2d = np.column_stack((x_valid, y_valid))

    tri = Delaunay(points_2d)
    triangles = tri.simplices

    vertices_top = np.column_stack((x_valid, y_valid, z_top))
    vertices_bottom = np.column_stack((x_valid, y_valid, z_bottom))

    triangles_top = filter_by_edge_length(
        vertices_top, triangles, max_edge_length
    )

    triangles_bottom = filter_by_edge_length(
        vertices_bottom, triangles, max_edge_length
    )

    boundary_edges = get_boundary_edges(triangles_top)

    wall_triangles = []
    n_vertices = len(vertices_top)

    for a, b in boundary_edges:
        a_bottom = a + n_vertices
        b_bottom = b + n_vertices

        wall_triangles.append([a, a_bottom, b])
        wall_triangles.append([b, a_bottom, b_bottom])

    triangles_bottom_shifted = triangles_bottom + n_vertices

    vertices = np.vstack((vertices_top, vertices_bottom))
    faces = np.vstack(
        (
            triangles_top,
            triangles_bottom_shifted,
            np.asarray(wall_triangles),
        )
    )

    save_stl(vertices, faces, output_path)


def main():
    input_file = Path("Aletsch.nc")
    output_dir = Path("STLS")
    output_dir.mkdir(exist_ok=True)

    with Dataset(input_file) as ds:
        print("Available variables:")
        print(ds.variables.keys())

        create_topography_stl(
            ds,
            output_dir / "topography.stl",
            vertical_offset=100,
        )

        create_glacier_stl(
            ds,
            output_dir / "glacier.stl",
            min_thickness=1.0,
            max_edge_length=200.0,
        )


if __name__ == "__main__":
    main()