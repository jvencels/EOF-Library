import gmsh
import sys


def main(out_msh: str):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("elmer_em_3d_quasi2d_safe")

    # ----------------------------
    # Geometry parameters
    # ----------------------------
    L, H, W = 1.0, 0.1, 0.1
    Nx, Ny, Nz = 40, 10, 10

    dom = max(L, H, W)
    # Be a bit generous; OCC bbox values can be slightly fuzzy
    tol = 1e-6 * dom

    occ = gmsh.model.occ

    # ----------------------------
    # Base YZ surface at x = 0
    # ----------------------------
    p1 = occ.addPoint(0.0, 0.0, 0.0)
    p2 = occ.addPoint(0.0, H,   0.0)
    p3 = occ.addPoint(0.0, H,   W)
    p4 = occ.addPoint(0.0, 0.0, W)

    l1 = occ.addLine(p1, p2)  # along +y
    l2 = occ.addLine(p2, p3)  # along +z
    l3 = occ.addLine(p3, p4)  # along -y
    l4 = occ.addLine(p4, p1)  # along -z

    cl = occ.addCurveLoop([l1, l2, l3, l4])
    s0 = occ.addPlaneSurface([cl])

    # ----------------------------
    # Extrude along X
    # ----------------------------
    ext = occ.extrude([(2, s0)], L, 0.0, 0.0, numElements=[Nx], recombine=True)
    occ.synchronize()

    # Extract volume
    vols = [tag for dim, tag in ext if dim == 3]
    assert len(vols) == 1, f"Expected 1 volume from extrude, got {len(vols)}"
    vol = vols[0]

    # ----------------------------
    # Get boundary surfaces of the volume (THIS is the robust list)
    # ----------------------------
    bnd = gmsh.model.getBoundary([(3, vol)], oriented=False, recursive=False)
    surf_tags = [tag for dim, tag in bnd if dim == 2]
    assert len(surf_tags) >= 6, f"Expected at least 6 boundary faces, got {len(surf_tags)}"

    # ----------------------------
    # Identify inlet/outlet + Y/Z walls by bbox planes
    # ----------------------------
    inlet = None
    outlet = None
    y0, yH, z0, zW = [], [], [], []

    for s in surf_tags:
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, s)

        # Planes x=0 and x=L
        if abs(xmin - 0.0) < tol and abs(xmax - 0.0) < tol:
            inlet = s
            continue
        if abs(xmin - L) < tol and abs(xmax - L) < tol:
            outlet = s
            continue

        # Planes y=0 and y=H
        if abs(ymin - 0.0) < tol and abs(ymax - 0.0) < tol:
            y0.append(s)
        elif abs(ymin - H) < tol and abs(ymax - H) < tol:
            yH.append(s)

        # Planes z=0 and z=W
        if abs(zmin - 0.0) < tol and abs(zmax - 0.0) < tol:
            z0.append(s)
        elif abs(zmin - W) < tol and abs(zmax - W) < tol:
            zW.append(s)

    assert inlet is not None, "Failed to identify inlet surface (x=0 plane)"
    assert outlet is not None, "Failed to identify outlet surface (x=L plane)"
    assert len(y0) == 1, f"Expected 1 y=0 face, got {len(y0)}: {y0}"
    assert len(yH) == 1, f"Expected 1 y=H face, got {len(yH)}: {yH}"
    assert len(z0) == 1, f"Expected 1 z=0 face, got {len(z0)}: {z0}"
    assert len(zW) == 1, f"Expected 1 z=W face, got {len(zW)}: {zW}"

    # ----------------------------
    # Physical groups (DETERMINISTIC IDS)
    # ----------------------------
    gmsh.model.addPhysicalGroup(2, [inlet], tag=20)
    gmsh.model.setPhysicalName(2, 20, "InletX")

    gmsh.model.addPhysicalGroup(2, [outlet], tag=21)
    gmsh.model.setPhysicalName(2, 21, "OutletX")

    gmsh.model.addPhysicalGroup(2, y0, tag=10)
    gmsh.model.setPhysicalName(2, 10, "FaradayMinusY")

    gmsh.model.addPhysicalGroup(2, yH, tag=11)
    gmsh.model.setPhysicalName(2, 11, "FaradayPlusY")

    gmsh.model.addPhysicalGroup(2, z0 + zW, tag=30)
    gmsh.model.setPhysicalName(2, 30, "SideWallsZ")

    gmsh.model.addPhysicalGroup(3, [vol], tag=1)
    gmsh.model.setPhysicalName(3, 1, "fluid")

    # ----------------------------
    # Transfinite meshing
    # ----------------------------
    gmsh.model.mesh.setTransfiniteCurve(l1, Ny + 1)
    gmsh.model.mesh.setTransfiniteCurve(l3, Ny + 1)
    gmsh.model.mesh.setTransfiniteCurve(l2, Nz + 1)
    gmsh.model.mesh.setTransfiniteCurve(l4, Nz + 1)

    # Curves created by extrusion that run along x will have bbox extent ~L
    for dim, c in gmsh.model.getEntities(1):
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(1, c)
        if abs((xmax - xmin) - L) < 10 * tol:
            gmsh.model.mesh.setTransfiniteCurve(c, Nx + 1)

    for s in surf_tags:
        gmsh.model.mesh.setTransfiniteSurface(s)
        gmsh.model.mesh.setRecombine(2, s)

    gmsh.model.mesh.setTransfiniteVolume(vol)

    # ----------------------------
    # Mesh + sanity check
    # ----------------------------
    gmsh.model.mesh.generate(3)

    bnd2 = gmsh.model.getBoundary([(3, vol)], oriented=False, recursive=True)
    bnd_surfs = [tag for dim, tag in bnd2 if dim == 2]
    missing = [s for s in bnd_surfs if not gmsh.model.getPhysicalGroupsForEntity(2, s)]
    assert len(missing) == 0, f"Missing boundary surfaces (no physical group): {missing}"

    # Global bbox print
    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    xs = coords[0::3]
    ys = coords[1::3]
    zs = coords[2::3]
    print(
        f"[gmsh] GLOBAL bbox: "
        f"({min(xs):.6g} {min(ys):.6g} {min(zs):.6g}) "
        f"({max(xs):.6g} {max(ys):.6g} {max(zs):.6g})"
    )

    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.write(out_msh)

    gmsh.finalize()


if __name__ == "__main__":
    main(sys.argv[1])
