
using Meshes, GeoIO, Rotations, GLMakie

# Example OPF of a palm tree:
opf = read_opf("/Users/rvezy/Documents/dev/VPalm_old/VPalm_Xpalm_coupling/new_format_biomass.opf")
# You can generate one using this project: https://github.com/PalmStudio/Vpalm

# Extract the reference meshes from the OPF:
ref_meshes = get_ref_meshes(opf)
# Extract a mesh of a snag:
snag_mesh = ref_meshes.meshes[2].mesh

# Normalize it:

# Rotate it:
new_mesh = snag_mesh |> Meshes.Rotate(AngleAxis(deg2rad(-35), 0.0, 1.0, 0.0))

# Translate it:
min_point = Meshes.boundingbox(new_mesh).min
zero_m = zero(coords(min_point).x)
new_mesh = new_mesh |> Meshes.Translate(-coords(min_point).x, zero_m, zero_m)
mesh_size = Meshes.boundingbox(new_mesh).max - Meshes.boundingbox(new_mesh).min
new_mesh = new_mesh |> Scale([1.0 / i.val for i in mesh_size.coords]...)

GeoIO.save("assets/snag.ply", GeoIO.georef(nothing, new_mesh))
# We open it in blender to remove duplicated vertices. To do so, we need to import it with Z axis forward and -Y up.
mesh_ = GeoIO.load("assets/snag.ply")

viz(mesh_.geometry)