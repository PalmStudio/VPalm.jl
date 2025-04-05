module VPalm

# For the random number generator:
import Random

# For managing the MTG:
import MultiScaleTreeGraph: Node, NodeMTG, traverse!, symbol, reparent!, addchild!, descendants, delete_nodes!

# IO:
import YAML, OrderedCollections

# For the 3D:
import PlantGeom
import Meshes
import TransformsBase: →
import Rotations: RotX, RotY, RotZ, RotYZ, RotXYZ, RotZY, RotYZX, RotZYX
import Rotations
import PlyIO
import Unitful: @u_str, ustrip, unit, NoUnits, uconvert, Quantity

# For the biomechanical model
import Interpolations: linear_interpolation

include("units.jl")
include("utils.jl")
include("IO/parameters_IO.jl")

# Entry point:
include("architecture/mtg_skeleton.jl")

# Biomechanical models
include("biomechanics/complete/xyz_dist_angles.jl")
include("biomechanics/complete/inertia_flex_rota.jl")
include("biomechanics/complete/interpolate_points.jl")
include("biomechanics/complete/bend.jl")
include("biomechanics/complete/unbend.jl")
include("biomechanics/simplified/young_modulus.jl")

# Allometries:
include("allometries/stem.jl")
include("allometries/internode.jl")
include("allometries/leaf.jl")
include("allometries/petiole.jl")

# Architecture:
include("architecture/leaf_rank.jl")
include("architecture/compute_properties_stem.jl")
include("architecture/compute_properties_internode.jl")
include("architecture/compute_properties_leaf.jl")
include("architecture/compute_properties_petiole.jl")
include("architecture/compute_properties_rachis.jl")

# Geometry:
include("geometry/read_ply.jl")
include("geometry/snag.jl")
include("geometry/cylinder.jl")
include("geometry/plane.jl")
include("geometry/add_geometry.jl")
include("geometry/sections.jl")
include("geometry/leaflets.jl")

# Instance (create an organ with architecture + geometry)
include("instance/petiole.jl")
include("instance/rachis.jl")
include("instance/leaflets.jl")

export read_parameters, write_parameters
end
