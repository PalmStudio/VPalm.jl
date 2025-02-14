module VPalm

# For the random number generator:
import Random

# For managing the MTG:
import MultiScaleTreeGraph: Node, NodeMTG, traverse!, symbol, reparent!

# IO:
import YAML, OrderedCollections

# For the 3D:
import PlantGeom
import Meshes
import TransformsBase: →
import Rotations: RotX, RotY, RotZ, RotXY
import Rotations
import PlyIO

# For the biomechanical model
import Interpolations: linear_interpolation

include("utils.jl")
include("IO/parameters_IO.jl")

# Entry point:
include("architecture/mtg_skeleton.jl")

# Biomechanical model
include("biomechanic/xyz_to_angles.jl")
include("biomechanic/angles_to_xyz.jl")
include("biomechanic/rotate_yz.jl")
include("biomechanic/rota_inverse_yz.jl")
include("biomechanic/inertia_flex_rota.jl")
include("biomechanic/interpolate_points.jl")
include("biomechanic/bend.jl")

# Allometries:
include("allometries/stem.jl")
include("allometries/internode.jl")
include("allometries/leaf.jl")
include("allometries/petiole.jl")
include("allometries/rachis.jl")

# Architecture:
include("architecture/compute_properties_stem.jl")
include("architecture/compute_properties_internode.jl")
include("architecture/compute_properties_leaf.jl")
include("architecture/compute_properties_petiole.jl")
include("architecture/compute_properties_rachis.jl")

# Geometry:
include("geometry/read_ply.jl")
include("geometry/snag.jl")
include("geometry/cylinder.jl")
include("geometry/add_geometry.jl")
include("geometry/petiole.jl")

# Instance (create an organ with architecture + geometry)
include("instance/petiole.jl")
include("instance/rachis.jl")

export read_parameters, write_parameters
end
