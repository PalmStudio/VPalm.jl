module VPalm

# For the random number generator:
import Random

# For managing the MTG:
import MultiScaleTreeGraph: Node, NodeMTG, traverse!, symbol

# IO:
import YAML, OrderedCollections

# For the 3D:
import PlantGeom
import Meshes
import TransformsBase: →
import Rotations: RotY, RotZ
import PlyIO

include("utils.jl")
include("IO/parameters_IO.jl")

# Entry point:
include("architecture/mtg_skeleton.jl")

# Allometries:
include("allometries/stem.jl")
include("allometries/internode.jl")
include("allometries/leaf.jl")

# Architecture:
include("architecture/compute_properties_stem.jl")
include("architecture/compute_properties_internode.jl")
include("architecture/compute_properties_leaf.jl")
include("architecture/compute_properties_petiole.jl")

# Geometry:
include("geometry/read_ply.jl")
include("geometry/snag.jl")
include("geometry/cylinder.jl")
include("geometry/add_geometry.jl")
include("geometry/petiole.jl")

export read_parameters, write_parameters
end
