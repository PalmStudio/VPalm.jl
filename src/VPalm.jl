module VPalm

# For the random number generator:
import Random

# For managing the MTG:
import MultiScaleTreeGraph: Node, NodeMTG

# IO:
import YAML, OrderedCollections

# For the 3D:
import PlantGeom
import Meshes

include("IO/parameters_IO.jl")
include("architecture/mtg_skeleton.jl")
include("allometries/stem.jl")
include("allometries/internode.jl")
include("allometries/leaf.jl")
include("static_mockup.jl")

export read_parameters, write_parameters
export static_mockup
end
