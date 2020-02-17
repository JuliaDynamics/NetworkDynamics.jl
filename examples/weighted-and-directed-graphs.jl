# You will find a step-by-step guide to this example in the docs and the
# corresponding jupyter notebook on our github repository.

using DelimitedFiles

G = readdlm("git/NetworkDynamics.jl/examples/Norm_G_DTI.txt", ',', Float64, '\n');
