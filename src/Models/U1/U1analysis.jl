
get_symmetry(::Type{U1PionCorrelator}) = SymmetricCorrelator
correlator_fit_function(::Type{U1PionCorrelator}, corrws::AbstractCorrelatorAnalysis) = SymCorrelator(T=corrws.T, s=1, c=0)


get_symmetry(::Type{U1PCACCorrelator}) = AntisymmetricCorrelator
correlator_fit_function(::Type{U1PCACCorrelator}, corrws::AbstractCorrelatorAnalysis) = SymCorrelator(T=corrws.T, s=1, c=0)
