@doc project_doc("THEMIS", "Time History of Events and Macroscale Interactions during Substorms")
module THEMIS
using ..SPEDAS: @load_project_config

@load_project_config "THEMIS.toml"

end