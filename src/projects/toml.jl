dmap(f, d::Dict) = Dict(k => f(v) for (k, v) in d)

# Function to load project configuration from TOML
function load_project_config(toml)
    config = TOML.parsefile(toml)

    # First load all datasets
    datasets = dmap(LDataSet, get(config, "datasets", Dict()))

    # Process instruments and associate datasets with them
    instruments = dmap(get(config, "instruments", Dict())) do v
        dataset_refs = get(v, "datasets", String[])
        v["datasets"] = filter(x -> in(x.first, dataset_refs), datasets)
        Instrument(v)
    end

    project = Project(;
        name=config["name"],
        metadata=get(config, "metadata", Dict()),
        instruments,
        datasets
    )

    dict = Dict{Symbol,Any}()
    dict[Symbol(lowercase(abbr(project)))] = project
    for (key, value) in pairs(datasets) ∪ pairs(instruments)
        dict[Symbol(key)] = value
    end

    return dict
end