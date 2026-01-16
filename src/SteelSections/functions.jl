#get all tabulated sections
allW() = [W(name) for name in names[Wrange]]
allC() = [C(name) for name in names[Crange]]
allL() = [L(name) for name in names[Lrange]]
allLL() = [LL(name) for name in names[LLrange]]
allWT() = [WT(name) for name in names[WTrange]]
allHSSRect() = [HSSRect(name) for name in names[HSSRectrange]]
allHSSRound() = [HSSRound(name) for name in names[HSSRoundrange]]

#get all names
Wnames = names[Wrange]
Cnames = names[Crange]
Lnames = names[Lrange]
LLnames = names[LLrange]
WTnames = names[WTrange]
HSSRectnames = names[HSSRectrange]
HSSRoundnames = names[HSSRoundrange]

# cache access to sections
const SECTION_CACHE = Dict{String, AbstractSection}()

# populate cache (fuzzy matching helper)
clean_name(n) = replace(uppercase(n), " " => "")

function populate_cache!()
    empty!(SECTION_CACHE)
    
    function cache_section!(s)
        SECTION_CACHE[clean_name(s.name)] = s
        SECTION_CACHE[clean_name(s.name_imperial)] = s
    end

    foreach(cache_section!, allW())
    foreach(cache_section!, allC())
    foreach(cache_section!, allL())
    foreach(cache_section!, allLL())
    foreach(cache_section!, allWT())
    foreach(cache_section!, allHSSRect())
    foreach(cache_section!, allHSSRound())
    
    return nothing
end

"""
get a steel section by name (metric or imperial)
is case-insensitive and space-insensitive.
"""
function get_section(name::String)
    # auto-populates the cache on first call.
    if isempty(SECTION_CACHE)
        populate_cache!()
    end
    
    key = clean_name(name)
    if haskey(SECTION_CACHE, key)
        return SECTION_CACHE[key]
    else
        error("Section '$name' (cleaned: '$key') not found in AISC database.")
    end
end

# convert units using Unitful
function toASAPframe(section::TorsionAllowed, E::Quantity, G::Quantity; unit = u"mm")
    return Section(
        Float64(ustrip(unit^2, section.A)),
        Float64(ustrip(u"N"/unit^2, E)),
        Float64(ustrip(u"N"/unit^2, G)),
        Float64(ustrip(unit^4, section.Ix)),
        Float64(ustrip(unit^4, section.Iy)),
        Float64(ustrip(unit^4, section.J))
    )
end

# only using section name
function toASAPframe(name::String, E::Quantity, G::Quantity; unit = u"mm")
    return toASAPframe(get_section(name), E, G; unit = unit)
end

# including standard defaults
function toASAPframe(name::String; E = 200u"GPa", G = 77u"GPa", unit = u"mm")
    return toASAPframe(name, E, G; unit = unit)
end

# fallback for Real
function toASAPframe(section::TorsionAllowed, E::Real, G::Real; unit = u"mm")
    return toASAPframe(section, E * u"N/mm^2", G * u"N/mm^2"; unit = unit)
end

function toASAPtruss(section::AbstractSection, E::Quantity; unit = u"mm")
    return TrussSection(
        Float64(ustrip(unit^2, section.A)),
        Float64(ustrip(u"N"/unit^2, E))
    )
end

# fallback for Real
function toASAPtruss(section::AbstractSection, E::Real; unit = u"mm")
    return toASAPtruss(section, E * u"N/mm^2"; unit = unit)
end
