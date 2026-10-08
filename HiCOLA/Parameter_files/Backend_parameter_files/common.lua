------------------------------------------------------------
-- Shared helpers for the comparison runs in this folder.
--
-- Layout:  <target>/base.lua                     settings shared by all runs matched to one
--                                                reference simulation (box, cosmology, ICs)
--          <target>/<model>_<realization>.lua    one run: sets `model` and `realization`,
--                                                then loads base.lua
--
-- model        LCDM | nDGP_N<rcH0>               e.g. nDGP_N1.2 -> DGP with rc*H0/c = 1.2
-- realization  s<seed>[-fixed[-inv]]             -fixed     : ic_fix_amplitude = true
--                                                -fixed-inv : ic_fix_amplitude = true,
--                                                             ic_reverse_phases = true
--
-- The run is written to $HICOLA_RUN_ROOT/<target>/<model>/<realization>/ with
-- simulation_name = <target>_<model>_<realization>.
-- HICOLA_RUN_ROOT defaults to /scratch/$USER/Hi-COLA_Output (see jobs/hicola_env.sh).
------------------------------------------------------------

local function absolute(path)
  if path:sub(1, 1) == "/" then return path end
  return (os.getenv("PWD") or ".") .. "/" .. path
end

-- Folder of this file, and the Hi-COLA repository it lives in
hicola_param_dir = absolute(debug.getinfo(1, "S").source:sub(2)):match("^(.*)/[^/]*$")
hicola_root = os.getenv("HICOLA_ROOT")
  or (hicola_param_dir:gsub("/HiCOLA/Parameter_files/Backend_parameter_files$", ""))
-- Frontend tables (expansion/force files), at their original Hi-COLA location
hicola_output_dir = hicola_root .. "/HiCOLA/Output"
-- Simulation outputs (snapshots, P(k), gravity model files)
hicola_run_root = os.getenv("HICOLA_RUN_ROOT")
  or ("/scratch/" .. (os.getenv("USER") or "unknown") .. "/Hi-COLA_Output")

function hicola_run(target, model, realization)
  assert(model and realization, "set `model` and `realization` before loading base.lua")
  local run = { target = target, model = model, realization = realization }

  if model == "LCDM" then
    run.gravity_model = "GR"
    run.dgp_rcH0overc = 1.0
  else
    local rcH0 = model:match("^nDGP_N(%d+%.?%d*)$")
    if not rcH0 then
      error("unknown model '" .. model .. "' (expected LCDM or nDGP_N<rcH0>)")
    end
    run.gravity_model = "DGP"
    run.dgp_rcH0overc = tonumber(rcH0)
  end

  local seed, variant = realization:match("^s(%d+)(.*)$")
  if not seed or not (variant == "" or variant == "-fixed" or variant == "-fixed-inv") then
    error("unknown realization '" .. realization .. "' (expected s<seed>[-fixed[-inv]])")
  end
  run.seed = math.tointeger(tonumber(seed))
  run.fix_amplitude = variant ~= ""
  run.reverse_phases = variant == "-fixed-inv"

  run.simulation_name = target .. "_" .. model .. "_" .. realization
  run.output_folder = hicola_run_root .. "/" .. target .. "/" .. model .. "/" .. realization .. "/"
  return run
end
