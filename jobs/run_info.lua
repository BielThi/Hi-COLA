-- Print the identity of a Hi-COLA run file as shell variable assignments.
-- usage: lua jobs/run_info.lua <run.lua>
--
--   RUN_LABEL      <target>/<model>/<realization>, or standalone/<file name>
--   TARGET MODEL REALIZATION   from where the file is: <target>/<model>_<realization>.lua,
--                              or standalone/<name>.lua (TARGET=standalone)
--   SIM_NAME OUTPUT_FOLDER RUN_ROOT BOXSIZE REDSHIFTS
--                  RUN_ROOT = output_folder without /<target>/<model>/<realization>/
--
-- Fails if simulation_name or output_folder do not follow that naming: the halo finder and
-- scripts/collect_products.py find a run's outputs by these names.
--   NPART NSTEPS Z_INIT GRAVITY COSMOLOGY IC_FILE IC_SEED IC_FIXED IC_INVERTED   (for the job logs)

local run_file = arg[1]
assert(run_file, "usage: lua run_info.lua <run.lua>")

local real_print = print
print = function() end
dofile(run_file)
print = real_print

local function quote(value)
  return "'" .. (tostring(value):gsub("'", "'\\''")) .. "'"
end

local redshifts = {}
for i, z in ipairs(output_redshifts or {}) do redshifts[i] = string.format("%.3f", z) end

local function fail(message)
  io.stderr:write(run_file .. ": " .. message .. "\n")
  os.exit(1)
end

local info = {}
local folder, stem = run_file:match("([^/]+)/([^/]+)%.lua$")
if not folder then fail("expected a path like <target>/<model>_<realization>.lua") end
local suffix
if folder == "standalone" then
  info.TARGET, info.MODEL, info.REALIZATION = "standalone", "", stem
  info.RUN_LABEL = "standalone/" .. stem
  suffix = "/standalone/" .. stem .. "/"
else
  local model, realization = stem:match("^(.+)_(s%d+[%w%-]*)$")
  if not model then fail("file name should be <model>_<realization>.lua, e.g. nDGP_N1.2_s1234-fixed.lua") end
  info.TARGET, info.MODEL, info.REALIZATION = folder, model, realization
  info.RUN_LABEL = folder .. "/" .. model .. "/" .. realization
  local expected = folder .. "_" .. model .. "_" .. realization
  if simulation_name ~= expected then
    fail('simulation_name is "' .. tostring(simulation_name) .. '", expected "' .. expected .. '"')
  end
  suffix = "/" .. info.RUN_LABEL .. "/"
end
if type(output_folder) ~= "string" or output_folder:sub(-#suffix) ~= suffix then
  fail('output_folder "' .. tostring(output_folder) .. '" should end with "' .. suffix .. '"')
end
info.SIM_NAME = simulation_name
info.OUTPUT_FOLDER = output_folder
info.RUN_ROOT = output_folder:sub(1, -#suffix - 1)
info.BOXSIZE = simulation_boxsize
info.REDSHIFTS = table.concat(redshifts, " ")

local steps = {}
for i, n in ipairs(timestep_nsteps or {}) do steps[i] = tostring(n) end
info.NPART = particle_Npart_1D
info.NSTEPS = table.concat(steps, " + ")
info.Z_INIT = ic_initial_redshift
info.GRAVITY = gravity_model
if gravity_model == "DGP" then
  info.GRAVITY = "DGP (rc*H0/c = " .. tostring(gravity_model_dgp_rcH0overc) .. ")"
end
info.COSMOLOGY = string.format("%s, h = %s, Omega_cdm = %s, Omega_b = %s, A_s = %s, n_s = %s",
  tostring(cosmology_model), tostring(cosmology_h), tostring(cosmology_OmegaCDM),
  tostring(cosmology_Omegab), tostring(cosmology_As), tostring(cosmology_ns))
info.IC_FILE = ic_input_filename
info.IC_SEED = ic_random_seed
info.IC_FIXED = ic_fix_amplitude
info.IC_INVERTED = ic_reverse_phases

for _, key in ipairs({ "RUN_LABEL", "TARGET", "MODEL", "REALIZATION", "SIM_NAME", "OUTPUT_FOLDER", "RUN_ROOT", "BOXSIZE",
                       "REDSHIFTS", "NPART", "NSTEPS", "Z_INIT", "GRAVITY", "COSMOLOGY", "IC_FILE", "IC_SEED",
                       "IC_FIXED", "IC_INVERTED" }) do
  print(key .. "=" .. quote(info[key]))
end
