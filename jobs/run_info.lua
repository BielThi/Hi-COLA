-- Print the identity of a Hi-COLA run file as shell variable assignments.
-- usage: lua jobs/run_info.lua <run.lua>
--
--   RUN_LABEL      <target>/<model>/<realization>, or standalone/<file name>
--   TARGET MODEL REALIZATION   (TARGET=standalone for Backend_parameter_files/standalone/*.lua)
--   SIM_NAME OUTPUT_FOLDER BOXSIZE REDSHIFTS
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

local info = {}
if target and model and realization then
  info.TARGET, info.MODEL, info.REALIZATION = target, model, realization
  info.RUN_LABEL = target .. "/" .. model .. "/" .. realization
else
  local stem = run_file:match("([^/]*)%.lua$")
  info.TARGET, info.MODEL, info.REALIZATION = "standalone", "", stem
  info.RUN_LABEL = "standalone/" .. stem
end
info.SIM_NAME = simulation_name
info.OUTPUT_FOLDER = output_folder
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

for _, key in ipairs({ "RUN_LABEL", "TARGET", "MODEL", "REALIZATION", "SIM_NAME", "OUTPUT_FOLDER", "BOXSIZE",
                       "REDSHIFTS", "NPART", "NSTEPS", "Z_INIT", "GRAVITY", "COSMOLOGY", "IC_FILE", "IC_SEED",
                       "IC_FIXED", "IC_INVERTED" }) do
  print(key .. "=" .. quote(info[key]))
end
