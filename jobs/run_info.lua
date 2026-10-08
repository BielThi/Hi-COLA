-- Print the identity of a Hi-COLA run file as shell variable assignments.
-- usage: lua jobs/run_info.lua <run.lua>
--
--   RUN_LABEL      <target>/<model>/<realization>, or standalone/<file name>
--   TARGET MODEL REALIZATION   (TARGET=standalone for Backend_parameter_files/standalone/*.lua)
--   SIM_NAME OUTPUT_FOLDER BOXSIZE REDSHIFTS

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

for _, key in ipairs({ "RUN_LABEL", "TARGET", "MODEL", "REALIZATION", "SIM_NAME", "OUTPUT_FOLDER", "BOXSIZE", "REDSHIFTS" }) do
  print(key .. "=" .. quote(info[key]))
end
