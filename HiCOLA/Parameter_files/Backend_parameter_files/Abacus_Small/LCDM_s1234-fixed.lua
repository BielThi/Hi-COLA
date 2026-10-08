-- Hi-COLA run: Abacus_Small / LCDM / s1234-fixed
-- All settings come from base.lua; put run-specific overrides below the dofile line.
model = "LCDM"
realization = "s1234-fixed"
dofile((debug.getinfo(1, "S").source:sub(2):match("^(.*)/") or ".") .. "/base.lua")
