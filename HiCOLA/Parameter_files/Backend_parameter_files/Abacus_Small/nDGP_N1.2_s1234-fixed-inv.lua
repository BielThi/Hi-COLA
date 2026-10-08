-- Hi-COLA run: Abacus_Small / nDGP_N1.2 / s1234-fixed-inv
-- All settings come from base.lua; put run-specific overrides below the dofile line.
model = "nDGP_N1.2"
realization = "s1234-fixed-inv"
dofile((debug.getinfo(1, "S").source:sub(2):match("^(.*)/") or ".") .. "/base.lua")
