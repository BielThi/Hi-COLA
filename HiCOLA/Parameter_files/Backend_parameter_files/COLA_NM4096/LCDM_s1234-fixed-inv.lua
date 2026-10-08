-- Hi-COLA run: COLA_NM4096 / LCDM / s1234-fixed-inv
-- All settings come from base.lua; put run-specific overrides below the dofile line.
model = "LCDM"
realization = "s1234-fixed-inv"
dofile((debug.getinfo(1, "S").source:sub(2):match("^(.*)/") or ".") .. "/base.lua")
