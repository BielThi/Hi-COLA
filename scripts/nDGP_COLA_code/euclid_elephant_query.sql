-- CosmoHub, table nonstandardeuclidsims: Elephant halo catalogues (Rockstar) for Step 3.
-- snapshot_id 2 and 5 are the two redshifts used in Stage 1a; Box0 / Box-1 are the paired boxes.
-- No parent id in the table: subhalos cannot be removed with a column (see the doc).
SELECT simulation, snapshot_id, id, m200c, m200c_all, m200b, r200c, r_vir, rs, rs_klypin,
       vmax, vrms, np, x, y, z, vx, vy, vz
FROM nonstandardeuclidsims
WHERE simulation IN ('Elephant|LCDM|Box0', 'Elephant|LCDM|Box-1',
                     'Elephant|N1.2|Box0', 'Elephant|N1.2|Box-1')
  AND snapshot_id IN (2, 5)
  AND m200c > 1e12
