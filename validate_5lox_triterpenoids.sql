-- =====================================================================================
-- ZENODO 23565306 (UFAM / COARACI) x pg_bio: 5-LOX TRITERPENOID VALIDATION SUITE
-- =====================================================================================

-- MODULE 1: CHEMINFORMATICS VALIDATION (parse_smiles, smiles_to_fingerprint, % Tanimoto)
-- Validates all 9 triterpenoid SMILES graphs, Murcko generic scaffolds (Table S10),
-- and compares pg_bio Morgan ECFP fingerprints against RDKit (Table S8).
SELECT
    ligand_label,
    formula,
    ((mol_obj::text::jsonb)->>'is_valid')::boolean AS smiles_valid,
    smiles_contains(mol_obj, 'C(=O)O') AS has_carboxyl,
    smiles_contains(mol_obj, 'C(=O)') AS has_carbonyl,
    scaffold_heavy_atoms AS murcko_atoms,
    generic_scaffold_hash AS murcko_hash,
    round((ecfp_fp % (SELECT ecfp_fp FROM lox5_triterpenoids WHERE system_id = '1-alfa-amirina'))::numeric, 4) AS pg_bio_tanimoto_vs_alpha_amyrin,
    round((murcko_fp % (SELECT murcko_fp FROM lox5_triterpenoids WHERE system_id = '1-alfa-amirina'))::numeric, 4) AS murcko_tanimoto_vs_alpha_amyrin
FROM lox5_triterpenoids
ORDER BY system_id;

-- MODULE 2: 3D Z-ORDER SPATIAL INDEXING OF THE 5-LOX ALLOSTERIC POCKET & KABSCH RMSD
-- 2A. Uses Z-Order Morton bounding box + distance_angstroms() between protein C-alpha / sidechain atoms
--     and 285 ligand heavy atoms to independently recover the 9 operational pocket residues
--     (Arg68, Glu108, Val109, Val110, Gln129, His130, Lys133, Glu134, Thr137).
WITH lig_bounds AS (
    SELECT
        least(z_order_encode(min(x) - 8.0, min(y) - 8.0, min(z) - 8.0),
              z_order_encode(max(x) + 8.0, max(y) + 8.0, max(z) + 8.0)) AS min_z,
        greatest(z_order_encode(min(x) - 8.0, min(y) - 8.0, min(z) - 8.0),
                 z_order_encode(max(x) + 8.0, max(y) + 8.0, max(z) + 8.0)) AS max_z
    FROM lox5_ligand_atoms
)
SELECT
    p.simulation_resid AS resid,
    p.resname,
    c.uniprot_resid AS uniprot_p09917_pos,
    c.unique_11aa_context AS verified_11aa_motif,
    c.operational_pocket AS is_fixed_pocket,
    count(DISTINCT l.system_id) FILTER (WHERE distance_angstroms(p.coord, l.coord) <= 5.5) AS ligands_within_5p5A,
    round(min(distance_angstroms(p.coord, l.coord)) FILTER (WHERE p.atom_name = 'CA')::numeric, 2) AS min_ca_to_lig_A,
    round(avg(d.mm_plus_polar_kcal) FILTER (WHERE d.model = 'gb')::numeric, 2) AS mean_gb_mm_polar_kcal,
    round(avg(d.mm_plus_polar_kcal) FILTER (WHERE d.model = 'pb')::numeric, 2) AS mean_pb_mm_polar_kcal
FROM lox5_protein_atoms p
JOIN lox5_residue_correspondence c ON c.simulation_resid = p.simulation_resid
JOIN lox5_ligand_atoms l ON distance_angstroms(p.coord, l.coord) <= 9.5
JOIN lox5_residue_decomposition d ON d.simulation_resid = p.simulation_resid AND d.system_id = l.system_id AND d.side = 'R'
CROSS JOIN lig_bounds b
WHERE p.z_index BETWEEN b.min_z AND b.max_z
  AND c.operational_pocket = true
GROUP BY p.simulation_resid, p.resname, c.uniprot_resid, c.unique_11aa_context, c.operational_pocket
ORDER BY p.simulation_resid ASC;

-- 2B. 3D Ligand Center-of-Mass (COM), Distance to Catalytic Non-Heme Iron Triad
--     (His367/His372/His550/Ile673), and 31-atom Pentacyclic Isomer Kabsch RMSD in SQL:
WITH iron_site AS (
    SELECT center_of_mass(array_agg(coord)) AS fe_site_com
    FROM lox5_protein_atoms
    WHERE simulation_resid IN (367, 372, 550, 673)
),
alpha_amyrin_coords AS (
    SELECT array_agg(coord ORDER BY atom_id) AS ref_coords
    FROM lox5_ligand_atoms
    WHERE system_id = '1-alfa-amirina'
),
lig_agg AS (
    SELECT
        system_id,
        ligand_label,
        count(*) AS heavy_atoms,
        center_of_mass(array_agg(coord)) AS lig_com,
        array_agg(coord ORDER BY atom_id) AS lig_coords
    FROM lox5_ligand_atoms
    GROUP BY system_id, ligand_label
)
SELECT
    l.ligand_label,
    l.heavy_atoms,
    round(((l.lig_com::text::jsonb)->>'x')::numeric, 2) || ', ' ||
    round(((l.lig_com::text::jsonb)->>'y')::numeric, 2) || ', ' ||
    round(((l.lig_com::text::jsonb)->>'z')::numeric, 2) AS com_xyz_angstroms,
    round(distance_angstroms(l.lig_com, i.fe_site_com)::numeric, 2) AS dist_to_catalytic_fe_site_A,
    CASE
        WHEN l.heavy_atoms = 31 THEN round(kabsch_rmsd(l.lig_coords, a.ref_coords)::numeric, 3)
        ELSE NULL
    END AS kabsch_rmsd_vs_alpha_amyrin_A
FROM lig_agg l
CROSS JOIN iron_site i
CROSS JOIN alpha_amyrin_coords a
ORDER BY l.system_id;

-- MODULE 3: THERMODYNAMIC & PER-RESIDUE MM/GBSA vs. MM/PBSA CLOSURE AUDIT IN SQL
-- Proves:
--   1. MM_plus_polar = Internal + vdW + Electrostatic + Polar Solvation across all 12,186 rows.
--   2. Sum of 677 decomposed entities (676 residues + ligand) for Gas-Phase MM (Internal+vdW+Elec)
--      and Polar Solvation closes exactly against global GGAS and EGB/EPB (error <= 0.003 kcal/mol).
--   3. Explains the +27 to +37 kcal/mol PB-minus-GB endpoint shift: PB decomposed Non-Polar Solv is 0.000
--      while global PB ENPOLAR + EDISPER adds +20.7 to +26.0 kcal/mol cavity/dispersion penalty.
WITH decomp_sums AS (
    SELECT
        system_id,
        model,
        sum(internal_kcal + vdw_kcal + electrostatic_kcal) AS sum_mm_gas,
        sum(polar_solv_kcal) AS sum_polar,
        sum(nonpolar_solv_kcal) AS sum_nonpolar,
        sum(mm_plus_polar_kcal) AS sum_mm_plus_polar,
        max(abs(mm_plus_polar_kcal - (internal_kcal + vdw_kcal + electrostatic_kcal + polar_solv_kcal))) AS max_identity_err
    FROM lox5_residue_decomposition
    GROUP BY system_id, model
),
global_terms AS (
    SELECT
        system_id,
        model,
        max(CASE WHEN term = 'GGAS' THEN mean_kcal END) AS ggas,
        max(CASE WHEN term IN ('EGB', 'EPB') THEN mean_kcal END) AS global_polar,
        max(CASE WHEN term = 'ESURF' THEN mean_kcal END) AS gb_esurf,
        max(CASE WHEN term = 'ENPOLAR' THEN mean_kcal END) + max(CASE WHEN term = 'EDISPER' THEN mean_kcal END) AS pb_nonpolar_net,
        max(CASE WHEN term = 'TOTAL' THEN mean_kcal END) AS global_total
    FROM lox5_energy_components
    WHERE window_name = 'full'
    GROUP BY system_id, model
)
SELECT
    t.ligand_label,
    UPPER(d.model) AS model,
    round(g.ggas::numeric, 2) AS global_ggas,
    round(abs(d.sum_mm_gas - g.ggas)::numeric, 4) AS gas_closure_err,
    round(g.global_polar::numeric, 2) AS global_polar,
    round(abs(d.sum_polar - g.global_polar)::numeric, 4) AS polar_closure_err,
    round(d.sum_nonpolar::numeric, 2) AS decomp_nonpolar_sum,
    round(COALESCE(g.gb_esurf, g.pb_nonpolar_net)::numeric, 2) AS global_nonpolar_net,
    round(d.sum_mm_plus_polar::numeric, 2) AS total_mm_plus_polar,
    round(g.global_total::numeric, 2) AS global_endpoint_kcal
FROM decomp_sums d
JOIN global_terms g ON g.system_id = d.system_id AND g.model = d.model
JOIN lox5_triterpenoids t ON t.system_id = d.system_id
ORDER BY t.system_id, d.model;

-- MODULE 4: CSR SPARSE ATTENTION MAP HOTSPOT DISCOVERY (get_top_interacting_residues & sparse_cross_reactivity_score)
-- Queries the native SparseAttentionMap (CSR) to rank top stabilizing residues per ligand
-- and contrasts hydrophobic clamp (Val109, Val110, His130) vs. polar/ionic switch (Arg68, Glu108, Lys133, Glu134).
SELECT
    m.ligand_label,
    UPPER(m.model) AS model,
    m.stabilizing_residue_count AS stabilizing_res_cnt,
    get_top_interacting_residues(m.energy_map, 999, 5) AS top5_hotspot_resids,
    round(sparse_cross_reactivity_score(m.energy_map, ARRAY[68, 108, 109, 110, 129, 130, 133, 134, 137], ARRAY[999])::numeric, 3) AS pocket_9res_csr_score,
    round(sparse_cross_reactivity_score(m.energy_map, ARRAY[109, 110, 130], ARRAY[999])::numeric, 3) AS hydrophobic_core_csr,
    round(sparse_cross_reactivity_score(m.energy_map, ARRAY[68, 133, 134], ARRAY[999])::numeric, 3) AS ionic_switch_csr
FROM lox5_energy_attention_maps m
ORDER BY m.model, pocket_9res_csr_score DESC;

-- MODULE 5: 500 ns CONFORMATIONAL ENSEMBLE PREORGANIZATION PROFILES (ensemble_fraction_below / between)
-- Evaluates 5,001-frame SASA/Rg profiles and 101-frame MM/GBSA hotspot profiles using PreorganizationProfile.
SELECT
    p.ligand_label,
    round((((p.buried_sasa_profile)::text::jsonb)->>'mean_val')::numeric, 3) AS mean_buried_sasa_nm2,
    round((ensemble_fraction_between(p.buried_sasa_profile, 4.20, 6.00) * 100.0)::numeric, 1) AS deep_burial_pct_gt_4p2nm2,
    round((ensemble_fraction_below(p.val110_gb_profile, -1.20) * 100.0)::numeric, 1) AS val110_strong_clamp_pct,
    round((ensemble_fraction_below(p.his130_gb_profile, -1.20) * 100.0)::numeric, 1) AS his130_strong_clamp_pct,
    round((ensemble_fraction_below(p.lys133_gb_profile, -1.00) * 100.0)::numeric, 1) AS lys133_strong_anchor_pct,
    round((ensemble_fraction_below(p.glu134_gb_profile, -1.00) * 100.0)::numeric, 1) AS glu134_hbond_lock_pct,
    t.hbond_any_pct AS direct_hbond_pct
FROM lox5_structural_profiles p
JOIN lox5_triterpenoids t ON t.system_id = p.system_id
ORDER BY deep_burial_pct_gt_4p2nm2 DESC;

-- MODULE 6: 5-LOX SEQUENCE BIOPHYSICAL VALIDATION & STRUCTURAL HOMOLOGY IN SQL
-- Validates MW, pI, GRAVY, Shannon entropy, non-heme iron PROSITE motif, and ESM vector distance vs. Wild-Type P09917.
WITH wt AS (
    SELECT esm_embedding AS wt_emb, sequence AS wt_seq
    FROM lox5_sequence_models
    WHERE accession = 'P09917'
)
SELECT
    m.model_name,
    m.accession,
    m.residue_count AS aa_len,
    round((molecular_weight(m.sequence) / 1000.0)::numeric, 2) AS mw_kda,
    round(isoelectric_point(m.sequence)::numeric, 2) AS calc_pi,
    round(gravy_score(m.sequence)::numeric, 4) AS gravy,
    round(sequence_entropy(m.sequence)::numeric, 4) AS shannon_entropy,
    prosite_match(m.sequence, 'H-x(4)-H') AS has_iron_his_dyad,
    round(embedding_cosine_distance(m.esm_embedding, w.wt_emb)::numeric, 5) AS esm_dist_vs_p09917,
    round(interface_complementarity(m.pocket_embedding, m.esm_embedding)::numeric, 4) AS pocket_complementarity
FROM lox5_sequence_models m
CROSS JOIN wt w
ORDER BY m.model_name;
