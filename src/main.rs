mod sim;

use sim::*;

use polars::prelude::*;
use rand::prelude::*;
use rayon::prelude::*;
use std::fs::File;
use std::iter;

fn main() {
    // Container for randomized simulation settings
    let mut settings: Vec<SimulationSettings> = Vec::new();

    // Uniform distribution for randomization
    let unif = rand::distr::Uniform::new(0.0, 1.0).unwrap();

    let mut rng: ThreadRng = rand::rng();
    for _ in 0..32 {
        let mut s: SimulationSettings = SimulationSettings::default();

        // Bounds for the parameter space
        let ub = 10.0;
        let lb = 0.1;
        
        // Scale reverse rates down to prevent implausible scenarios(polymerases going backwards on average)
        let s_lb = 1.0;
        let s_ub = 0.01;
        
        // Randomize settings ----
        // TF2D
        s.tf2d_conc = unif.sample(&mut rng) * (ub - lb) + lb;
        s.tf2d_kon = unif.sample(&mut rng) * (ub - lb) + lb;
        s.tf2d_koff = s.tf2d_kon * (unif.sample(&mut rng) * (s_ub - s_lb) + s_lb);

        // TF2A
        s.tf2a_conc = unif.sample(&mut rng) * (ub - lb) + lb;
        s.tf2a_kon = unif.sample(&mut rng) * (ub - lb) + lb;
        s.tf2a_koff = s.tf2a_kon * (unif.sample(&mut rng) * (s_ub - s_lb) + s_lb);

        // TF2B
        s.tf2b_conc = unif.sample(&mut rng) * (ub - lb) + lb;
        s.tf2b_kon = unif.sample(&mut rng) * (ub - lb) + lb;
        s.tf2b_koff = s.tf2b_kon * (unif.sample(&mut rng) * (s_ub - s_lb) + s_lb);

        // RNAP-TF2F
        s.rnap_tf2f_conc = unif.sample(&mut rng) * (ub - lb) + lb;
        s.rnap_tf2f_kon = unif.sample(&mut rng) * (ub - lb) + lb;
        s.rnap_tf2f_koff = s.rnap_tf2f_kon * (unif.sample(&mut rng) * (s_ub - s_lb) + s_lb);

        // TF2E
        s.tf2e_conc = unif.sample(&mut rng) * (ub - lb) + lb;
        s.tf2e_kon = unif.sample(&mut rng) * (ub - lb) + lb;
        s.tf2e_koff = s.tf2e_kon * (unif.sample(&mut rng) * (s_ub - s_lb) + s_lb);

        // TF2H
        s.tf2h_conc = unif.sample(&mut rng) * (ub - lb) + lb;
        s.tf2h_kon = unif.sample(&mut rng) * (ub - lb) + lb;
        s.tf2h_koff = s.tf2h_kon * (unif.sample(&mut rng) * (s_ub - s_lb) + s_lb);

        s.tf2h_phospho = unif.sample(&mut rng) * (ub - lb) + lb;
        s.exchange = unif.sample(&mut rng) * (ub - lb) + lb;

        // Transcription Elongation Parameters
        // Translocation Step (1)
        s.translocate_kfor = unif.sample(&mut rng) * (ub - lb) + lb;
        s.translocate_krev = s.translocate_kfor * (unif.sample(&mut rng) * (s_ub - s_lb) + s_lb);

        // NTP Binding Step (2)
        s.binding_kfor = unif.sample(&mut rng) * (ub - lb) + lb;
        s.binding_krev = s.binding_kfor * (unif.sample(&mut rng) * (s_ub - s_lb) + s_lb);

        // NTP Alignment Step (3)
        s.alignment_kfor = unif.sample(&mut rng) * (ub - lb) + lb;
        s.alignment_krev = s.alignment_kfor * (unif.sample(&mut rng) * (s_ub - s_lb) + s_lb);

        // Catalysis Step (4)
        s.catalysis_kfor = unif.sample(&mut rng) * (ub - lb) + lb;
        s.catalysis_krev = s.catalysis_kfor * (unif.sample(&mut rng) * (s_ub - s_lb) + s_lb);

        // PPi Release Step (5)
        s.ppi_release = unif.sample(&mut rng) * (ub - lb) + lb;

        // Paused State
        // Pausing Step (-1)
        s.pausing_krev = unif.sample(&mut rng) * (ub - lb) + lb;
        s.pausing_kfor = s.pausing_krev * (unif.sample(&mut rng) * (s_ub - s_lb) + s_lb);

        // Half-Backtrack Step (-2)
        s.p_halftranslocate_krev = unif.sample(&mut rng) * (ub - lb) + lb;
        s.p_halftranslocate_kfor =
            s.p_halftranslocate_krev * (unif.sample(&mut rng) * (s_ub - s_lb) + s_lb);

        // Fraying Step (-3)
        s.p_fray_krev = unif.sample(&mut rng) * (ub - lb) + lb;
        s.p_fray_kfor = s.p_fray_krev * (unif.sample(&mut rng) * (s_ub - s_lb) + s_lb);

        // Full Backtrack Step (-4)
        s.p_backtrack_krev = unif.sample(&mut rng) * (ub - lb) + lb;
        s.p_backtrack_kfor = s.p_backtrack_krev * (unif.sample(&mut rng) * (s_ub - s_lb) + s_lb);

        // Elongation Regulating Parameters ----
        s.ntp_conc = unif.sample(&mut rng) * (ub - lb) + lb;
        s.topo_econc = unif.sample(&mut rng) * (ub - lb) + lb;

        // RNA Death Rate
        s.death_rate = 0.01;

        // Append to settings
        settings.push(s);
    }

    // Evaluate all settings in parallel
    let data: Vec<(Vec<u32>, Vec<u32>)> = settings
        .par_iter()
        .enumerate()
        .map(|(i, s)| {
            // Construct simulation
            let mut sim: Simulation = Simulation::default();

            // Update settings
            sim.settings = *s;

            // Run simulation and sample stationary distribution
            let sample: Vec<u32> = sim.run();
            let sim_id: Vec<u32> = iter::repeat(i as u32).take(sample.len()).collect();
            println!("{i}");

            (sample, sim_id)
        })
        .collect();

    // Concatenate vectors
    let mut samples: Vec<u32> = Vec::new();
    let mut sim_ids: Vec<u32> = Vec::new();
    for (mut sample, mut sim_id) in data {
        samples.append(&mut sample);
        sim_ids.append(&mut sim_id);
    }

    // Format into dataframe
    let mut data: DataFrame = df!(
        "count" => samples,
        "sim_id" => sim_ids
    )
    .unwrap();

    let mut file: File = File::create("data.csv").unwrap();

    CsvWriter::new(&mut file)
        .include_header(true)
        .with_separator(b',')
        .finish(&mut data)
        .unwrap();
}
