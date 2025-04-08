use rand::prelude::*;

// Structures ----
#[derive(Clone, Copy)]
pub struct Polymerase {
    phase: i8,
    supercoil_us: f32,
    supercoil_ds: f32,
    pos: i32,
}

#[derive(Clone, Copy)]
pub struct SimulationSettings {
    // Misc. Variables
    pub gene_length: i32, // Length of the gene
    pub pol_size: i32, // Size of the polymerase
    pub dt: f32, // Stepsize
    pub update: f32, // Time between samples from the stationary distribution
    pub n: u32, // Number of samples from the stationary distribution

    // Transcription Initiation Parameters
    // TF2D
    pub tf2d_conc: f32,
    pub tf2d_kon: f32,
    pub tf2d_koff: f32,

    // TF2A
    pub tf2a_conc: f32,
    pub tf2a_kon: f32,
    pub tf2a_koff: f32,

    // TF2B
    pub tf2b_conc: f32,
    pub tf2b_kon: f32,
    pub tf2b_koff: f32,

    // RNAP-TF2F
    pub rnap_tf2f_conc: f32,
    pub rnap_tf2f_kon: f32,
    pub rnap_tf2f_koff: f32,

    // TF2E
    pub tf2e_conc: f32,
    pub tf2e_kon: f32,
    pub tf2e_koff: f32,

    // TF2H
    pub tf2h_conc: f32,
    pub tf2h_kon: f32,
    pub tf2h_koff: f32,

    pub tf2h_phospho: f32, // TF2H Phosphorylation Rate
    pub exchange: f32,     // GTF - GEF Exchange Rate

    // Transcription Elongation Parameters
    // Translocation Step (1)
    pub translocate_kfor: f32,
    pub translocate_krev: f32,

    // NTP Binding Step (2)
    pub binding_kfor: f32,
    pub binding_krev: f32,

    // NTP Alignment Step (3)
    pub alignment_kfor: f32,
    pub alignment_krev: f32,

    // Catalysis Step (4)
    pub catalysis_kfor: f32,
    pub catalysis_krev: f32,

    // PPi Release Step (5)
    pub ppi_release: f32,

    // Pausing Step (-1)
    pub pausing_kfor: f32,
    pub pausing_krev: f32,

    // Half-Backtrack Step (-2)
    pub p_halftranslocate_kfor: f32,
    pub p_halftranslocate_krev: f32,

    // Fraying Step (-3)
    pub p_fray_kfor: f32,
    pub p_fray_krev: f32,

    // Full Backtrack Step (-4)
    pub p_backtrack_kfor: f32,
    pub p_backtrack_krev: f32,

    // Elongation Regulating Parameters
    pub ntp_conc: f32,
    pub topo_econc: f32,

    // RNA Death Rate
    pub death_rate: f32,
}

pub struct Simulation {
    pub settings: SimulationSettings,
    pub time: f32,
    pub tick: u64,
    pub rng: ThreadRng,
    pub polymerases: Vec<Polymerase>,
    pub promoter_state: u8,
    pub rna_count: u32,
}

// Methods ----
impl Default for Polymerase {
    fn default() -> Self {
        Polymerase {
            phase: 0,
            supercoil_us: 0.0,
            supercoil_ds: 0.0,
            pos: 1,
        }
    }
}

impl Default for SimulationSettings {
    fn default() -> Self {
        SimulationSettings {
            // Misc. Variables
            gene_length: 500,
            pol_size: 15,
            dt: 0.01,
            update: 1e5,
            n: 100000,

            // Transcription Initiation Parameters
            // TF2D
            tf2d_conc: 1.0,
            tf2d_kon: 0.01,
            tf2d_koff: 0.01,

            // TF2A
            tf2a_conc: 1.0,
            tf2a_kon: 10.0,
            tf2a_koff: 0.01,

            // TF2B
            tf2b_conc: 1.0,
            tf2b_kon: 0.07,
            tf2b_koff: 0.01,

            // RNAP-TF2F
            rnap_tf2f_conc: 1.0,
            rnap_tf2f_kon: 0.07,
            rnap_tf2f_koff: 0.01,

            // TF2E
            tf2e_conc: 1.0,
            tf2e_kon: 0.07,
            tf2e_koff: 0.01,

            // TF2H
            tf2h_conc: 1.0,
            tf2h_kon: 0.01,
            tf2h_koff: 0.01,

            tf2h_phospho: 1.0, // TF2H Phosphorylation Rate
            exchange: 1.0,     // GTF - GEF Exchange Rate

            // Transcription Elongation Parameters
            // Translocation Step (1)
            translocate_kfor: 10.0,
            translocate_krev: 0.01,

            // NTP Binding Step (2)
            binding_kfor: 10.0,
            binding_krev: 0.01,

            // NTP Alignment Step (3)
            alignment_kfor: 10.0,
            alignment_krev: 0.01,

            // Catalysis Step (4)
            catalysis_kfor: 10.0,
            catalysis_krev: 0.01,

            // PPi Release Step (5)
            ppi_release: 10.0,

            // Pausing Step (-1)
            pausing_kfor: 1.00,
            pausing_krev: 0.2, // fixed

            // Half-Backtrack Step (-2)
            p_halftranslocate_kfor: 0.05,
            p_halftranslocate_krev: 0.1,

            // Fraying Step (-3)
            p_fray_kfor: 0.02,
            p_fray_krev: 0.1,

            // Full Backtrack Step (-4)
            p_backtrack_kfor: 0.01,
            p_backtrack_krev: 0.01,

            // Elongation Regulating Parameters
            ntp_conc: 10.0,
            topo_econc: 10.0,

            // RNA Death Rate
            death_rate: 0.0001,
        }
    }
}

impl Default for Simulation {
    fn default() -> Self {
        Simulation {
            settings: SimulationSettings::default(),
            rng: rand::rng(),
            time: 0.0,
            tick: 0,
            polymerases: Vec::<Polymerase>::new(),
            promoter_state: 0,
            rna_count: 0,
        }
    }
}

impl Simulation {
    pub fn run(&mut self) -> Vec<u32> {
        // Number of ticks for each update to minimize autocorrelation
        let update: u64 = (self.settings.update / self.settings.dt) as u64;

        // Sample for stationary distribution
        let mut sample: Vec<u32> = Vec::new();

        // Run simulation
        let start = 1e5; // When to start sampling
        while self.time < start + self.settings.update * self.settings.n as f32 {
            // Iterate
            self.step();

            if start < self.time && self.tick % update == 0 {
                sample.push(self.rna_count);
            }
        }

        sample
    }

    pub fn step(&mut self) {
        // Simulate transcription initiation
        if self.polymerases.is_empty() {
            self.initiation();
        } else {
            if self
                .polymerases
                .iter()
                .all(|pol| pol.pos > self.settings.pol_size)
            {
                self.initiation();
            }
        }

        if !self.polymerases.is_empty() {
            // Simulate elongation
            self.elongation();
        }

        if !self.polymerases.is_empty() {
            // Cancel out torsional strain between polymerases
            self.cancel_strain();

            // Simulate the action of topoisomerase
            self.topoisomerase();
        }

        // Simulate RNA death
        self.death_process();

        // Update time
        self.tick += 1;
        self.time = (self.tick as f32) * self.settings.dt;
    }

    // Simulate transcription initiation
    pub fn initiation(&mut self) {
        let s: &SimulationSettings = &self.settings;

        match self.promoter_state {
            0 => {
                // Unbound promoter
                // forward | TF2D binds
                if self
                    .rng
                    .random_bool((s.tf2d_conc * s.tf2d_kon * s.dt) as f64)
                {
                    self.promoter_state = 1;
                }
            }
            1 => {
                // TF2D bound promoter
                if self
                    .rng
                    .random_bool((s.tf2a_conc * s.tf2a_kon * s.dt) as f64)
                {
                    self.promoter_state = 1 + 1;
                }
                if self.rng.random_bool((s.tf2d_koff * s.dt) as f64) {
                    self.promoter_state = 1 - 1;
                }
            }
            2 => {
                // TF2D, TF2A bound promoter
                // forward | TF2B binds
                if self
                    .rng
                    .random_bool((s.tf2b_conc * s.tf2b_kon * s.dt) as f64)
                {
                    self.promoter_state = 2 + 1;
                }
                // reverse | TF2A dissociates
                if self.rng.random_bool((s.tf2a_koff * s.dt) as f64) {
                    self.promoter_state = 2 - 1;
                }
            }
            3 => {
                // TF2D, TF2A, TF2B bound promoter
                // forward | RNAP-TF2F binds
                if self
                    .rng
                    .random_bool((s.rnap_tf2f_conc * s.rnap_tf2f_kon * s.dt) as f64)
                {
                    self.promoter_state = 3 + 1;
                }
                // reverse | TF2B dissociates
                if self.rng.random_bool((s.tf2b_koff * s.dt) as f64) {
                    self.promoter_state = 3 - 1;
                }
            }
            4 => {
                // TF2D, TF2A, TF2B, RNAP-TF2F bound promoter
                // forward | TF2E binds
                if self
                    .rng
                    .random_bool((s.tf2e_conc * s.tf2e_kon * s.dt) as f64)
                {
                    self.promoter_state = 4 + 1;
                }
                // reverse | RNAP-TF2F dissociates
                if self.rng.random_bool((s.rnap_tf2f_koff * s.dt) as f64) {
                    self.promoter_state = 4 - 1;
                }
            }
            5 => {
                // TF2D, TF2A, TF2B, RNAP-TF2F, TF2E bound promoter
                // forward | TF2H binds
                if self
                    .rng
                    .random_bool((s.tf2h_conc * s.tf2h_kon * s.dt) as f64)
                {
                    self.promoter_state = 5 + 1;
                }
                // reverse | TF2E dissociates
                if self.rng.random_bool((s.tf2e_koff * s.dt) as f64) {
                    self.promoter_state = 5 - 1;
                }
            }
            6 => {
                // TF2D, TF2A, TF2B, RNAP-TF2F, TF2E, TF2H bound promoter
                // forward | TF2H phosphorylates (ligand)
                if self.rng.random_bool((s.tf2h_phospho * s.dt) as f64) {
                    self.promoter_state = 6 + 1;
                }
                // reverse | TF2H dissociates
                if self.rng.random_bool((s.tf2h_koff * s.dt) as f64) {
                    self.promoter_state = 6 - 1;
                }
            }
            7 => {
                // TF2H-phosphorylated (ligand) promoter
                // forward | exchange of GTFs for ETFs (successful escape)
                if self.rng.random_bool((s.exchange * s.dt) as f64) {
                    // Everything unbinds from promoter
                    self.promoter_state = 0;

                    // Spawn new `polymerase` agent
                    self.polymerases.push(Polymerase::default());
                }
            }
            _ => {}
        }
    }

    // Transcription Elongation(polymerase-specific)
    pub fn elongation(&mut self) {
        let s: &SimulationSettings = &self.settings;
        let mut i: usize = 0;
        let mut removed: bool = false;

        while i < self.polymerases.len() {
            match self.polymerases[i].phase {
                0 => {
                    // Pre-translocated State
                    // Check if physically blocked
                    let blocked: bool = self.polymerases
                        .iter()
                        .any(|other: &Polymerase| {
                            // Ensure polymerases are distinct
                            other.pos != self.polymerases[i].pos &&

                            // Check if `other` is right in front of `pol`
                            0 <= (other.pos - self.polymerases[i].pos) && (other.pos - self.polymerases[i].pos) <= s.pol_size
                        });

                    // Attempt translocation if unblocked
                    if !blocked {
                        // Calculate effect of torsional strain
                        let strain_diff: f32 =
                            self.polymerases[i].supercoil_ds - self.polymerases[i].supercoil_us;
                        let effect: f32 = 1.0 / (1.0 + f32::exp(strain_diff - 100.0));

                        if self
                            .rng
                            .random_bool((s.translocate_kfor * effect * s.dt) as f64)
                        {
                            // Move forward along DNA strand
                            self.polymerases[i].phase = 0 + 1;
                            self.polymerases[i].pos += 1;

                            // Generate torsional strain
                            self.polymerases[i].supercoil_ds += 1.0;
                            self.polymerases[i].supercoil_us -= 1.0;

                            // Check if at termination site
                            if i == 0 {
                                if self.polymerases[i].pos >= self.settings.gene_length {
                                    // Produce RNA
                                    self.rna_count += 1;

                                    // Remove polymerase
                                    self.polymerases.remove(i);

                                    removed = true;
                                    continue;
                                }
                            }
                        }
                    }
                }
                1 => {
                    // Post-Translocated State
                    // Bind correct nucleotide
                    if self
                        .rng
                        .random_bool((s.binding_kfor * s.ntp_conc * s.dt) as f64)
                    {
                        self.polymerases[i].phase = 1 + 1;
                    }

                    // Move back to pre-translocated state
                    if self.polymerases[i].phase == 1
                        && self.rng.random_bool((s.translocate_krev * s.dt) as f64)
                    {
                        self.polymerases[i].phase = 1 - 1;
                    }

                    // Move to ePEC state(paused polymerase)
                    if self.polymerases[i].phase == 1
                        && self.rng.random_bool((s.pausing_kfor * s.dt) as f64)
                    {
                        self.polymerases[i].phase = -1;
                    }
                }
                2 => {
                    // NTP-bound State
                    // Align nucleotide
                    if self.rng.random_bool((s.alignment_kfor * s.dt) as f64) {
                        self.polymerases[i].phase = 2 + 1;
                    }

                    // Nucleotide dissociates
                    if self.polymerases[i].phase == 2
                        && self.rng.random_bool((s.binding_krev * s.dt) as f64)
                    {
                        self.polymerases[i].phase = 2 - 1;
                    }
                }
                3 => {
                    // NTP-aligned State
                    // Catalysis occurs
                    if self.rng.random_bool((s.catalysis_kfor * s.dt) as f64) {
                        self.polymerases[i].phase = 3 + 1;
                    }

                    // NTP unaligns
                    if self.polymerases[i].phase == 3
                        && self.rng.random_bool((s.alignment_krev * s.dt) as f64)
                    {
                        self.polymerases[i].phase = 3 - 1;
                    }
                }
                4 => {
                    // NTP-aligned State
                    // PPi release -> pre-translocated state
                    if self.rng.random_bool((s.ppi_release * s.dt) as f64) {
                        self.polymerases[i].phase = 0;
                    }

                    // Catalysis reverses
                    if self.polymerases[i].phase == 4
                        && self.rng.random_bool((s.catalysis_krev * s.dt) as f64)
                    {
                        self.polymerases[i].phase = 4 - 1;
                    }
                }
                // Paused Polymerase ----
                -1 => {
                    // ePEC State
                    // Backtrack to paused pre-translocated state
                    if self
                        .rng
                        .random_bool((s.p_halftranslocate_kfor * s.dt) as f64)
                    {
                        self.polymerases[i].phase = -1 - 1;
                    }

                    // Unpause to moving pre-translocated state
                    if self.polymerases[i].phase == -1
                        && self.rng.random_bool((s.pausing_krev * s.dt) as f64)
                    {
                        self.polymerases[i].phase = -1 + 1;
                    }
                }
                -2 => {
                    // Paused Pre-Translocated State
                    // Move to frayed paused state
                    if self.rng.random_bool((s.p_fray_kfor * s.dt) as f64) {
                        self.polymerases[i].phase = -2 - 1;
                    }

                    // Move back to ePEC state
                    if self.polymerases[i].phase == -2
                        && self
                            .rng
                            .random_bool((s.p_halftranslocate_krev * s.dt) as f64)
                    {
                        self.polymerases[i].phase = -2 + 1;
                    }
                }
                -3 => {
                    // Paused Frayed State
                    if self.rng.random_bool((s.p_backtrack_kfor * s.dt) as f64) {
                        self.polymerases[i].phase = -3 - 1;
                    }

                    // Move to pre-translocated state
                    if self.polymerases[i].phase == -3
                        && self.rng.random_bool((s.p_fray_krev * s.dt) as f64)
                    {
                        self.polymerases[i].phase = -3 + 1;
                    }
                }
                -4 => {
                    // Backtracked State
                    // Move to frayed state
                    if self.rng.random_bool((s.p_backtrack_krev * s.dt) as f64) {
                        self.polymerases[i].phase = -4 + 1;
                    }
                }
                _ => {}
            }

            // Increment index
            if !removed {
                i += 1;
            }
            removed = false;
        }
    }

    /// RNA death process
    pub fn death_process(&mut self) {
        let s: &SimulationSettings = &self.settings;

        if self
            .rng
            .random_bool((s.death_rate * (self.rna_count as f32) * s.dt) as f64)
        {
            self.rna_count -= 1;
        }
    }

    /// Cancel out strain between polymerases
    pub fn cancel_strain(&mut self) {
        // Iterate over pairs of polymerases (downstream, upstream)

        for i in 0..(self.polymerases.len() - 1) {
            // Split vector into distinct mutable slices
            let (f, s) = self.polymerases.split_at_mut(i + 1);

            // Identify pair of polymerases
            let pol_ds: &mut Polymerase = &mut f[i];
            let pol_us: &mut Polymerase = &mut s[0];

            // Cancel out torsional strain
            let avg: f32 = (pol_ds.supercoil_us + pol_us.supercoil_ds) / 2.0;

            pol_ds.supercoil_us = avg;
            pol_us.supercoil_ds = avg;
        }
    }

    /// Simulate if topoisomerase unwinds DNA
    pub fn topoisomerase(&mut self) {
        let s: &SimulationSettings = &self.settings;

        if self.rng.random_bool((s.topo_econc * s.dt) as f64) {
            let hit: i32 = self.rng.random_range(0..s.gene_length);

            // Find the index of the closest polymerase ahead
            if let Some((pol_index, _)) = self
                .polymerases
                .iter()
                .enumerate()
                .filter(|(_, pol)| pol.pos > hit && (pol.pos - hit).abs() < 50) // Check if nearby
                .min_by_key(|(_, pol)| pol.pos - hit)
            {
                // Remove torsional strain
                self.polymerases[pol_index].supercoil_us = 0.0;
            }

            // Find the index of the closest polymerase behind
            if let Some((pol_index, _)) = self
                .polymerases
                .iter()
                .enumerate()
                .filter(|(_, pol)| pol.pos < hit && (pol.pos - hit).abs() < 50) // Check if nearby
                .min_by_key(|(_, pol)| hit - pol.pos) 
            {
                // Remove torsional strain
                self.polymerases[pol_index].supercoil_ds = 0.0;
            }
        }
    }
}
