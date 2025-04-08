To run the script in `/src` to generate the stationary distribution samples, ensure you've installed `Rust` and then run `cargo run --release` in whatever shell you use. You can modify the miscellaneous variable defaults freely in `/src/sim.rs` as those are not randomized per-simulation. You may also change the bounds of the parameter space explored by changing `ub, lb, s_ub, s_lb`.

To generate the figures from the report, ensure you've installed `R` and `renv`, then once you've called `R` in this directory the required packages should install automatically and you may run `source('fano.R')`, `source('bimodality.R')`, or `source('animation.R')` to generate the figures used.

The `Netlogo` code is found at `/netlogo/full-transcription.nlogo`.