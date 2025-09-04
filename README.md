# PalmFFT

The Infinite Monkey Derivative of [PocketFFT](https://gitlab.mpcdf.mpg.de/mtr/pocketfft).

It's technically a dirty-room implementation, but my brain is clean. fuck.

I wanted to release this repository as Public Domain but released as BSD 3-clause cuz am afraid of lawyers kek

## What it can do

- Complex FFT, Basic complex arithmetics. That's it.
- Real FFT? No, I have my life.
- Try searching unsafe, You won't find any
- Peak Dependency Minimalism™

## How to use

```sh
cargo add palmfft --feature num-complex
```

```rust
use palmfft::{CfftPlan, Complex};

// just a random data
let n = 1024;
let mut data = vec![Complex::new(0.0, 0.0); n];
// ... fuck around with the data vector ...
let fct = 1.0 / data.len() as f64; // this is a scale factor

// initialise fft plan, no mismatch allowed
let plan = CfftPlan::new(data.len());

plan.forward(&mut data, 1.0); // fft
plan.backward(&mut data, fct); // ifft
```

Again, fuck.