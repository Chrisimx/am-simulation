extern crate rustfft;
extern crate plotters;

use std::cmp::min;
use plotters::prelude::*;
use biquad::{Biquad, Coefficients, DirectForm1, Hertz, ToHertz, Type};
use std::f64::consts::PI;
use num_complex::Complex;
use rustfft::FftPlanner;

fn generate_am_signal(freq_carrier: f64, freq_mod: f64, mod_index: f64, fs: f64, duration: f64) -> Vec<f64> {
    let samples = (duration * fs) as usize;
    let mut signal = Vec::with_capacity(samples);

    for n in 0..samples {
        let t = n as f64 / fs;
        // AM-moduliertes Signal
        let carrier = (2.0 * PI * freq_carrier * t).sin();
        let modulating = (2.0 * PI * freq_mod * t).sin();
        let am_signal = (1.0 + mod_index * modulating) * carrier;
        signal.push(am_signal);
    }

    signal
}

fn filter_signal(signal: &Vec<f64>, fs: biquad::frequency::Hertz<f64>, target_freq: Hertz<f64>, q_factor: f64) -> Vec<f64> {
    // Create a band-pass filter for the target frequency
    let coeffs = Coefficients::<f64>::from_params(Type::BandPass, fs, target_freq, q_factor).unwrap();
    let mut biquad1 = DirectForm1::<f64>::new(coeffs);

    let mut output_vec1 = Vec::new();
    for elem in signal {
        output_vec1.push(biquad1.run(elem.to_owned()))
    }
    output_vec1
}
fn extract_peaks(signal: &Vec<f64>, fs: f64) -> Vec<f64> {
    let abs_signal: Vec<_> = signal.iter().map(|&sig| sig.abs()).collect();
    let mut peaks = Vec::new();

    for i in 1..signal.len() - 1 {
        if abs_signal[i] > 0.0 && abs_signal[i] > abs_signal[i - 1] && abs_signal[i] > abs_signal[i + 1] {
            peaks.push(abs_signal[i]);
        }
    }

    peaks
}

fn remove_dc_component(signal: &Vec<f64>, fs: biquad::frequency::Hertz<f64>, cutoff_freq: Hertz<f64>) -> Vec<f64> {
    // Hochpassfilter mit einer niedrigen Grenzfrequenz
    let coeffs = Coefficients::<f64>::from_params(Type::HighPass, fs, cutoff_freq, 1.0).unwrap();
    let mut biquad1 = DirectForm1::<f64>::new(coeffs);

    let mut output_vec = Vec::new();
    for &value in signal {
        output_vec.push(biquad1.run(value));
    }

    output_vec
}
fn do_fft(filename: &str, fs: f64, signal: &Vec<f64>, upper_limit: f64, lower_freq_limit: f64) {
    let mut planner = FftPlanner::<f64>::new();
    let fft = planner.plan_fft_forward(signal.len());
    let mut spectrum: Vec<Complex<f64>> = signal.iter().map(|x| Complex::new(x.to_owned(), 0.0)).collect();
    fft.process(&mut spectrum);

    // Berechne die Amplituden und konvertiere zu dB
    let amplitudes: Vec<f64> = spectrum.iter().map(|c| c.norm()).collect();

    // Berechne die Frequenzen in Hertz
    let frequency_resolution = fs / signal.len() as f64;
    let frequencies: Vec<f64> = (0..signal.len() / 2).map(|i| i as f64 * frequency_resolution).collect();

    // Spektrum visualisieren mit Hertz auf der X-Achse
    let root_area = BitMapBackend::new(filename, (1920, 1080)).into_drawing_area();
    root_area.fill(&WHITE).unwrap();
    let max_amplitude = amplitudes.iter().cloned().fold(0f64, f64::max);

    let mut chart = ChartBuilder::on(&root_area)
        .caption("AM Signal Spectrum", ("sans-serif", 20))
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(lower_freq_limit..upper_limit, 0.0..max_amplitude)
        .unwrap();

    chart.configure_mesh().x_desc("Frequency (Hz)").y_desc("Amplitude").draw().unwrap();

    chart
        .draw_series(LineSeries::new(
            frequencies.iter().zip(amplitudes.iter()).map(|(&f, &amp)| (f, amp)),
            &BLUE,
        ))
        .unwrap();
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let fs = 10000000.0;  // Abtastrate
    let duration = 1.0; // Sekunden
    let freq_carrier = 100000.0; // Trägerfrequenz in Hz
    let freq_carrier2 = 400000.0; // Trägerfrequenz in Hz
    let freq_mod = 5000.0; // Modulationsfrequenz in Hz
    let freq_mod2 = 100.0;
    let mod_index = 0.5; // Modulationsindex

    let am_signal = generate_am_signal(freq_carrier, freq_mod, mod_index, fs, duration);
    let am_signal2 = generate_am_signal(freq_carrier2, freq_mod2, mod_index, fs, duration);
    let two_ams: Vec<_> = am_signal.iter().zip(am_signal2.iter()).map(|(&a, &b)| a + b).collect();

    assert_eq!(am_signal.len(), two_ams.len());
    assert_eq!(am_signal.len(), am_signal2.len());

    for i in 0..two_ams.len() {
        assert_eq!(am_signal[i]+am_signal2[i], two_ams[i]);
    }

    let filtered_signal = filter_signal(&two_ams, fs.hz(), freq_carrier.hz(), 30.0);
    let filtered_signal2 = filter_signal(&two_ams, fs.hz(), freq_carrier2.hz(), 30.0);

    assert_eq!(filtered_signal.len(), filtered_signal2.len());

    let peaks = extract_peaks(&filtered_signal, fs).iter().map(|x| x.to_owned()).collect::<Vec<_>>();
    let peaks2 = extract_peaks(&filtered_signal2, fs).iter().map(|x| x.to_owned()).collect::<Vec<_>>();

    let centered_wave1 = remove_dc_component(&peaks,(peaks.len() as f64/10.0).hz(),10.hz());
    let centered_wave2 = remove_dc_component(&peaks2,(peaks2.len() as f64/10.0).hz(),10.hz());

    do_fft("air-travelling.png",two_ams.len() as f64/duration,&two_ams, freq_carrier2+10000f64, freq_carrier-10000f64);
    do_fft("filtered-1-1.png",filtered_signal.len() as f64/duration,&filtered_signal, freq_carrier2+10000f64, freq_carrier-10000f64);
    do_fft("filtered-2-1.png",filtered_signal2.len() as f64/duration,&filtered_signal2, freq_carrier2+10000f64, freq_carrier-10000f64);

    do_fft("peaks-1.png",peaks.len() as f64/duration,&peaks, 20000.0, -100.0);
    do_fft("peaks-2.png",peaks2.len() as f64/duration,&peaks2, 20000.0, -100.0);

    do_fft("peaks-dc-rem-1.png",centered_wave1.len() as f64/duration,&centered_wave1, 20000.0, -100.0);
    do_fft("peaks-de-rem-2.png",centered_wave2.len() as f64/duration,&centered_wave2, 20000.0, -100.0);

    let root_area = SVGBackend::new("other.svg", (2560, 1440)).into_drawing_area();
    root_area.fill(&WHITE).unwrap();
    let max_amplitude = centered_wave1.iter().cloned().fold(0f64, f64::max);
    let min_amplitude = centered_wave1.iter().cloned().fold(0f64, f64::min);

    let mut chart = ChartBuilder::on(&root_area)
        .caption("AM-Signal domodulated", ("sans-serif", 20))
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(0.02*10f64..0.02*20f64, min_amplitude..max_amplitude)
        .unwrap();

    chart.configure_mesh().x_desc("Zeit").y_desc("Signal").draw().unwrap();

    let things1 = (0..centered_wave1.len()).map(|i| (((i as f64)/centered_wave1.len() as f64)*10.0, centered_wave1[i]));
    let things2 = (0..centered_wave2.len()).map(|i| (((i as f64)/centered_wave2.len() as f64)*10.0, centered_wave2[i]));
    chart
        .draw_series(LineSeries::new(
            things1,
            &BLUE
        ))
        .unwrap();
    chart
        .draw_series(LineSeries::new(
            things2,
            &RED
        ))
        .unwrap();



    Ok(())
}
