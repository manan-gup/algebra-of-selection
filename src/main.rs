use symbolica::{atom::Atom, atom::AtomCore, parse, printer::PrintOptions, symbol, LicenseManager};

fn main() {
    LicenseManager::set_license_key("").unwrap();

    // Symbols
    let (Q1, Q3, Z1, Z2, Z3) = symbol!("Q1", "Q3", "Z1", "Z2", "Z3");
    let Q2 = (Q1 + Q3) * (-1) + 1;
    // println!("{}", &Q2.printer(PrintOptions::latex()));

    // Generation 1

    // Pre-selection population state

    // Mean phenotype
    let EZ = Atom::var(Q1) * Atom::var(Z1) + &Q2 * Atom::var(Z2) + Atom::var(Q3) * Atom::var(Z3);

    // Mean allele dosage (a1 = 2, a2 = 1, a3 = 0)
    let Ea = Atom::var(Q1) * 2 + &Q2;
    let Eaa = Atom::var(Q1) * 4 + &Q2;
    let EaZ = Atom::var(Q1) * Atom::var(Z1) * 2 + Atom::var(Z2) * &Q2;
    // println!("EaZ = {}", &EaZ.printer(PrintOptions::latex()));

    // Z vs a regression
    let Beta_Za = ((&EaZ - &Ea * &EZ) / (&Eaa - &Ea * &Ea)).expand();
    let K_Za = (&EZ - &Beta_Za * &Ea).expand();
    // println!(
    //     "Beta_Za = {}\nK_Za = {}",
    //     &Beta_Za.printer(PrintOptions::sympy()),
    //     &K_Za.printer(PrintOptions::sympy())
    // );

    // Additive genetic values Ai
    let A1 = &Beta_Za * &K_Za * 2;
    let A2 = &Beta_Za * &K_Za;
    // A3 = 0;

    // Mean additive value E(A)
    let EA = Atom::var(Q1) * &A1 + &Q2 * &A2;

    // Selection

    // Mean fitness
    // Use EZ for EW
    let EWA = Atom::var(Q1) * Atom::var(Z1) * &A1 + &Q2 * Atom::var(Z2) * &A2;

    // Post selection frequency
    let Q11 = Atom::var(Q1) * (Atom::var(Z1) / &EZ);
    let Q21 = &Q2 * (Atom::var(Z2) / &EZ);
    let Q31 = Atom::var(Q3) * (Atom::var(Z3) / &EZ);

    // Fisherian selection response (Cov(g,w)/EW)
    let SR_F = (&EWA - &EZ * &EA) / &EZ;
    // println!("SR_F = {}", &SR_F.together().printer(PrintOptions::latex()));

    // Mating

    // Normalized mating frequencies (fij)
    let f11 = &Q11 * &Q11;
    let f12 = &Q11 * &Q21;
    let f13 = &Q11 * &Q31;
    let f21 = &Q21 * &Q11;
    let f22 = &Q21 * &Q21;
    let f23 = &Q21 * &Q31;
    let f31 = &Q31 * &Q11;
    let f32 = &Q31 * &Q21;
    let f33 = &Q31 * &Q31;

    // Generation 2

    // Post mating / next generation frequencies
    let Q111 = &f11 + &f12 / 2 + &f21 / 2 + &f22 / 4;
    let Q211 = (&f12 + &f13 + &f21 + &f22 + &f23 + &f32) / 2 + &f31;
    let Q311 = &f22 / 4 + (&f23 + &f32) / 2 + &f33;

    // Generation 2 phenotypic mean
    let Ez = Atom::var(Z1) * &Q111 + Atom::var(Z2) * &Q211 + Atom::var(Z1) * &Q311;

    // Total selection response
    let SR_T = &Ez - &EZ;

    // Discrepancy: TSR - FSR
    let Ds_1 = &SR_T - &SR_F;
    // Full Ds expression
    let Ds = Ds_1.expand();
    println!(
        "Raw Ds expression:\n{}",
        &Ds.printer(PrintOptions::sympy().hide_namespace("satya_symbolica"))
    );

    // Full foctorization of Ds
    let Ds_factored = &Ds.factor();
    println!(
        "Full foctorized Ds:\n{}",
        &Ds_factored.printer(PrintOptions::sympy().hide_namespace("satya_symbolica"))
    );
}
