use tauvslagrange::{
    prover::Prover,
    serialize::SerializedSRS,
    srs::generate_srs,
    utils::{random_fr, random_poly, to_lagrange_basis},
};

#[macro_export]
macro_rules! time_it {
    ($label:expr, $block:expr) => {{
        println!("{} ...", $label);
        let start = std::time::Instant::now();
        let result = $block;
        let elapsed = start.elapsed();
        println!("{} - Elapsed: {:?}", $label, elapsed);
        result
    }};
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("*******************************");
    println!("*                             *");
    println!("*  Powers of Tau vs Lagrange  *");
    println!("*                             *");
    println!("*******************************");

    let mut rl = rustyline::DefaultEditor::new()?;
    let n = 2_usize.pow(17);

    loop {
        // Display options to the user
        println!("\n\nOptions:");
        println!("1. Run commitment with pre-generated SRS");
        println!("2. Generate new SRS");
        println!("3. Exit");

        let readline = rl.readline("> ");
        match readline {
            Ok(line) => match line.trim() {
                "1" => {
                    println!("\n\n------------ Setup ------------");
                    let tau_srs =
                        time_it!("Loading powers of tau", SerializedSRS::load("srs.json")?);
                    let lagrange_srs = time_it!(
                        "Loading powers of tau in Lagrange basis",
                        SerializedSRS::load("lagrange_srs.json")?
                    );

                    // generate a random polynomial of degree n-1
                    let poly = time_it!("Polynomial Generation", { random_poly(n - 1) });
                    let prover = Prover::new(poly)?;

                    println!("\n\n------------ Prover ------------");
                    let witness = time_it!("Witness Generation", random_poly(n - 1));
                    let commitment1 = time_it!("Commitment Calculation (Powers of Tau)", {
                        prover.commit_polynomial(&witness, tau_srs.to_ec_points().as_slice())
                    })?;

                    let commitment2 = time_it!("Commitment Calculation (Lagrange)", {
                        prover.commit_lagrange(&witness, &lagrange_srs.to_ec_points().as_slice())
                    })?;

                    println!("\n\n------------ Result ------------");
                    println!(
                        "Commitment[t] G1: ({},{})",
                        commitment1.to_affine().x(),
                        commitment1.to_affine().y()
                    );
                    println!(
                        "Commitment[l] G1: ({},{})",
                        commitment2.to_affine().x(),
                        commitment2.to_affine().y()
                    );
                }
                "2" => {
                    println!("\n\n------------ Setup ------------");
                    let srs = time_it!("SRS Generation", { generate_srs(2 * n, random_fr()) });

                    let lagrange_srs =
                        time_it!("Lagrange SRS Generation", { to_lagrange_basis(&srs)? });

                    SerializedSRS::from(srs).dump("srs.json")?;
                    SerializedSRS::from(lagrange_srs).dump("lagrange_srs.json")?;
                }
                "3" => {
                    println!("Bye!");
                    break;
                }
                _ => {
                    println!("Invalid option. Try again.");
                    continue;
                }
            },
            Err(rustyline::error::ReadlineError::Interrupted) => {
                println!("Input interrupted by user. Exiting.");
                break;
            }
            Err(_) => {
                println!("Bye!");
                break;
            }
        }
    }

    Ok(())
}
